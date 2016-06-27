package edu.tcnj.oligos.ext;

import com.google.common.collect.Collections2;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.library.Design;
import edu.tcnj.oligos.library.Fragment;
import jep.Jep;
import jep.JepException;

import java.util.*;

public class PythonHandler {

    private final String protein;
    private final String acidsOfInterest;
    private final int oligoLength;
    private final int overlapLength;
    private final double[] minPercentages;
    private final double[] maxPercentages;
    private final int[] numFreqLevels;

    public PythonHandler(String protein, String acidsOfInterest, int oligoLength, int overlapLength,
                         double[] minPercentages, double[] maxPercentages, int[] numFreqLevels) {
        this.protein = protein;
        this.acidsOfInterest = acidsOfInterest;
        this.oligoLength = oligoLength;
        this.overlapLength = overlapLength;
        this.minPercentages = minPercentages;
        this.maxPercentages = maxPercentages;
        this.numFreqLevels = numFreqLevels;
    }

    @SuppressWarnings("unchecked")
    public Map<AminoAcid, Design> run() {
        Object res = getDesigns(protein, acidsOfInterest, oligoLength, overlapLength,
                minPercentages, maxPercentages, numFreqLevels);
        if (!(res instanceof List)) {
            throw new IllegalStateException();
        }
        Map<AminoAcid, Map<Fragment.Range, Integer>> map = Maps.newLinkedHashMap();
        List<Object> resList = ((List<Object>) res);
        for (Object o : resList) {
            List<Object> result = ((List<Object>) o);
            String acid = String.valueOf(result.get(0));
            List<Object> ranges = ((List<Object>) result.get(1));
            List<Integer> occurrencesList = ((List<Integer>) result.get(2));
            Map<Fragment.Range, Integer> acidMap = Maps.newHashMap();
            for (int i = 0; i < ranges.size(); i++) {
                List<Integer> rangeList = ((List<Integer>) ranges.get(i));
                Fragment.Range range = new Fragment.Range(rangeList.get(0), rangeList.get(1));
                int numOccurrences = occurrencesList.get(i);
                acidMap.put(range, numOccurrences);
            }
            map.put(AminoAcid.getAcidForSymbol(acid),acidMap);
        }
        int i = 0;
        Map<AminoAcid, Design> designs = Maps.newLinkedHashMap();
        for (Map.Entry<AminoAcid, Map<Fragment.Range, Integer>> entry : map.entrySet()) {
            int numpts = numFreqLevels[i];
            double max = maxPercentages[i];
            double min = minPercentages[i];
            int numGlobalOccurences = protein.length() - protein.replace(entry.getKey().getCh(), "").length();
            double delta = (max - min)/(numpts - 1) * numGlobalOccurences;
            designs.put(entry.getKey(), calculateDesign(numpts, delta, entry.getValue()));
            i++;
        }
        return designs;
    }

    private Object getDesigns(String protein, String acids, int segmentLength, int overlapLength, double[] mins, double[] maxs, int[] numLevels) {
        return runScript(protein, acids, segmentLength, overlapLength, mins, maxs, numLevels);
    }

    private Object runScript(String protein, String acids, int segmentLength, int overlapLength, double[] mins, double[] maxs, int[] numLevels) {
        try {
            Jep jep = new Jep(false);
            jep.runScript("design.py");
            Object res = jep.invoke("compute_best_design", protein, acids, segmentLength, overlapLength, mins, maxs, numLevels);
            return res;
        } catch (JepException e) {
            e.printStackTrace();
        }
        return null;
    }

    private static List<Integer> primeFactors(int number) {
        int n = number;
        List<Integer> factors = new ArrayList<Integer>();
        for (int i = 2; i <= n; i++) {
            while (n % i == 0) {
                factors.add(i);
                n /= i;
            }
        }
        return factors;
    }

    private static Comparator<Map.Entry<Fragment.Range, Integer>> rangeComparator = new Comparator<Map.Entry<Fragment.Range, Integer>>() {
        @Override
        public int compare(Map.Entry<Fragment.Range, Integer> r1, Map.Entry<Fragment.Range, Integer> r2) {
            return r1.getValue().compareTo(r2.getValue());
        }
    };

    private static Design calculateDesign(int numpts, double delta, Map<Fragment.Range, Integer> numInRegion) {
        //List of ranges, ordered by numOccurrences
        List<Map.Entry<Fragment.Range, Integer>> ranges = Lists.newArrayList(numInRegion.entrySet());

        //When the delta lists are made, the first delta list has a smaller max than the second delta list, etc.
        //So, the ranges should be sorted according to the number of occurrences in order to make sure the delta
        //lists correspond to  the proper ranges.
        Collections.sort(ranges, rangeComparator);

        List<Integer> factorList = primeFactors(numpts);
        double bestOver = Double.POSITIVE_INFINITY;
        Map<Fragment.Range, List<Integer>> bestDesign = Maps.newHashMap();
        //Try every possible design
        for (List<Integer> factors : Collections2.orderedPermutations(factorList)) {
            Map<Fragment.Range, List<Integer>> levelsForRange = new HashMap<>();
            int stepSize = 1;
            double amtOver = 0;
            //For every factor: compute the corresponding delta levels
            for (int i = 0; i < factors.size(); i++) {
                int factor = factors.get(i);

                Fragment.Range range = ranges.get(i).getKey();
                Integer numOccurrences = ranges.get(i).getValue();

                List<Integer> deltas = new ArrayList<>();
                //Compute the delta levels
                for (int j = 0; j < factor * stepSize; j += stepSize) {
                    deltas.add(j);
                }
                //See how far the exact max is over the max available (for comparing designs; the farther the
                //max available is under the exact max, the worse the design can approximate that level).

                //deltas.get(deltas.size() - 1) is the max number of deltas needed in this fragment.
                //That times delta gives the max number of occurrences of the acid needed in this fragment.
                double maxNeeded = delta * deltas.get(deltas.size() - 1);
                double diff = maxNeeded - numOccurrences;
                //If there aren't enough occurrences in the region (diff > 0), track how far off the num occurrences is.
                if (diff > 0) {
                    amtOver += diff;
                }
                //Increase the stepSize according to the methods described in the paper
                stepSize *= factor;
                //Add the deltas for this factor to its corresponding range
                levelsForRange.put(range, deltas);
            }
            //Track the best design found so far
            if (amtOver < bestOver) {
                bestOver = amtOver;
                bestDesign = levelsForRange;
            }
        }

        //Given the best design, convert that into the deltasForRange needed for a Design object
        Map<Fragment.Range, List<Integer>> deltasForRange = new HashMap<>();
        for (Map.Entry<Fragment.Range, List<Integer>> entry : bestDesign.entrySet()) {
            //Get the max number possible to cap the approximations
            int max = numInRegion.get(entry.getKey());
            List<Integer> theseCounts = new ArrayList<>();
            //for every delta level, get the best whole-number approximate that is <= max
            for (int level : entry.getValue()) {
                double exactNeeded = level * delta;
                int approximateNeeded = (int) (exactNeeded + 0.5);
                if (approximateNeeded > max) {
                    approximateNeeded = max;
                }
                theseCounts.add(approximateNeeded);
            }
            deltasForRange.put(entry.getKey(), theseCounts);
        }
        return new Design(deltasForRange);
    }

    // test
    public static void main(String[] args) {
        String P = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGVQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK";
        String AOI = "TEV";
        int SEGSIZE = 30;
        int OVERLAPSIZE = 6;
        double[] MINS = {.1, .1, .1};
        double[] MAXS = {.85, .85, .85};
        int[] NUMPTS = {4, 6, 4};
        Map<AminoAcid, Design> res = (new PythonHandler(P, AOI, SEGSIZE, OVERLAPSIZE, MINS, MAXS, NUMPTS)).run();
    }

    private Object getExampleDesigns() {
        String P = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGVQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK";
        String AOI = "TEV";
        int SEGSIZE = 30;
        int OVERLAPSIZE = 6;
        double[] MINS = {.1, .1, .1};
        double[] MAXS = {.85, .85, .85};
        int[] NUMPTS = {4, 4, 6};
        return runScript(P, AOI, SEGSIZE, OVERLAPSIZE, MINS, MAXS, NUMPTS);
    }
}