package edu.tcnj.oligos.ext;

import com.google.common.base.Joiner;
import com.google.common.collect.Collections2;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.library.Design;
import edu.tcnj.oligos.library.Fragment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Use Ryan's python script to figure out the optimal library, and then builds that library
 */
public class PythonHandler {

    private final String scriptName;
    private final String protein;
    private final String acidsOfInterest;
    private final int oligoLength;
    private final int overlapLength;
    private final double[] minPercentages;
    private final double[] maxPercentages;
    private final int[] numFreqLevels;

    public PythonHandler(String scriptName, String protein, String acidsOfInterest, int oligoLength, int overlapLength,
                         double[] minPercentages, double[] maxPercentages, int[] numFreqLevels) {
        this.scriptName = scriptName;
        this.protein = protein;
        this.acidsOfInterest = acidsOfInterest;
        this.oligoLength = oligoLength;
        this.overlapLength = overlapLength;
        this.minPercentages = minPercentages;
        this.maxPercentages = maxPercentages;
        this.numFreqLevels = numFreqLevels;
    }

    /**
     * Run the script itself and interpret the results
     * @return A map of AminoAcids to the best designs for them
     */
    @SuppressWarnings("unchecked")
    public Map<AminoAcid, Design> run() {
        Object res;
        try {
            res = runScript(protein, acidsOfInterest, oligoLength, overlapLength,
                    minPercentages, maxPercentages, numFreqLevels);
        } catch (IOException e) {
            throw new IllegalStateException("Couldn't run python design script.", e);
        }
        if (!(res instanceof List)) {
            throw new IllegalStateException();
        }
        Map<AminoAcid, Map<Fragment.Range, Integer>> map = Maps.newLinkedHashMap();
        List<Object> resList = ((List<Object>) res);
        for (Object o : resList) {
            List<Object> result = (List<Object>) o;
            //Result: <string (acid name), list of list of ints (ranges), list of ints (number of occurrences)>
            String acid = String.valueOf(result.get(0));
            List<Object> ranges = ((List<Object>) result.get(1));
            List<Integer> occurrencesList = ((List<Integer>) result.get(2));
            //From the list of list of ints, build Fragment.Range objects
            //and map them to the number of occurrences in the interval
            Map<Fragment.Range, Integer> acidMap = Maps.newHashMap();
            for (int i = 0; i < ranges.size(); i++) {
                List<Integer> rangeList = ((List<Integer>) ranges.get(i));
                Fragment.Range range = new Fragment.Range(rangeList.get(0), rangeList.get(1));
                int numOccurrences = occurrencesList.get(i);
                acidMap.put(range, numOccurrences);
            }
            //Build a map of amino acids to ranges to number of occurrences
            map.put(AminoAcid.getAcidForSymbol(acid), acidMap);
        }
        int i = 0;
        Map<AminoAcid, Design> designs = Maps.newLinkedHashMap();
        //Turn the map of acids->ranges->number of occurrences into a map of acids->designs
        for (Map.Entry<AminoAcid, Map<Fragment.Range, Integer>> entry : map.entrySet()) {
            int numpts = numFreqLevels[i];
            double max = maxPercentages[i];
            double min = minPercentages[i];
            int numGlobalOccurrences = protein.length() - protein.replace(entry.getKey().getCh(), "").length();
            //Figure out the number of oligos needed for one delta step
            double delta = (max - min) / (numpts - 1) * numGlobalOccurrences;
            //For the given acid, figure out the design/deltas that optimally fit the data
            designs.put(entry.getKey(), calculateDesign(numpts, delta, entry.getValue()));
            i++;
        }
        return designs;
    }

    private static Object runScript(String protein, String acids, int segmentLength, int overlapLength, double[] mins, double[] maxs, int[] numLevels) throws IOException {
        ProcessBuilder pb = new ProcessBuilder("python", "design.py");
        pb.redirectErrorStream(true);
        // set up process and stdin stream
        Process proc = pb.start();
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(proc.getOutputStream()));

        // write all the args to stdin for python
        bw.write(protein + "\n");
        bw.write(acids + "\n"); //Joiner.on(" ").join(Chars.asList(acids.toCharArray())) + "\n");
        bw.write(segmentLength + "\n");
        bw.write(overlapLength + "\n");
        bw.write(Joiner.on(" ").join(Doubles.asList(mins)) + "\n");
        bw.write(Joiner.on(" ").join(Doubles.asList(maxs)) + "\n");
        bw.write(Joiner.on(" ").join(Ints.asList(numLevels)) + "\n");
        bw.close();

        BufferedReader br = new BufferedReader(new InputStreamReader(proc.getInputStream()));
        String line;
        List<String> rawOut = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            rawOut.add(line);
        }
        try {
            proc.waitFor();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        br.close();

        // parse design script output into all the lists
        // this is all extremely fragile because sharing data between processes through stdin and stdout sucks
        // make sure you have the right design.py or everything will break
        int numOfCodons = Integer.valueOf(rawOut.get(0));
        int lNum = 1;
        List<Object> resList = new LinkedList<>();
        for (int i = 0; i < numOfCodons; i++) {
            List<Object> result = new LinkedList<>();
            String codonCode = rawOut.get(lNum);
            result.add(codonCode);
            lNum++;
            int numOfIntervals = Integer.valueOf(rawOut.get(lNum));
            lNum++;
            List<List<Integer>> ranges = new LinkedList<>();
            for (int j = 0; j < numOfIntervals; j++) {
                int start = Integer.valueOf(rawOut.get(lNum));
                lNum++;
                int end = Integer.valueOf(rawOut.get(lNum));
                lNum++;
                List<Integer> range = new LinkedList<>();
                range.add(start);
                range.add(end);
                ranges.add(range);
            }
            result.add(ranges);
            List<Integer> occurenceList = new LinkedList<>();
            for (int j = 0; j < numOfIntervals; j++) {
                int occurences = Integer.valueOf(rawOut.get(lNum));
                occurenceList.add(occurences);
                lNum++;
            }
            result.add(occurenceList);
            resList.add(result);
        }

        return resList;
    }

    //Helper method used in calculateDesign;
    //designs are calculated based on the prime factors of the number of regions
    private static List<Integer> primeFactors(int number) {
        int n = number;
        List<Integer> factors = new ArrayList<>();
        for (int i = 2; i <= n; i++) {
            while (n % i == 0) {
                factors.add(i);
                n /= i;
            }
        }
        return factors;
    }

    //Comparator used for calculateDesign; the map needs to be sorted in ascending order
    //according to the number of occurrences (as the delta ranges are made from smallest to largest)
    private static Comparator<Map.Entry<Fragment.Range, Integer>> rangeComparator = new Comparator<Map.Entry<Fragment.Range, Integer>>() {
        @Override
        public int compare(Map.Entry<Fragment.Range, Integer> r1, Map.Entry<Fragment.Range, Integer> r2) {
            return r1.getValue().compareTo(r2.getValue());
        }
    };

    private static Design calculateDesign(int numpts, double delta, Map<Fragment.Range, Integer> numInRegion) {
        //List of ranges, to be ordered by numOccurrences
        List<Map.Entry<Fragment.Range, Integer>> ranges = Lists.newArrayList(numInRegion.entrySet());

        //When the delta lists are made, the first delta list has a smaller max than the second delta list, etc.
        //So, the ranges should be sorted according to the number of occurrences in order to make sure the delta
        //lists correspond to  the proper ranges.
        Collections.sort(ranges, rangeComparator);

        List<Integer> factorList = primeFactors(numpts);
        double bestOver = Double.POSITIVE_INFINITY;
        Map<Fragment.Range, List<Integer>> bestDesign = Maps.newHashMap();
        //Try every possible design (every distinct permutation of the prime factors)
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

}