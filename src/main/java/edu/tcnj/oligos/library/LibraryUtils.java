package edu.tcnj.oligos.library;

import com.google.common.collect.Maps;
import edu.tcnj.oligos.data.Codon;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class LibraryUtils {
    /**
     * Determine whether any Sequence in the list contains a restriction enzyme site
     *
     * @param possiblePermutations a list of Sequences to check
     * @param restrictions         a list of BaseSequences to be avoided
     * @return true iff any one of the restrictions appears in any of the possiblePermutations
     */
    static boolean containsRestrictionEnzyme(List<Sequence> possiblePermutations,
                                                    List<BaseSequence> restrictions) {
        if (restrictions == null || restrictions.isEmpty()) return false;
        for (Sequence possiblePermutation : possiblePermutations) {
            BaseSequence baseSequence = possiblePermutation.asBases();
            for (BaseSequence restriction : restrictions) {
                if (baseSequence.indexOf(restriction) != -1) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Determine whether this sequence contains a restriction enzyme site
     *
     * @param sequence     the sequence to check
     * @param restrictions a list of BaseSequences to be avoided
     * @return true iff the sequence contains any of the restrictions
     */
    static boolean containsRestrictionEnzyme(Sequence sequence, List<BaseSequence> restrictions) {
        if (restrictions == null || restrictions.isEmpty()) return false;
        for (BaseSequence restriction : restrictions) {
            if (sequence.asBases().indexOf(restriction) != -1) return true;
        }
        return false;
    }

    /**
     * Build all possible permutations of oligos in a given range
     *
     * @param range         the range of oligo positions to be permuted (inclusive)
     * @param oligos        the oligos to be permuted
     * @param oligoLength   the length of an oligo
     * @param overlapLength the length of an overlap
     * @return a list of all possible BaseSequences made from these oligos
     */
    public static List<Sequence> buildPermutations(Fragment.Range range, Map<Integer, List<Oligo>> oligos,
                                                   int oligoLength, int overlapLength) {
        return buildPermutationsRecursive(range.getStartPosition(), range, oligos, oligoLength, overlapLength);
    }

    //Recursively build all permutations
    private static List<Sequence> buildPermutationsRecursive(int position, Fragment.Range range,
                                                             Map<Integer, List<Oligo>> oligos, int oligoLength, int overlapLength) {
        checkInterrupt();
        //Base case: for the last position, just return a list of the oligos as BaseSequences
        List<Sequence> finalPermutations = new ArrayList<>();
        if (position == range.getEndPosition()) {
            finalPermutations.addAll(oligos.get(position));
            return finalPermutations;
        }
        //Recurse: find all possible sequences made of all oligos in later positions than this current one
        List<Sequence> partialPermutations = buildPermutationsRecursive(position + 1, range, oligos,
                oligoLength, overlapLength);
        //For every oligo in this position
        for (Oligo oligo : oligos.get(position)) {
            //Get the ending overlap base sequence for this oligo
            Sequence overlap = oligo.subList(oligoLength - overlapLength, oligoLength);
            for (Sequence partialPermutation : partialPermutations) {
                //If any of the permutations start with this overlap, that means that this oligo can connect to it
                if (Sequence.regionsMatch(partialPermutation, 0, overlapLength, overlap, 0, overlapLength)) {
                    //Make the sequence produced by putting this oligo in
                    //front of the permutation, add it to fullSequences
                    Sequence fullSequence = new Sequence(new ArrayList<>(oligo));
                    fullSequence.addAll(partialPermutation.subList(overlapLength, partialPermutation.size()));
                    finalPermutations.add(fullSequence);
                }
            }
        }

        return finalPermutations;
    }

    /**
     * get the best approximate counts to fill the spots with the given relative frequencies
     *
     * @param frequencies the relative frequencies of all codons to be used
     * @param totalCount  the number of spots to be filled
     * @return a map from codon to the number of times it should appear (counts will sum to totalCount)
     */
    static Map<Codon, Integer> findCodonCounts(Map<Codon, Double> frequencies, int totalCount) {
        Map<Codon, Integer> counts = Maps.newHashMap();
        int totalSoFar = 0;
        //Get approximate counts for every codon, rounding down always
        for (Map.Entry<Codon, Double> entry : frequencies.entrySet()) {
            int thisCount = (int) (entry.getValue() * totalCount);
            totalSoFar += thisCount;
            counts.put(entry.getKey(), thisCount);
        }
        //Hit the total needed by continuously adding one to the count which
        //is furthest under its expected (exact, non-integral) count
        while (totalSoFar < totalCount) {
            double maxDelta = Double.NEGATIVE_INFINITY;
            Codon codonToIncrease = Codon.PAD;
            for (Map.Entry<Codon, Double> entry : frequencies.entrySet()) {
                double thisDelta = entry.getValue() * totalCount - counts.get(entry.getKey());
                if (thisDelta > maxDelta) {
                    maxDelta = thisDelta;
                    codonToIncrease = entry.getKey();
                }
            }
            counts.put(codonToIncrease, counts.get(codonToIncrease) + 1);
            totalSoFar++;
        }
        return counts;
    }

    static void checkInterrupt() {
        if (Thread.currentThread().isInterrupted()) {
            Exception e = new InterruptedException("Execution was cancelled.");
            throw new RuntimeException(e.getMessage(), e);
        }
    }
}
