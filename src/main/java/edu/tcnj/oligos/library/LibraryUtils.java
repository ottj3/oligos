package edu.tcnj.oligos.library;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static edu.tcnj.oligos.library.Library.checkInterrupt;

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
            BaseSequence overlap = oligo.subList(oligoLength - overlapLength, oligoLength).asBases();
            for (Sequence partialPermutation : partialPermutations) {
                //If any of the permutations start with this overlap, that means that this oligo can connect to it
                if (partialPermutation.asBases().indexOf(overlap) == 0) {
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
}
