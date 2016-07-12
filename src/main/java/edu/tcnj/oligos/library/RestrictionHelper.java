package edu.tcnj.oligos.library;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class RestrictionHelper {
    public static boolean containsRestrictionEnzyme(List<BaseSequence> possiblePermutations, List<BaseSequence> restrictions) {
        if (restrictions == null || restrictions.isEmpty()) return false;
        for (BaseSequence possiblePermutation : possiblePermutations) {
            for (BaseSequence restriction : restrictions) {
                if (possiblePermutation.contains(restriction) != -1) {
                    System.out.println("Found a restriction enzyme, fixing it.");
                    return true;
                }
            }
        }
        return false;
    }

    public static boolean containsRestrictionEnzyme(Sequence sequence, List<BaseSequence> restrictions) {
        if (restrictions == null || restrictions.isEmpty()) return false;
        for (BaseSequence restriction : restrictions) {
            if (sequence.asBases().contains(restriction) != -1) return true;
        }
        return false;
    }

    public static List<BaseSequence> buildPermutations(Fragment.Range range, Map<Integer, List<Oligo>> oligos, int oligoLength, int overlapLength) {
        return buildPermutationsRecursive(range.getStartPosition(), range, oligos, oligoLength, overlapLength);
    }

    private static List<BaseSequence> buildPermutationsRecursive(int position, Fragment.Range range, Map<Integer, List<Oligo>> oligos, int oligoLength, int overlapLength) {
        List<BaseSequence> finalPermutations = new ArrayList<>();
        if (position == range.getEndPosition()) {
            for (Oligo oligo : oligos.get(position)) {
                finalPermutations.add(oligo.asBases());
            }
        } else {
            List<BaseSequence> partialPermutations = buildPermutationsRecursive(position + 1, range, oligos, oligoLength, overlapLength);
            for (Oligo oligo : oligos.get(position)) {
                //TODO check off-by-one potential in overlap and beginningOverlap
                BaseSequence overlap = oligo.subList(oligoLength - overlapLength, oligoLength).asBases();
                for (BaseSequence partialPermutation : partialPermutations) {
//                    BaseSequence beginningOverlap = partialPermutation.subList(0, (overlapLength * 3));
                    if (partialPermutation.contains(overlap) == 0) {
                        BaseSequence fullSequence = oligo.asBases();
                        fullSequence.addAll(partialPermutation.subList(overlapLength * 3, partialPermutation.size()));
                        finalPermutations.add(fullSequence);
                    }
                }
            }
        }
        return finalPermutations;
    }
}
