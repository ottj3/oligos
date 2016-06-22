package edu.tcnj.oligos;

import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Codon;
import edu.tcnj.oligos.data.RNASequence;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class oligos {
    private static List<List<RNASequence>> sequence = new ArrayList<>();
    private static List<List<RNASequence>> overlaps = new ArrayList<>();

    //{0:{oligoA, oligoB},1:{oligoC, oligoD}}
    public static void main(String[] args) {
        //do stuff, put base sequence into sequence[0..n][0]
        int numDesigns = 1;
        int numRegions = 4;
        int overlapLength = 2;
        for (int i = 0; i < numRegions; i++) {
            ArrayList<RNASequence> list = new ArrayList<>();
            RNASequence oli1 = new RNASequence(Integer.toString(i));
            oli1.data.add(Codon.CTC);
            oli1.data.add(Codon.CTT);
            oli1.data.add(Codon.TTA);
            oli1.data.add(Codon.CTG);
            oli1.data.add(Codon.TTG);
//            oli1.data.add(Codon.ACC);
//            oli1.data.add(Codon.GTC);
//            oli1.data.add(Codon.GAA);
            list.add(oli1);
            sequence.add(list);
        }
        Codon[] codons = {Codon.CTA};
        Integer[][][] fragments = {{{0, 1}, {2, 3}}};
        Integer[][][] designs = {{{0, 1}, {0, 2}}};
        String[] names = new String[numDesigns];
//        Codon[] codons = {Codon.ACT, Codon.GTG, Codon.GAG};
//        Integer[][][] fragments = {{{1, 2}, {7, 9}}, {{0, 0}, {2, 7}}, {{0, 1}, {3, 5}}};//new Integer[numDesigns][2][];
//        Integer[][][] designs = {{{0, 2}, {0, 1}}, {{0, 1}, {0, 2}}, {{0, 1, 2}, {0, 3}}};//new Integer[numDesigns][][];
        Integer[][][] designCopy = new Integer[numDesigns][][];
        for (int i = 0; i < designs.length; i++) {
            designCopy[i] = new Integer[designs[i].length][];
            for (int j = 0; j < designs[i].length; j++) {
                designCopy[i][j] = new Integer[designs[i][j].length];
                System.arraycopy(designs[i][j], 0, designCopy[i][j], 0, designs[i][j].length);
            }
        }

        //For every amino acid in the library
        for (int acidI = 0; acidI < names.length; acidI++) {
            //For every fragment used by this acid
            for (int fragmentI = 0; fragmentI < fragments[acidI].length; fragmentI++) {
                //For every position in this fragment
                for (int positionI = fragments[acidI][fragmentI][0]; positionI <= fragments[acidI][fragmentI][1]; positionI++) {
                    int size = sequence.get(positionI).size();
                    //For every delta needed in this fragment
                    for (int deltaI = 0; deltaI < designs[acidI][fragmentI].length; deltaI++) {
                        //Make more copies of the current sequences if necessary
                        while (sequence.get(positionI).size() < (deltaI + 1) * size) {
                            for (int copyI = 0; copyI < size; copyI++) {
                                sequence.get(positionI).add(sequence.get(positionI).get(copyI).clone());
                            }
                        }
                        int numLeft = 0;
                        //For every copy of the current sequence
                        for (int copyI = (deltaI * size); copyI < (deltaI + 1) * size; copyI++) {
                            //Fill it with as many of the desired codon as possible
                            numLeft = sequence.get(positionI).get(copyI)
                                    .fill(codons[acidI], designs[acidI][fragmentI][deltaI]);
                            //Use designCopy to mark this oligo with the corresponding delta (used for overlaps)
                            sequence.get(positionI).get(copyI)
                                    .deltas.put(codons[acidI].getAminoAcid(), designCopy[acidI][fragmentI][deltaI]);
                        }
                        designs[acidI][fragmentI][deltaI] = numLeft;
                    }
                }
            }
        }
        //TODO: fill in other codons

        for (int positionI = 0; positionI < numRegions - 1; positionI++) {
            //TODO: fill arrayList with proper template codons for the overlap
            List<RNASequence> list = new ArrayList<>();
            RNASequence thisOverlap = new RNASequence("");
            thisOverlap.data = sequence.get(positionI).get(0).data.subList(sequence.get(positionI).get(0).data.size() - overlapLength, sequence.get(positionI).get(0).data.size());
            list.add(thisOverlap);
            overlaps.add(list);
        }

        //For every amino acid in the library
        for (int acidI = 0; acidI < names.length; acidI++) {
            //For every fragment used by this acid
            for (int fragmentI = 0; fragmentI < fragments[acidI].length; fragmentI++) {
                //For every position in this fragment
                for (int positionI = fragments[acidI][fragmentI][0]; positionI < fragments[acidI][fragmentI][1]; positionI++) {
                    int size = overlaps.get(positionI).size();
                    //For every delta needed in this fragment
                    for (int deltaI = 0; deltaI < designs[acidI][fragmentI].length; deltaI++) {
                        //Make more copies of the current sequences if necessary
                        while (overlaps.get(positionI).size() < (deltaI + 1) * size) {
                            for (int copyI = 0; copyI < size; copyI++) {
                                overlaps.get(positionI).add(overlaps.get(positionI).get(copyI).clone());
                            }
                        }
                        //For every copy of the current sequence
                        for (int copyI = (deltaI * size); copyI < (deltaI + 1) * size; copyI++) {
                            //Use designCopy to mark this overlap with the corresponding delta
                            overlaps.get(positionI).get(copyI)
                                    .deltas.put(codons[acidI].getAminoAcid(), designCopy[acidI][fragmentI][deltaI]);
                        }
                    }
                }
            }
        }

        for (int overlapI = 0; overlapI < numRegions - 1; overlapI++) {
            for (RNASequence overlap : overlaps.get(overlapI)) {
                for (RNASequence oligo : sequence.get(overlapI)) {
                    boolean edgeExists = true;
                    for (Map.Entry<AminoAcid, Integer> entry : overlap.deltas.entrySet()) {
                        if (oligo.deltas.get(entry.getKey()).compareTo(entry.getValue()) != 0) {
                            edgeExists = false;
                            break;
                        }
                    }
                    if (edgeExists) {
                        overlap.backEdges.add(oligo);
                        oligo.forwardEdges.add(overlap);
                    }
                }
                for (RNASequence oligo : sequence.get(overlapI + 1)) {
                    boolean edgeExists = true;
                    for (Map.Entry<AminoAcid, Integer> entry : overlap.deltas.entrySet()) {
                        if (oligo.deltas.get(entry.getKey()).compareTo(entry.getValue()) != 0) {
                            edgeExists = false;
                            break;
                        }
                    }
                    if (edgeExists) {
                        overlap.forwardEdges.add(oligo);
                        oligo.backEdges.add(overlap);
                    }
                }
            }
        }

        //For every region of overlaps, make sure that every overlap is unique
        for (int overlapRegionI = 0; overlapRegionI < overlaps.size(); overlapRegionI++) {
                List<RNASequence> overlapRegion = overlaps.get(overlapRegionI);

            //Get all possible substitutions for every codon in the overlap,
            //permute those possibilities to generate unique overlaps
            //TODO: try "do nothing" as possible permutation?
            List<List<Integer>> possibleSubstitutions = new ArrayList<>(overlapRegion.get(0).data.size());
            RNASequence preOligo = overlapRegion.get(0).backEdges.get(0);
            for (int overlapPos = 0; overlapPos < overlapRegion.get(0).data.size(); overlapPos++) {
                //For every position in the overlaps, get indices of all unique and
                //non-controlled codons that correspond to the given amino acid. If
                //any acid at the given index is controlled, that index is ignored.
                possibleSubstitutions.add(new ArrayList<Integer>());
                //TODO: search globally to find possible swaps?
                for (int oligoPos = 0; oligoPos < preOligo.data.size() - overlapRegion.get(0).data.size(); oligoPos++) {
                    //TODO: use master protein sequence to avoid looking at specific oligo/overlap
                    if (preOligo.data.get(oligoPos).getAminoAcid() ==
                            overlapRegion.get(0).data.get(overlapPos).getAminoAcid()) {
                        //For potential substitutions, make sure that no controlled codon shares this index
                        for (RNASequence overlap : overlapRegion) {
                            boolean isImportant = false;
                            for (RNASequence backOligo : overlap.backEdges) {
                                for (Codon codon : codons) {
                                    if (codon == backOligo.data.get(oligoPos)) {
                                        isImportant = true;
                                    }
                                }
                            }
                            if (!isImportant) {
                                //make sure that this codon is not redundantly added to the list
                                boolean isRepeated = false;
                                for (Integer pos : possibleSubstitutions.get(overlapPos)) {
                                    if (preOligo.data.get(pos) == preOligo.data.get(oligoPos)) {
                                        isRepeated = true;
                                        break;
                                    }
                                }
//                                if (preOligo.data.get(oligoPos) == overlap.data.get(overlapPos)) {
//                                    isRepeated = true;
//                                }
                                if (!isRepeated) {
                                    possibleSubstitutions.get(overlapPos).add(oligoPos);
                                }
                            }
                        }
                    }
                }
            }
            List<Integer> substitutionIndices = new ArrayList<>();
            for (int i = 0; i < possibleSubstitutions.size(); i++) {
                substitutionIndices.add(0);
            }
            for (int i = 0; i < overlapRegion.size(); i++) {
//                for (Integer substitutionIndex : substitutionIndices) {
//                    System.out.print(substitutionIndex + " ");
//                }
//                System.out.println("");
                RNASequence overlap = overlapRegion.get(i);
                for (int overlapPos = 0; overlapPos < overlap.data.size(); overlapPos++) {
                    int substitutionIndex = substitutionIndices.get(overlapPos);
                    int preI = possibleSubstitutions.get(overlapPos).get(substitutionIndex);
                    Codon swapCodon = preOligo.data.get(preI);
                    for (RNASequence backOligo : overlap.backEdges) {
                        backOligo.data.set(preI, overlap.data.get(overlapPos));
                    }
                    overlap.data.set(overlapPos, swapCodon);
                }
                int overlapPos = 0;
                substitutionIndices.set(overlapPos, substitutionIndices.get(overlapPos) + 1);
                while (substitutionIndices.get(overlapPos) >= possibleSubstitutions.get(overlapPos).size()) {
                    substitutionIndices.set(overlapPos, 0);
                    overlapPos++;
                    if (overlapPos >= substitutionIndices.size()) {
                        if (i + 1 < overlapRegion.size()) {
                            //TODO: If the backEdge changes don't work, try forwardEdges?
                            System.out.println("Ran out of options to make unique overlap regions.");
                            return;
                        } else break;
                    }
                    substitutionIndices.set(overlapPos, substitutionIndices.get(overlapPos) + 1);
                }
                //Safeguard: if two codons of the same acid are both varied, then the result may not be unique
                //in some cases. If this permutation matches a previous one, redo it with the new indices.
                if (matchesAnyPrevious(overlaps, overlapRegionI, i)) {
                    i--;
//                    System.out.println("Matches previous permutation, trying a different one instead.");
                }
            }
//            System.out.println("");
        }
        return;
    }
    private static final int NUM_DIFFERENCES_NEEDED = 1;
    private static boolean matches(RNASequence overlap1, RNASequence overlap2) {
        int numDifferences = 0;
        for (int codonI = 0; codonI < overlap1.data.size() && codonI < overlap2.data.size(); codonI++) {
            if (overlap1.data.get(codonI) != overlap2.data.get(codonI)) {
                numDifferences++;

            }
        }
        return (numDifferences < NUM_DIFFERENCES_NEEDED);
    }

    private static boolean matchesAnyPrevious(List<List<RNASequence>> overlaps, int index1, int index2) {
        //look at every overlap in every earlier overlap region
        for (int overlapRegionI = 0; overlapRegionI < index1; overlapRegionI++) {
            for (int prevOverlapI = 0; prevOverlapI < overlaps.get(overlapRegionI).size(); prevOverlapI++) {
                if (matches(overlaps.get(index1).get(index2), overlaps.get(overlapRegionI).get(prevOverlapI))) {
                    return true;
                }
            }
        }
        //for current overlap region, only look up to the current index
        for (int prevOverlapI = 0; prevOverlapI < index2; prevOverlapI++) {
            if (matches(overlaps.get(index1).get(index2), overlaps.get(index1).get(prevOverlapI))) {
                return true;
            }
        }
        return false;
    }
}
