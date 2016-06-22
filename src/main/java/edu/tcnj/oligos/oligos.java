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
        int numDesigns = 3;
        int numRegions = 10;
        int overlapLength = 1;
        for (int i = 0; i < numRegions; i++) {
            ArrayList<RNASequence> list = new ArrayList<>();
            RNASequence oli1 = new RNASequence(Integer.toString(i));
//            oli1.data.add(Codon.CTC);
//            oli1.data.add(Codon.CTT);
//            oli1.data.add(Codon.TTA);
//            oli1.data.add(Codon.CTG);
//            oli1.data.add(Codon.TTG);
            oli1.data.add(Codon.ACC);
            oli1.data.add(Codon.GTC);
            oli1.data.add(Codon.GAA);
            list.add(oli1);
            sequence.add(list);
        }
//        Codon[] codons = {Codon.CTA};
//        Integer[][][] fragments = {{{0, 1}, {2, 3}}};
//        Integer[][][] designs = {{{0, 1}, {0, 2}}};
        String[] names = new String[numDesigns];
        Codon[] codons = {Codon.ACT, Codon.GTG, Codon.GAG};
        Integer[][][] fragments = {{{1, 2}, {7, 9}}, {{0, 0}, {2, 7}}, {{0, 1}, {3, 5}}};//new Integer[numDesigns][2][];
        Integer[][][] designs = {{{0, 2}, {0, 1}}, {{0, 1}, {0, 2}}, {{0, 1, 2}, {0, 3}}};//new Integer[numDesigns][][];
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
        //TODO: fill in other codons, check overlaps

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
        for (List<RNASequence> overlapRegion : overlaps) {
            //For every overlap except the first (first will stay unchanged)
            for (int overlapI = 1; overlapI < overlapRegion.size(); overlapI++) {
                RNASequence thisOverlap = overlapRegion.get(overlapI);
                //For every codon in this overlap, try to change it to make this overlap unique
                for (int codonI = 0; codonI < thisOverlap.data.size(); codonI++) {
                    //If it's already unique, good job. done.
                    if (!matchesAnyPrevious(overlapRegion, overlapI)) break;

                    //Find a codon in thisOverlap.backEdges that can swap with the current codon
                    //TODO: If the backEdge changes don't work, try forwardEdges
                    RNASequence preOligo = thisOverlap.backEdges.get(0);
                    for (int preI = 0; preI < preOligo.data.size() - thisOverlap.data.size(); preI++) {
                        if (preOligo.data.get(preI).getAminoAcid() == thisOverlap.data.get(codonI).getAminoAcid()) {
                            //Makes sure that none of the oligos have a controlled codon in this position
                            boolean isImportant = false;
                            for (RNASequence backOligo : thisOverlap.backEdges) {
                                for (Codon codon : codons) {
                                    if (codon == backOligo.data.get(preI)) {
                                        isImportant = true;
                                    }
                                }
                            }
                            if (isImportant) continue;
                            //Swap the codon from every backEdge with thisOverlap's current codon
                            Codon swapCodon = preOligo.data.get(preI);
                            for (RNASequence backOligo : thisOverlap.backEdges) {
                                backOligo.data.set(preI, thisOverlap.data.get(codonI));
                            }
                            thisOverlap.data.set(codonI, swapCodon);
                            //If the swap made the overlap unique, good. If not, try the next swap for this codon.
                            if (!matchesAnyPrevious(overlapRegion, overlapI)) break;
                        }
                    }
                }
                for (RNASequence backOligo : thisOverlap.backEdges) {
                    int startPos = backOligo.data.size() - thisOverlap.data.size();
                    for (int i = startPos; i < backOligo.data.size(); i++) {
                        backOligo.data.set(i, thisOverlap.data.get(i - startPos));
                    }
                }
                for (RNASequence forwardOligo : thisOverlap.forwardEdges) {
                    for (int i = 0; i < thisOverlap.data.size(); i++) {
                        forwardOligo.data.set(i, thisOverlap.data.get(i));
                    }
                }
            }
        }
        return;
    }

    private static boolean matches(RNASequence overlap1, RNASequence overlap2) {
        for (int codonI = 0; codonI < overlap1.data.size() && codonI < overlap2.data.size(); codonI++) {
            if (overlap1.data.get(codonI) != overlap2.data.get(codonI)) {
                return false;
            }
        }
        return true;
    }

    private static boolean matchesAnyPrevious(List<RNASequence> overlapRegion, int index) {
        for (int prevOverlapI = 0; prevOverlapI < index; prevOverlapI++) {
            if (matches(overlapRegion.get(index), overlapRegion.get(prevOverlapI))) {
                return true;
            }
        }
        return false;
    }
}
