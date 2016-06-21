package edu.tcnj.oligos;

import edu.tcnj.oligos.data.Codon;
import edu.tcnj.oligos.data.RNASequence;

import java.util.ArrayList;
import java.util.List;

public class oligos {
    private static List<ArrayList<RNASequence>> sequence = new ArrayList<>();
    //{0:{oligoA, oligoB},1:{oligoC, oligoD}}
    public static void main(String[] args) {
        //do stuff, put base sequence into sequence[0..n][0]
        int numDesigns = 3;
        int numRegions = 10;
        for (int i = 0; i < numRegions; i++) {
            ArrayList<RNASequence> list = new ArrayList<>();
            RNASequence oli1 = new RNASequence(Integer.toString(i));
            oli1.data.add(Codon.AAA);
            oli1.data.add(Codon.CAC);
            list.add(oli1);
            sequence.add(list);
        }
        Codon[] codons = {Codon.AAA, Codon.CAT, Codon.CCC};
        String[] names = new String[numDesigns];
        Integer[][][] fragments = {{{1,2},{7,9}}, {{0,0},{2,7}}, {{0,1},{3,5}}};//new Integer[numDesigns][2][];
        Integer[][][] designs =   {{{0,2},{0,1}}, {{0,1},{0,2}}, {{0,1,2},{0,3}}};//new Integer[numDesigns][][];
        Integer[][][] designCopy = new Integer[numDesigns][][];
        for (int i = 0; i < designs.length; i++) {
            designCopy[i] = new Integer[designs[i].length][];
            for (int j = 0; j < designs[i].length; j++) {
                designCopy[i][j] = new Integer[designs[i][j].length];
                System.arraycopy(designs[i][j], 0, designCopy[i][j], 0, designs[i][j].length);
            }
        }

        for (int i = 0; i < names.length; i++) {
            //AminoAcid thisAcid = AminoAcid.valueOf(names[i]);
            for (int j = 0; j < fragments[i].length; j++) {
                for (int k = fragments[i][j][0]; k <= fragments[i][j][1]; k++) {
                    int size = sequence.get(k).size();
                    for (int l = 0; l < designs[i][j].length; l++) {
                        while (sequence.get(k).size() < (l+1)*size) {
                            for (int m = 0; m < size; m++) {
                                sequence.get(k).add(sequence.get(k).get(m).clone());
                            }
                        }
                        int numLeft = 0;
                        for (int m = (l*size); m < (l+1)*size; m++) {
                            sequence.get(k);
                            numLeft = sequence.get(k).get(m).fill(codons[i], designs[i][j][l]);
                            sequence.get(k).get(m).deltas.put(codons[i].getAminoAcid(), designCopy[i][j][l]);
                        }
                        designs[i][j][l] = numLeft;
                    }
                }
            }
        }
        //TODO: make and check overlaps, fill in other codons
        return;
    }
}
