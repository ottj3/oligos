package edu.tcnj.oligos.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class RNASequence {
    public List<Codon> data = new ArrayList<>();
    public Map<AminoAcid, Integer> deltas = new HashMap<>();

    String label;

    public List<RNASequence> forwardEdges = new ArrayList<>();
    public List<RNASequence> backEdges = new ArrayList<>();

    public RNASequence(String label) {
        this.label = label;
    }

    public RNASequence clone() {
        RNASequence newRNASequence = new RNASequence(this.label);
        newRNASequence.data = new ArrayList<>(this.data);
        newRNASequence.deltas = new HashMap<>(this.deltas);
        return newRNASequence;
    }

    public int fill(Codon newCodon, int numNeeded) {
        int numOccurrences = 0;
        for (Codon codon : data) {
            if (codon == newCodon) {
                numOccurrences++;
            }
        }
        for (int i = 0; i < data.size(); i++) {
            if (data.get(i).getAminoAcid() == newCodon.getAminoAcid()) {
                if (numOccurrences < numNeeded && data.get(i) != newCodon) {
                    data.set(i, newCodon);
                    numOccurrences++;
                } else if (numOccurrences > numNeeded) {
                    data.set(i, data.get(i).getAminoAcid().getWildcard());
                }
            }
        }
        return numNeeded - numOccurrences;
    }
}
