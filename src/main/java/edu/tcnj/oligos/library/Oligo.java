package edu.tcnj.oligos.library;

import edu.tcnj.oligos.data.Codon;

import java.util.Collections;
import java.util.List;
import java.util.Map;

public class Oligo extends Sequence {

    private Map<Codon, Integer> deltas;

    public Oligo(List<Codon> codons, Map<Codon, Integer> deltas) {
        super(codons);
        this.deltas = Collections.unmodifiableMap(deltas);
    }

    public Map<Codon, Integer> getDeltas() {
        return deltas;
    }

}
