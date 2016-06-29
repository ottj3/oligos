package edu.tcnj.oligos.library;

import edu.tcnj.oligos.data.Codon;

import java.util.Collections;
import java.util.Map;

/**
 * A sequence set up for specific codon->delta mappings
 */
public class Oligo extends Sequence {

    private Map<Codon, Integer> deltas;

    public Oligo(Sequence codons, Map<Codon, Integer> deltas) {
        super(codons);
        this.deltas = Collections.unmodifiableMap(deltas);
    }

    public Map<Codon, Integer> getDeltas() {
        return deltas;
    }

}
