package edu.tcnj.oligos.library;

import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Codon;

import java.util.Collections;
import java.util.Map;

public class Gene extends Sequence {
    private Map<AminoAcid, Map<Codon, Double>> freqs;

    public Gene(Sequence codons, Map<AminoAcid, Map<Codon, Double>> freqs) {
        super(codons);
        this.freqs = Collections.unmodifiableMap(freqs);
    }

    public Map<AminoAcid, Map<Codon, Double>> getFreqs() {
        return freqs;
    }

}
