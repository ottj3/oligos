package edu.tcnj.oligos.library;

import com.google.common.collect.Maps;
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

    public static Gene fromSequence(Sequence gene, Map<AminoAcid, Codon> coi) {
        Map<AminoAcid, Map<Codon, Integer>> counts = Maps.newHashMap();
        for (AminoAcid acid : coi.keySet()) {
            counts.put(acid, Maps.<Codon, Integer>newHashMap());
        }
        for (Codon codon : gene) {
            AminoAcid acid = codon.getAminoAcid();
            if (coi.containsKey(acid)) {
                Map<Codon, Integer> codonMap = counts.get(acid);
                if (codonMap.containsKey(codon)) {
                    codonMap.put(codon, codonMap.get(codon) + 1);
                } else {
                    codonMap.put(codon, 1);
                }
            }
        }
        Map<AminoAcid, Map<Codon, Double>> freqs = Maps.newHashMap();
        for (Map.Entry<AminoAcid, Map<Codon, Integer>> acidEntry : counts.entrySet()) {
            freqs.put(acidEntry.getKey(), Maps.<Codon, Double>newHashMap());
            double totalForAcid = 0;
            for (Map.Entry<Codon, Integer> codonEntry : acidEntry.getValue().entrySet()) {
                totalForAcid += codonEntry.getValue();
            }
            for (Map.Entry<Codon, Integer> codonEntry : acidEntry.getValue().entrySet()) {
                freqs.get(acidEntry.getKey()).put(codonEntry.getKey(), codonEntry.getValue() / totalForAcid);
            }
        }
        for (Map.Entry<AminoAcid, Codon> entry : coi.entrySet()) {
            if (!freqs.get(entry.getKey()).containsKey(entry.getValue())) {
                freqs.get(entry.getKey()).put(entry.getValue(), 0.0D);
            }
        }
        return new Gene(gene, freqs);
    }
}
