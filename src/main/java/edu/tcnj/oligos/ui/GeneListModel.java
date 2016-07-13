package edu.tcnj.oligos.ui;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Codon;
import edu.tcnj.oligos.library.Fragment;
import edu.tcnj.oligos.library.Gene;
import edu.tcnj.oligos.library.Library;
import edu.tcnj.oligos.library.LibraryUtils;
import edu.tcnj.oligos.library.Sequence;

import java.util.List;
import java.util.Map;

class GeneListModel extends SequenceListModel {
    private List<Gene> genes;

    GeneListModel() {
    }

    void addGenes(Library lib) {
        super.removeAllElements();

        List<Sequence> geneSeqs = LibraryUtils.buildPermutations(new Fragment.Range(0, lib.getSize() - 1),
                lib.getOligos(), lib.getOligoLength(), lib.getOverlapLength());
        this.genes = Lists.newArrayListWithCapacity(geneSeqs.size());

        for (Sequence gene : geneSeqs) {
            Map<AminoAcid, Map<Codon, Integer>> counts = Maps.newHashMap();
            for (AminoAcid acid : lib.getCodonsOfInterest().keySet()) {
                counts.put(acid, Maps.<Codon, Integer>newHashMap());
            }
            for (Codon codon : gene) {
                AminoAcid acid = codon.getAminoAcid();
                if (lib.getCodonsOfInterest().containsKey(acid)) {
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
            this.genes.add(new Gene(gene, freqs));
        }

        for (int i = 0; i < genes.size(); i++) {
            Gene gene = genes.get(i);
            String name = (i + 1) + " [";
            List<String> codonFreqs = Lists.newArrayList();
            for (Map.Entry<AminoAcid, Map<Codon, Double>> acidEntry : gene.getFreqs().entrySet()) {
                for (Map.Entry<Codon, Double> codonEntry : acidEntry.getValue().entrySet()) {
                    if (lib.getCodonsOfInterest().get(acidEntry.getKey()) != codonEntry.getKey()) {
                        continue;
                    }
                    codonFreqs.add(codonEntry.getKey() + " (" + acidEntry.getKey() + "): "
                            + ((int) (codonEntry.getValue() * 100 + 0.5) + "%"));
                }
            }
            name += Joiner.on(", ").join(codonFreqs) + "]";
            super.addElement(name);
        }
    }

    @Override
    Gene getActualAt(int i) {
        return genes.get(i);
    }
}