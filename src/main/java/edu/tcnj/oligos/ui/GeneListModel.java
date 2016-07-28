package edu.tcnj.oligos.ui;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Codon;
import edu.tcnj.oligos.library.*;

import java.util.List;
import java.util.Map;

class GeneListModel extends SequenceListModel {
    private List<Gene> genes;

    void addGenes(Library lib) {
        super.removeAllElements();

        List<Sequence> geneSeqs = LibraryUtils.buildPermutations(new Fragment.Range(0, lib.getSize() - 1),
                lib.getOligos(), lib.getOligoLength(), lib.getOverlapLength());
        this.genes = Lists.newArrayListWithCapacity(geneSeqs.size());

        Map<AminoAcid, Codon> coi = lib.getCodonsOfInterest();
        for (Sequence gene : geneSeqs) {
            this.genes.add(Gene.fromSequence(gene, coi));
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