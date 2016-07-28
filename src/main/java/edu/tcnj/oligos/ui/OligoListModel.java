package edu.tcnj.oligos.ui;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import edu.tcnj.oligos.data.Codon;
import edu.tcnj.oligos.library.Oligo;

import java.util.List;
import java.util.Map;

class OligoListModel extends SequenceListModel {
    private Map<Integer, List<Oligo>> oligos;

    void setOligos(Map<Integer, List<Oligo>> oligos) {
        super.removeAllElements();
        this.oligos = oligos;
        for (Map.Entry<Integer, List<Oligo>> entry : oligos.entrySet()) {
            char var = 'a';
            for (Oligo oligo : entry.getValue()) {
                List<String> deltaList = Lists.newArrayList();
                for (Map.Entry<Codon, Integer> deltas : oligo.getDeltas().entrySet()) {
                    deltaList.add(deltas.getKey() + " (" + deltas.getKey().getAminoAcid() + "): " + deltas.getValue());
                }
                super.addElement(entry.getKey() + " (" + var++ + "): [" + Joiner.on(", ").join(deltaList) + "]");
            }
        }
    }

    @Override
    Oligo getActualAt(int i) {
        int j = 0;
        for (Map.Entry<Integer, List<Oligo>> entry : oligos.entrySet()) {
            for (Oligo oligo : entry.getValue()) {
                if (j == i) return oligo;
                j++;
            }
        }
        return null;
    }
}
