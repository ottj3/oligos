package edu.tcnj.oligos.library;

import edu.tcnj.oligos.data.Codon;

import java.util.ArrayList;
import java.util.List;

import static com.google.common.base.Preconditions.checkState;

public class Overlap extends Sequence {
    private List<Oligo> preAttachements;
    private List<Oligo> postAttachments;

    public Overlap(List<Codon> codons) {
        super(codons);
        this.preAttachements = new ArrayList<>();
        this.postAttachments = new ArrayList<>();
    }

    public void linkPreAttachment(Oligo preattach) {
        preAttachements.add(preattach);
    }

    public void linkPostAttachment(Oligo postattach) {
        postAttachments.add(postattach);
    }

    public List<Oligo> getPreAttachements() {
        return preAttachements;
    }

    public List<Oligo> getPostAttachments() {
        return postAttachments;
    }

    @Override
    public Codon set(int index, Codon codon) {
        Codon prev = sequence.set(index, codon);
        for (Oligo oligo : preAttachements) {
            Codon prevO = oligo.set(oligo.size() - this.size() + index, codon);
            checkState(prevO == prev, "Sequence was changed directly after adding to Overlap region.");
        }
        for (Oligo oligo : postAttachments) {
            Codon prevO = oligo.set(index, codon);
            checkState(prevO == prev, "Sequence was changed directly after adding to Overlap region.");
        }
        return prev;
    }
}
