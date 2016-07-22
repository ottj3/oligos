package edu.tcnj.oligos.library;

import edu.tcnj.oligos.data.Codon;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;

import static com.google.common.base.Preconditions.checkState;

/**
 * A sequence that is part of neighboring oligos. Stores the oligos that should connect to it on both sides.
 */
public class Overlap extends Oligo {
    private List<Oligo> preAttachments;
    private List<Oligo> postAttachments;
    private int overlapLength = -1;

    public Overlap(Sequence codons, Map<Codon, Integer> deltas) {
        super(codons, deltas);
        this.preAttachments = new ArrayList<>();
        this.postAttachments = new ArrayList<>();
    }

    void linkPreAttachment(Oligo preAttach) {
        preAttachments.add(preAttach);
    }

    void linkPostAttachment(Oligo postAttach) {
        postAttachments.add(postAttach);
    }

    //Set the contents of the overlap, and make sure that it matches all of the oligos linked to it
    void finalizeLink(int overlapLength) {
        checkState(this.overlapLength == -1);
        checkState(!this.preAttachments.isEmpty());
        checkState(!this.postAttachments.isEmpty());

        this.overlapLength = overlapLength;
        Oligo oligo = preAttachments.get(0);
        this.setSequence(oligo.subList(oligo.size() - overlapLength, oligo.size()));
        for (Oligo attachment : preAttachments) {
            checkState(Sequence.regionsMatch(this, 0, this.size(),
                    attachment, attachment.size() - overlapLength, attachment.size()),
                    "Linked to non-matching pre oligo: %s, %s", this, attachment);
        }
        for (Oligo attachment : postAttachments) {
            checkState(Sequence.regionsMatch(this, 0, this.size(),
                    attachment, 0, this.size()),
                    "Linked to non-matching post oligo: %s, %s", this, attachment);
        }
    }

    //Set the codon in all oligos that connect to this overlap (to keep the overlap regions matching)
    @Override
    public Codon set(int index, Codon codon) {
        Codon prev = super.set(index, codon);
        for (Oligo oligo : preAttachments) {
            Codon prevO = oligo.set(oligo.size() - this.size() + index, codon);
            checkState(prevO == prev);
        }
        for (Oligo oligo : postAttachments) {
            Codon prevO = oligo.set(index, codon);
            checkState(prevO == prev);
        }
        return prev;
    }

    /**
     * Swap a codon from an overlap into either its pre or post attachments.
     */
    public void swapWithAttachments(int overlapIndex, int attachIndex, List<Oligo> attachments) {
        Codon current = get(overlapIndex);
        Codon target = attachments.get(0).get(attachIndex);
        checkState(current.getAminoAcid() == target.getAminoAcid());
        for (Oligo attachment : attachments) {
            if (attachment.get(attachIndex) != target) {
                System.out.println("Codon at pre attach index does not match in swap. " +
                        "Relative frequencies may be skewed.");
            }
            attachment.set(attachIndex, current);
        }
        set(overlapIndex, target);
    }

    public List<Oligo> getPreAttachments() {
        return preAttachments;
    }

    public List<Oligo> getPostAttachments() {
        return postAttachments;
    }

    @Override
    public String toString() {
        String ret = "";
        for (Map.Entry<Codon, Integer> entry : this.getDeltas().entrySet()) {
            ret += entry.getKey() + " (" + entry.getKey().getAminoAcid() + ") -> " + entry.getValue() + "\n";
        }
        return ret;
    }

    /**
     * Helps to iterate through all overlaps for making sure they are unique.
     * Iterates through every overlap position,
     * and through every overlap in that position.
     */
    static class OverlapIterator implements Iterator<Overlap> {
        private Iterator<Map.Entry<Integer, List<Overlap>>> positionIterator;
        private int currentPosition;
        private Iterator<Overlap> overlapIterator;

        public OverlapIterator(Map<Integer, List<Overlap>> overlaps) {
            this.positionIterator = overlaps.entrySet().iterator();
            checkState(this.positionIterator.hasNext());
            Map.Entry<Integer, List<Overlap>> nextPosition = this.positionIterator.next();
            this.overlapIterator = nextPosition.getValue().iterator();
            checkState(this.overlapIterator.hasNext());
        }

        @Override
        public boolean hasNext() {
            return overlapIterator.hasNext() || positionIterator.hasNext();
        }

        @Override
        public Overlap next() {
            if (overlapIterator.hasNext()) {
                return overlapIterator.next();
            } else {
                if (positionIterator.hasNext()) {
                    //go to the next position, set up the overlap iterator for this position
                    Map.Entry<Integer, List<Overlap>> nextPos = positionIterator.next();
                    this.overlapIterator = nextPos.getValue().iterator();

                    this.currentPosition = nextPos.getKey();
                    if (overlapIterator.hasNext()) {
                        return overlapIterator.next();
                    } else {
                        throw new NoSuchElementException();
                    }
                } else {
                    throw new NoSuchElementException();
                }
            }
        }

        public int getCurrentPosition() {
            return currentPosition;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }
}
