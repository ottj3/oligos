package edu.tcnj.oligos.library;

import com.google.common.collect.Iterables;
import edu.tcnj.oligos.data.Codon;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;

import static com.google.common.base.Preconditions.checkState;

public class Overlap extends Oligo {
    private List<Oligo> preAttachments;
    private List<Oligo> postAttachments;
    private int overlapLength = -1;

    public Overlap(Sequence codons, Map<Codon, Integer> deltas) {
        super(codons, deltas);
        this.preAttachments = new ArrayList<>();
        this.postAttachments = new ArrayList<>();
    }

    public void linkPreAttachment(Oligo preattach) {
        preAttachments.add(preattach);
    }

    public void linkPostAttachment(Oligo postattach) {
        postAttachments.add(postattach);
    }

    public void finalizeLink(int overlapLength) {
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

    @Override
    public Codon set(int index, Codon codon) {
        Codon prev = sequence.set(index, codon);
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

    public void swapWithPreAttachments(int overlapIndex, int preAttachIndex) {
        Codon current = get(overlapIndex);
        Codon target = preAttachments.get(0).get(preAttachIndex);
        checkState(current.getAminoAcid() == target.getAminoAcid());
        for (Oligo preAttachment : preAttachments) {
            preAttachment.set(preAttachIndex, current);
        }
        set(overlapIndex, target);
    }

    public List<Oligo> getPreAttachments() {
        return preAttachments;
    }

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
