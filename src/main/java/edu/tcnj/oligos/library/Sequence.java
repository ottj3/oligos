package edu.tcnj.oligos.library;

import com.google.common.base.MoreObjects;
import com.google.common.base.Preconditions;
import edu.tcnj.oligos.data.Codon;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.List;

import static com.google.common.base.Preconditions.checkElementIndex;

public class Sequence extends AbstractList<Codon> {
    protected List<Codon> sequence;

    Sequence() {}

    @Override
    public int size() {
        return sequence.size();
    }

    @Override
    public Codon set(int index, Codon element) {
        return sequence.set(index, element);
    }

    @Override
    public void add(int index, Codon element) {
        sequence.add(index, element);
    }

    @Override
    public Codon remove(int index) {
        return sequence.remove(index);
    }

    @Override
    public Codon get(int i) {
        return sequence.get(i);
    }

    @Override
    public Sequence subList(int fromIndex, int toIndex) {
        return new Sequence(sequence.subList(fromIndex, toIndex));
    }

    public Sequence(String codons) {
        if (codons.length() % 3 != 0) {
            throw new IllegalArgumentException("Sequence constructed with invalid number of bases.");
        }
        char[] codonArray = codons.toCharArray();
        List<Codon> tempSeq = new ArrayList<>(codonArray.length / 3);
        for (int i = 0; i < codonArray.length; i += 3) {
            Codon codon = Codon.valueOf(String.valueOf(codonArray, i, 3));
            tempSeq.add(codon);
        }
        this.sequence = tempSeq;
    }

    public Sequence(List<Codon> codons) {
        this.sequence = codons;
    }

    public void setSequence(List<Codon> seq) {
        this.sequence = seq;
    }

    public List<Codon> getSequence() {
        return sequence;
    }

    public static boolean regionsMatch(Sequence seq1, Sequence seq2) {
        return regionsMatch(seq1, 0, seq1.size(), seq2, 0, seq2.size());
    }

    // TODO potential non-1 number of differences needed
    public static boolean regionsMatch(Sequence seq1, int start1, int end1, Sequence seq2, int start2, int end2) {
        int length = end1 - start1;
        if (length != end2 - start2) {
            return false;
        }
        for (int i = 0; i < length; i++) {
            if (seq1.get(start1 + i) != seq2.get(start2 + i)) {
                return false;
            }
        }
        return true;
    }


    @Override
    public String toString() {
        return MoreObjects.toStringHelper(this)
                .add("seq", sequence.toArray())
                .toString();
    }
}
