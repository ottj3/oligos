package edu.tcnj.oligos.library;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import edu.tcnj.oligos.data.Base;
import edu.tcnj.oligos.data.Codon;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.List;

/**
 * A list of codons, used to store RNA sequences.
 * Also contains helper methods to determine if sequences match.
 */
public class Sequence extends AbstractList<Codon> {
    protected List<Codon> sequence;
    private BaseSequence bases;

    Sequence() {
    }

    //Operations on the underlying list
    @Override
    public int size() {
        return sequence.size();
    }

    @Override
    public Codon set(int index, Codon element) {
        int start = index * 3;
        List<Base> codonBase = element.toBases();
        for (int i = 0; i < 3; i++) {
            bases.set(i + start, codonBase.get(i));
        }
        return sequence.set(index, element);
    }

    @Override
    public void add(int index, Codon element) {
        int start = index * 3;
        List<Base> codonBase = element.toBases();
        for (int i = 0; i < 3; i++) {
            bases.add(start + i, codonBase.get(i));
        }
        sequence.add(index, element);
    }

    @Override
    public Codon remove(int index) {
        int start = index * 3;
        for (int i = start + 2; i >= start; i--) {
            bases.remove(i);
        }
        return sequence.remove(index);
    }

    @Override
    public Codon get(int i) {
        return sequence.get(i);
    }

    @Override
    public Sequence subList(int fromIndex, int toIndex) {
        return new Sequence(Lists.newArrayList(sequence.subList(fromIndex, toIndex)));
    }

    /**
     * Construct a sequence of codons from a string of RNA
     *
     * @param codons a string of bases/RNA
     */
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
        if (codons != null) asBases();
    }

    public void setSequence(List<Codon> seq) {
        this.sequence = seq;
        this.bases = null;
        if (seq != null) asBases();
    }

    public List<Codon> getSequence() {
        return sequence;
    }

    BaseSequence asBases() {
        if (this.bases == null) {
            List<Base> baseList = Lists.newArrayList();
            for (Codon codon : sequence) {
                baseList.addAll(codon.toBases());
            }
            this.bases = new BaseSequence(baseList);
        }
        return this.bases;
    }

    /**
     * Determine whether two sequences match
     *
     * @param seq1 the first sequence
     * @param seq2 the second sequence
     * @return boolean, true iff the sequences are identical
     */
    public static boolean regionsMatch(Sequence seq1, Sequence seq2) {
        return regionsMatch(seq1, 0, seq1.size(), seq2, 0, seq2.size());
    }

    /**
     * Determine whether two subsequences match
     *
     * @param seq1   the first sequence
     * @param start1 the first sequence's start position
     * @param end1   the first sequence's end position (exclusive)
     * @param seq2   the second sequence
     * @param start2 the second sequence's start position
     * @param end2   the second sequence's end position (exclusive)
     * @return boolean, true iff the sequences are identical
     */
    public static boolean regionsMatch(Sequence seq1, int start1, int end1, Sequence seq2, int start2, int end2) {
        // TODO potential non-1 number of differences needed
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
        return Joiner.on("").join(sequence);
    }
}
