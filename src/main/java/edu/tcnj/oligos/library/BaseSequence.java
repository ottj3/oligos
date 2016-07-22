package edu.tcnj.oligos.library;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import edu.tcnj.oligos.data.Base;

import java.util.AbstractList;
import java.util.List;

public class BaseSequence extends AbstractList<Base> {

    private List<Base> sequence;

    public BaseSequence(List<Base> sequence) {
        this.sequence = sequence;
    }

    public int posOfMatch(BaseSequence other) {
        for (int i = 0; i < this.sequence.size(); i++) {
            for (int j = 0; j < other.sequence.size(); j++) {
                if (i + j >= this.sequence.size()) return -1;
                if (this.get(i + j).matches(other.get(j))) {
                    if (j == other.sequence.size() - 1) {
                        return i;
                    }
                } else {
                    break;
                }
            }
        }
        return -1;
    }

    //Operations on the underlying list
    @Override
    public int size() {
        return sequence.size();
    }

    @Override
    public Base set(int index, Base element) {
        return sequence.set(index, element);
    }

    @Override
    public void add(int index, Base element) {
        sequence.add(index, element);
    }

    @Override
    public Base remove(int index) {
        return sequence.remove(index);
    }

    @Override
    public Base get(int i) {
        return sequence.get(i);
    }

    @Override
    public BaseSequence subList(int fromIndex, int toIndex) {
        return new BaseSequence(Lists.newArrayList(sequence.subList(fromIndex, toIndex)));
    }

    @Override
    public String toString() {
        return Joiner.on("").join(sequence);
    }

    public static int numDifferences(BaseSequence seq1, BaseSequence seq2) {
        return numDifferences(seq1, 0, seq1.size(), seq2, 0, seq2.size());
    }
    public static int numDifferences(BaseSequence seq1, int start1, int end1, BaseSequence seq2, int start2, int end2) {
        int numDifferences;
        int length;
        if (end1 - start1 > end2 - start2) {
            length = end2 - start2;
            numDifferences = (end1 - start1) - length;
        } else {
            length = end1 - start1;
            numDifferences = (end2 - start2) - length;
        }
        for (int i = 0; i < length; i++) {
            if (seq1.get(start1 + i) != seq2.get(start2 + i)) {
                numDifferences++;
            }
        }
        return numDifferences;
    }
}
