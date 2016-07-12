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

    public int contains(BaseSequence other) {
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
}
