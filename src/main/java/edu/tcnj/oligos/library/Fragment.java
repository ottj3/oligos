package edu.tcnj.oligos.library;

import com.google.common.base.MoreObjects;
import com.google.common.collect.Lists;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Codon;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkState;

public class Fragment extends Sequence {

    private final int overlapLength;
    private final int oligoLength;
    private final Codon codon;
    private final Range range;
    private final Integer delta;
    private final Map<Integer, List<Oligo>> oligos;

    public Map<Integer, List<Oligo>> getOligos() {
        return oligos;
    }

    public Fragment(Sequence protein, Codon codon, Range range, Integer delta, Map<Integer,
                List<Oligo>> oligos, int oligoLength, int overlapLength) {
        super(protein);
        this.codon = codon;
        this.range = range;
        this.delta = delta;
        this.oligos = oligos;
        this.oligoLength = oligoLength;
        this.overlapLength = overlapLength;
    }

    @Override
    public Codon set(int index, Codon elem) {
        checkArgument(elem == codon);
        Codon temp = elem;
        AminoAcid acid = elem.getAminoAcid();
        for (int i = range.getStartPosition(); i < range.getEndPosition(); i++) {
            int offset = i - range.getStartPosition();
            int startOfPositionI = offset * (oligoLength - overlapLength);
            int positionInOligo = index - startOfPositionI;
            if (0 <= positionInOligo && positionInOligo < oligoLength) {
                for (Oligo oligo : oligos.get(i)) {
                    temp = oligo.set(positionInOligo, elem);
                    checkState(temp.getAminoAcid() == acid);
                }
            }
        }
        return temp;
    }

    void fill() {
        AminoAcid acid = codon.getAminoAcid();
        List<Integer> positionsOfInterest = Lists.newArrayList();
        for (int i = 0; i < sequence.size(); i++) {
            if (sequence.get(i).getAminoAcid() == acid) {
                positionsOfInterest.add(i);
            }
        }
        Collections.shuffle(positionsOfInterest);
        for (int i = 0; i < delta; i++) {
            set(i, codon);
        }
        Codon wildcard = acid.getWildcard();
        for (int i = delta; i < positionsOfInterest.size(); i++) {
            set(i, wildcard);
        }
    }

    static class Range {
        private int startPosition;
        private int endPosition;

        public Range(int startPosition, int endPosition) {
            this.startPosition = startPosition;
            this.endPosition = endPosition;
        }

        Sequence subSequence(Sequence sequence, int oligoLength, int overlapLength) {
            int start = startPosition * (oligoLength - overlapLength);
            int end = endPosition * (oligoLength - overlapLength) + oligoLength;
            return sequence.subList(start, end);
        }

        int getStartPosition() {
            return startPosition;
        }

        int getEndPosition() {
            return endPosition;
        }

        boolean contains(int pos) {
            return startPosition <= pos && pos <= endPosition;
        }

        @Override
        public String toString() {
            return MoreObjects.toStringHelper(this).add("start", startPosition).add("end", endPosition).toString();
        }

        @Override
        public boolean equals(Object other) {
            return (other instanceof Range
                    && ((Range) other).getStartPosition() == this.getStartPosition()
                    && ((Range) other).getEndPosition() == this.getEndPosition());
        }

        @Override
        public int hashCode() {
            return startPosition << 16 ^ endPosition;
        }
    }
}
