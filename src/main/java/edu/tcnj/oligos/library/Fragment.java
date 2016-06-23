package edu.tcnj.oligos.library;

import com.google.common.base.MoreObjects;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Codon;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkState;

public class Fragment extends Sequence {

    private final int overlapLength;
    private final int oligoLength;
    private final int size;
    private final Codon codon;
    private final Range range;
    private final Integer delta;
    private final Map<Integer, List<Oligo>> oligos;

    public Map<Integer, List<Oligo>> getOligos() {
        return oligos;
    }

    private Fragment(Sequence protein, Codon codon, Range range, Integer delta, Map<Integer,
            List<Oligo>> oligos, int oligoLength, int overlapLength, int size) {
        super(protein);
        this.codon = codon;
        this.range = range;
        this.delta = delta;
        this.oligos = oligos;
        this.oligoLength = oligoLength;
        this.overlapLength = overlapLength;
        this.size = size;
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
        int start = (range.getStartPosition() == 0 ? 0 : overlapLength);
        int end = (range.getEndPosition() == size - 1 ? sequence.size() : sequence.size() - overlapLength);
        for (int i = start; i < end; i++) {
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

        public boolean contains(Range other) {
            return (this.startPosition <= other.getStartPosition() && other.getEndPosition() <= this.endPosition);
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

    static class FragmentIterator implements Iterator<Fragment> {

        private final Iterator<Map.Entry<Codon, Design>> designs;
        private final Map<Integer, List<Oligo>> oligos;
        private final Protein protein;
        private final int oligoLength;
        private final int overlapLength;
        private final int size;

        private Iterator<Map.Entry<Range, List<Integer>>> ranges;
        private Iterator<Integer> deltas;

        private int delta;
        private Range range;
        private Codon codon;

        FragmentIterator(Map<Codon, Design> designs, Map<Integer, List<Oligo>> oligos, Protein protein,
                         int oligoLength, int overlapLength, int size) {
            this.designs = designs.entrySet().iterator();
            checkState(this.designs.hasNext());
            this.ranges = this.designs.next().getValue().iterator();
            checkState(this.ranges.hasNext());
            this.deltas = this.ranges.next().getValue().iterator();
            checkState(this.deltas.hasNext());

            this.oligos = oligos;
            this.protein = protein;
            this.oligoLength = oligoLength;
            this.overlapLength = overlapLength;
            this.size = size;
        }

        @Override
        public boolean hasNext() {
            return deltas.hasNext() || ranges.hasNext() || designs.hasNext();
        }

        @Override
        public Fragment next() {
            if (deltas.hasNext()) {
                // next delta!
                delta = deltas.next();
            } else {
                if (ranges.hasNext()) {
                    // move ranges forward, set new deltas
                    Map.Entry<Range, List<Integer>> nextRange = ranges.next();
                    deltas = nextRange.getValue().iterator();

                    // set range and delta
                    range = nextRange.getKey();
                    delta = deltas.next();
                } else {
                    if (designs.hasNext()) {
                        // move designs forward, set new ranges and deltas
                        Map.Entry<Codon, Design> nextDesign = designs.next();
                        ranges = nextDesign.getValue().iterator();
                        Map.Entry<Range, List<Integer>> nextRange = ranges.next();
                        deltas = nextRange.getValue().iterator();

                        // set codon, range, and delta
                        codon = nextDesign.getKey();
                        range = nextRange.getKey();
                        delta = deltas.next();
                    } else {
                        throw new NoSuchElementException();
                    }
                }
            }

            return new Fragment(range.subSequence(protein, oligoLength, overlapLength), codon, range, delta,
                    filter(oligos, range, codon, delta), oligoLength, overlapLength, size);
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

        private static Map<Integer, List<Oligo>> filter(Map<Integer, List<Oligo>> oligoMap,
                                                        Range range, Codon codon, int delta) {
            Map<Integer, List<Oligo>> map = Maps.newHashMap();
            for (int i = range.getStartPosition(); i < range.getEndPosition(); i++) {
                List<Oligo> oligoList = oligoMap.get(i);
                for (Oligo oligo : oligoList) {
                    Map<Codon, Integer> deltas = oligo.getDeltas();
                    if (deltas.containsKey(codon)) {
                        if (deltas.get(codon) == delta) {
                            if (map.get(i) == null) {
                                List<Oligo> oligos = Lists.newArrayList();
                                oligos.add(oligo);
                                map.put(i, oligos);
                            } else {
                                map.get(i).add(oligo);
                            }
                        }
                    }
                }
            }
            return map;
        }
    }
}
