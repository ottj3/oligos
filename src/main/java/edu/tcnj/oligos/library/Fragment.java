package edu.tcnj.oligos.library;

import com.google.common.base.MoreObjects;
import com.google.common.collect.Collections2;
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

/**
 * A Fragment corresponds to a given range, and a specific
 * delta for a codons of interest varied in that range.
 * Also tracks and works on all oligos that have that
 * specific delta for that codon. Used to fill these
 * oligos so that they contribute towards the desired
 * frequencies for the codons of interest
 */
public class Fragment extends Sequence {

    private final int overlapLength;
    private final int oligoLength;
    private final int size;
    private final Codon codon;
    private final Range range;
    private final Integer delta;
    private final Map<Integer, List<Oligo>> oligos;
    private final Map<Codon, Double> codonFreqs;
    private final List<BaseSequence> restrictions;

    public Map<Integer, List<Oligo>> getOligos() {
        return oligos;
    }

    private Fragment(Sequence protein, Codon codon, Range range, Integer delta,
                     Map<Integer, List<Oligo>> oligos, int oligoLength, int overlapLength,
                     int size, Map<Codon, Double> codonFreqs, List<BaseSequence> restrictions) {
        super(protein);
        this.codon = codon;
        this.range = range;
        this.delta = delta;
        this.oligos = oligos;
        this.oligoLength = oligoLength;
        this.overlapLength = overlapLength;
        this.size = size;
        this.codonFreqs = codonFreqs;
        this.restrictions = restrictions;
    }

    @Override
    public Codon set(int index, Codon elem) {
        checkArgument(elem.getAminoAcid() == codon.getAminoAcid());
        Codon temp = elem;
        AminoAcid acid = elem.getAminoAcid();
        for (int i = range.getStartPosition(); i <= range.getEndPosition(); i++) {
            int offset = i - range.getStartPosition();
            int smalligo = oligoLength - overlapLength;
            int startOfPositionI = offset * smalligo;
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

    /**
     * Fills a fragment with delta of the codon of interest,
     * and proportional numbers of all other codons for that
     * acid (based on the frequency table built by the Library)
     * arranged in random order
     */
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
        Map<Codon, Integer> counts = Library.findCodonCounts(codonFreqs, positionsOfInterest.size() - delta);
        Collections.shuffle(positionsOfInterest);
        Iterator<List<Integer>> perm = Collections2.permutations(positionsOfInterest).iterator();
        do {
            //Permute the positions and fill them; if this happened to
            //make any restriction enzyme sites try a different permutation
            if (!perm.hasNext()) {
                throw new RuntimeException(new OutOfSwapsException("Ran out of permutations when filling fragment: "
                        + this.toString()));
            }
            positionsOfInterest = perm.next();
            for (int i = 0; i < delta; i++) {
                set(positionsOfInterest.get(i), codon);
            }
            int index = delta;
            for (Map.Entry<Codon, Integer> entry : counts.entrySet()) {
                for (int j = 0; j < entry.getValue(); j++) {
                    set(positionsOfInterest.get(index), entry.getKey());
                    index++;
                }
            }
        }
        while (LibraryUtils.containsRestrictionEnzyme(
                LibraryUtils.buildPermutations(range, oligos, oligoLength, overlapLength), restrictions));
    }

    @Override
    public String toString() {
        return codon.toString() + " (" + range.getStartPosition() + ","
                + range.getEndPosition() + ") -> " + delta;// + "\n" + super.toString();
    }

    /**
     * Stores a starting and ending position (both inclusive),
     * each position corresponds to an oligo position in the overall protein.
     */
    public static class Range {
        private int startPosition;
        private int endPosition;

        public Range(int startPosition, int endPosition) {
            this.startPosition = startPosition;
            this.endPosition = endPosition;
        }

        //Get the subsequence of the original RNA sequence that corresponds to this range
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

    /**
     * Helps to iterate through every fragment for filling them.
     * From designs, it iterates through every codon. For every codon,
     * it iterates through every range. For every range, it iterates
     * through every delta level.
     */
    static class FragmentIterator implements Iterator<Fragment> {

        // internal iterators
        private final Iterator<Map.Entry<Codon, Design>> designs;
        private Iterator<Map.Entry<Range, List<Integer>>> ranges;
        private Iterator<Integer> deltas;
        // iterated
        private int delta;
        private Codon codon;
        private Range range;


        // pass through to fragment constructor
        private final Protein protein;
        private final Map<Integer, List<Oligo>> oligos;
        private final int oligoLength;
        private final int overlapLength;
        private final int size;
        private final Map<AminoAcid, Map<Codon, Double>> codonFrequencies;
        private final List<BaseSequence> restrictions;

        //Constructs the iterator and sets up the first fragment's info
        FragmentIterator(Map<Codon, Design> designs, Map<Integer, List<Oligo>> oligos, Protein protein,
                         int oligoLength, int overlapLength, int size,
                         Map<AminoAcid, Map<Codon, Double>> codonFrequencies, List<BaseSequence> restrictions) {
            this.designs = designs.entrySet().iterator();
            checkState(this.designs.hasNext());
            Map.Entry<Codon, Design> firstDesign = this.designs.next();
            codon = firstDesign.getKey();

            this.ranges = firstDesign.getValue().iterator();
            checkState(this.ranges.hasNext());
            Map.Entry<Range, List<Integer>> firstRange = this.ranges.next();
            range = firstRange.getKey();

            this.deltas = firstRange.getValue().iterator();
            checkState(this.deltas.hasNext());

            this.oligos = oligos;
            this.protein = protein;
            this.oligoLength = oligoLength;
            this.overlapLength = overlapLength;
            this.size = size;

            this.codonFrequencies = codonFrequencies;
            this.restrictions = restrictions;
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
                    filter(oligos, range, codon, delta), oligoLength, overlapLength, size,
                    codonFrequencies.get(codon.getAminoAcid()), restrictions);
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

        //Filter out the overall position->list<oligo> map to only oligos from
        //positions in the given range with the given codon->delta pairing
        private static Map<Integer, List<Oligo>> filter(Map<Integer, List<Oligo>> oligoMap,
                                                        Range range, Codon codon, int delta) {
            Map<Integer, List<Oligo>> map = Maps.newHashMap();
            for (int i = range.getStartPosition(); i <= range.getEndPosition(); i++) {
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
