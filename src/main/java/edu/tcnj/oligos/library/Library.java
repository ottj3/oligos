package edu.tcnj.oligos.library;

import com.google.common.base.Strings;
import com.google.common.collect.EnumBiMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Codon;

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkNotNull;
import static com.google.common.base.Preconditions.checkState;

public class Library {

    private Protein protein;

    private int size; // number of positions
    private int oligoLength;
    private int overlapLength;

    private Map<Codon, Design> designs;
    private EnumBiMap<AminoAcid, Codon> codonsOfInterest = EnumBiMap.create(AminoAcid.class, Codon.class);
    private Map<Integer, List<Oligo>> oligos;

    private Library(Protein protein, int size, int oligoLength, int overlapLength, Map<Codon, Design> designs,
                    EnumBiMap<AminoAcid, Codon> codonsOfInterest) {
        this.protein = protein;
        this.size = size;
        this.oligoLength = oligoLength;
        this.overlapLength = overlapLength;
        this.designs = designs;
        this.codonsOfInterest = codonsOfInterest;
    }

    public Design getDesignForCodon(Codon codon) {
        return designs.get(codon);
    }

    public Design getDesignForAminoAcid(AminoAcid acid) {
        return designs.get(codonsOfInterest.get(acid));
    }

    public void createOligos() {
        if (oligos != null) return;
        oligos = Maps.newHashMap();
        for (int i = 0; i < size; i++) {
            createOligosForPosition(i);
        }
    }

    public FragmentIterator fragments() {
        return new FragmentIterator(designs, oligos, protein, oligoLength, overlapLength);
    }

    private void createOligosForPosition(int pos) {
        checkArgument(0 <= pos && pos < size);
        List<Map<Codon, Integer>> oligos = Lists.newArrayList();
        oligos.add(Maps.<Codon, Integer>newHashMap());
        for (Map.Entry<Codon, Design> entry : designs.entrySet()) {
            Design design = entry.getValue();
            List<Integer> deltasForPos = design.getDeltasForPosition(pos);
            if (deltasForPos != null) {
                int size = oligos.size();
                for (int i = 1; i < deltasForPos.size(); i++) {
                    for (Map<Codon, Integer> oligo : oligos) {
                        oligos.add(Maps.newHashMap(oligo));
                    }
                }
                for (int i = 0; i < deltasForPos.size(); i++) {
                    for (int j = i * size; j < (i + 1) * size; j++) {
                        Map<Codon, Integer> oligo = oligos.get(j);
                        oligo.put(entry.getKey(), deltasForPos.get(i));
                    }
                }
            }
        }
        List<Oligo> retOligos = Lists.newArrayList();
        for (Map<Codon, Integer> oligo : oligos) {
            int start = pos * (oligoLength - overlapLength);
            retOligos.add(new Oligo(protein.subList(start, start + oligoLength), oligo));
        }

        this.oligos.put(pos, retOligos);
    }

    private static class FragmentIterator implements Iterator<Fragment> {

        private final Iterator<Map.Entry<Codon, Design>> designs;
        private final Map<Integer, List<Oligo>> oligos;
        private final Protein protein;
        private final int oligoLength;
        private final int overlapLength;

        private Iterator<Map.Entry<Fragment.Range, List<Integer>>> ranges;
        private Iterator<Integer> deltas;

        private int delta;
        private Fragment.Range range;
        private Codon codon;

        FragmentIterator(Map<Codon, Design> designs, Map<Integer, List<Oligo>> oligos, Protein protein,
                         int oligoLength, int overlapLength) {
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
                    Map.Entry<Fragment.Range, List<Integer>> nextRange = ranges.next();
                    deltas = nextRange.getValue().iterator();

                    // set range and delta
                    range = nextRange.getKey();
                    delta = deltas.next();
                } else {
                    if (designs.hasNext()) {
                        // move designs forward, set new ranges and deltas
                        Map.Entry<Codon, Design> nextDesign = designs.next();
                        ranges = nextDesign.getValue().iterator();
                        Map.Entry<Fragment.Range, List<Integer>> nextRange = ranges.next();
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
                    filter(oligos, range, codon, delta), oligoLength, overlapLength);
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

        private static Map<Integer, List<Oligo>> filter(Map<Integer, List<Oligo>> oligoMap,
                                                        Fragment.Range range, Codon codon, int delta) {
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

    public static class Builder {
        private String proteinRNA = "";
        private int seqStart = -1;
        private int seqEnd;
        private int size = -1;
        private int oligoLength = -1;
        private int overlapSize = -1;
        private Map<Codon, Design> designs;
        private EnumBiMap<AminoAcid, Codon> codonsOfInterest;

        public Builder withSequenceLength(int start, int end) {
            checkArgument(start >= 0 && start < end,
                    "Invalid start and end positions for sequence: %s, %s", start, end);
            this.seqStart = start;
            this.seqEnd = end;
            return this;
        }
        public Builder withProteinFromRNA(String rna) {
            checkArgument(Strings.isNullOrEmpty(rna), "Empty input RNA.");
            checkArgument(rna.length() % 3 != 0, "Input RNA length not multiple of 3.");
            this.proteinRNA = rna;
            return this;
        }
        public Builder withPositions(int positions) {
            checkArgument(positions > 0, "Must have positive number of positions.");
            this.size = positions;
            return this;
        }
        public Builder withOligoSize(int length, int overlap) {
            checkArgument(length > 0 && length % 3 == 0, "Invalid oligo length %s", length);
            checkArgument(overlap > 0 && overlap < length && overlap % 3 == 0, "Invalid overlap size %s", overlap);
            this.oligoLength = length;
            this.overlapSize = overlap;
            return this;
        }
        public Builder withDesigns(Map<Codon, Design> designs) {
            checkNotNull(designs, "Can't have null designs.");
            this.designs = designs;
            return this;
        }
        public Builder withCodonsOfInterest(Map<AminoAcid, Codon> codons) {
            checkNotNull(codons, "Can't have null codons");
            this.codonsOfInterest = EnumBiMap.create(codons);
            return this;
        }
        public Library build() {
            checkState(!proteinRNA.isEmpty());
            checkState(designs != null);
            checkState(codonsOfInterest != null);
            checkState(oligoLength != -1);
            checkState(overlapSize != -1);
            Protein protein = (seqStart == -1)
                    ? new Protein(new Sequence(proteinRNA))
                    : new Protein(new Sequence(proteinRNA), seqStart, seqEnd);

            return new Library(protein, size, oligoLength, overlapSize, designs, codonsOfInterest);
        }
    }

}