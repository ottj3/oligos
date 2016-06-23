package edu.tcnj.oligos.library;

import com.google.common.base.Strings;
import com.google.common.collect.EnumBiMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Codon;
import edu.tcnj.oligos.library.Fragment.FragmentIterator;

import java.util.List;
import java.util.Map;
import java.util.Objects;

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
    private Map<Integer, List<Overlap>> overlaps;

    private Library(Protein protein, int size, int oligoLength, int overlapLength, Map<Codon, Design> designs,
                    EnumBiMap<AminoAcid, Codon> codonsOfInterest) {
        this.protein = protein;
        this.size = size;
        this.oligoLength = oligoLength;
        this.overlapLength = overlapLength;
        this.designs = designs;
        this.codonsOfInterest = codonsOfInterest;
    }

    public void createOligos() {
        if (oligos != null) return;
        oligos = Maps.newHashMap();
        for (int i = 0; i < size; i++) {
            createOligosForPosition(i);
        }
    }

    public FragmentIterator fragmentIterator() {
        return new FragmentIterator(designs, oligos, protein, oligoLength, overlapLength, size);
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

    private void createOverlapsForPosition(int pos) {
        checkArgument(0 <= pos && pos < size - 1);
        List<Map<Codon, Integer>> overlaps = Lists.newArrayList();
        overlaps.add(Maps.<Codon, Integer>newHashMap());
        for (Map.Entry<Codon, Design> entry : designs.entrySet()) {
            Design design = entry.getValue();
            List<Integer> deltasForPos = design.getDeltasForRange(new Fragment.Range(pos, pos + 1));
            if (deltasForPos != null) {
                int size = overlaps.size();
                for (int i = 1; i < deltasForPos.size(); i++) {
                    for (Map<Codon, Integer> overlap : overlaps) {
                        overlaps.add(Maps.newHashMap(overlap));
                    }
                }
                for (int i = 0; i < deltasForPos.size(); i++) {
                    for (int j = i * size; j < (i + 1) * size; j++) {
                        Map<Codon, Integer> overlap = overlaps.get(j);
                        overlap.put(entry.getKey(), deltasForPos.get(i));
                    }
                }
            }
        }
        List<Overlap> retOverlaps = Lists.newArrayList();
        for (Map<Codon, Integer> deltas : overlaps) {
            Overlap overlap = new Overlap(null, deltas);
            for (Oligo oligo : oligos.get(pos)) {
                for (Map.Entry<Codon, Integer> entry : deltas.entrySet()) {
                    if (Objects.equals(deltas.get(entry.getKey()), oligo.getDeltas().get(entry.getKey()))) {
                        overlap.linkPreAttachment(oligo);
                    }
                }
            }
            for (Oligo oligo : oligos.get(pos + 1)) {
                for (Map.Entry<Codon, Integer> entry : deltas.entrySet()) {
                    if (Objects.equals(deltas.get(entry.getKey()), oligo.getDeltas().get(entry.getKey()))) {
                        overlap.linkPostAttachment(oligo);
                    }
                }
            }
            overlap.finalizeLink(overlapLength);
        }

        this.overlaps.put(pos, retOverlaps);
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