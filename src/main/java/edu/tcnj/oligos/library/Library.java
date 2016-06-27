package edu.tcnj.oligos.library;

import com.google.common.base.Strings;
import com.google.common.collect.EnumBiMap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multiset;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Codon;
import edu.tcnj.oligos.library.Fragment.FragmentIterator;
import edu.tcnj.oligos.library.Overlap.OverlapIterator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkNotNull;
import static com.google.common.base.Preconditions.checkState;

public class Library {

    private Protein protein;

    private final int size; // number of positions
    private final int oligoLength;
    private final int overlapLength;
    private final int smalligo;

    private Map<Codon, Design> designs;
    private EnumBiMap<AminoAcid, Codon> codonsOfInterest = EnumBiMap.create(AminoAcid.class, Codon.class);
    private final Map<AminoAcid, Map<Codon, Double>> codonFrequencies;
    private Map<Integer, List<Oligo>> oligos;
    private Map<Integer, List<Overlap>> overlaps;

    private Library(Protein protein, int size, int oligoLength, int overlapLength, Map<Codon, Design> designs,
                    EnumBiMap<AminoAcid, Codon> codonsOfInterest) {
        this.protein = protein;
        this.size = size;
        this.oligoLength = oligoLength;
        this.overlapLength = overlapLength;
        this.smalligo = oligoLength - overlapLength;
        this.designs = designs;
        this.codonsOfInterest = codonsOfInterest;

        this.codonFrequencies = calcFrequencies();
    }

    private Map<AminoAcid, Map<Codon, Double>> calcFrequencies() {
        Map<AminoAcid, Multiset<Codon>> counts = Maps.newHashMap();
        for (Codon codon : protein) {
            AminoAcid acid = codon.getAminoAcid();
            if (codonsOfInterest.containsKey(acid)) {
                Multiset<Codon> codons = counts.get(acid);
                if (codons == null) {
                    codons = HashMultiset.create();
                    counts.put(acid, codons);
                    codons.add(codon);
                } else {
                    codons.add(codon);
                }
            }
        }

        Map<AminoAcid, Map<Codon, Double>> frequencies = Maps.newHashMap();
        for (Map.Entry<AminoAcid, Multiset<Codon>> entry : counts.entrySet()) {
            AminoAcid acid = entry.getKey();
            frequencies.put(acid, Maps.<Codon, Double>newHashMap());
            Multiset<Codon> codons = entry.getValue();
            int total = codons.size();
            for (Codon possible : Codon.getCodonsForAcid(acid)) {
                int number = codons.count(possible);
                double freq = number / total;
                frequencies.get(acid).put(possible, freq);
            }
        }
        return frequencies;
    }

    public void createOligos() {
        if (oligos != null) return;
        oligos = Maps.newHashMap();
        for (int i = 0; i < size; i++) {
            createOligosForPosition(i);
        }
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
                    for (int j = 0; j < size; j++) {
                        oligos.add(Maps.newHashMap(oligos.get(j)));
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
            int start = pos * smalligo;
            retOligos.add(new Oligo(protein.subList(start, start + oligoLength), oligo));
        }

        this.oligos.put(pos, retOligos);
    }

    public FragmentIterator fragmentIterator() {
        return new FragmentIterator(designs, oligos, protein, oligoLength, overlapLength, size, codonFrequencies);
    }

    public void fillFragments() {
        FragmentIterator it = fragmentIterator();
        while (it.hasNext()) {
            it.next().fill();
        }
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
                    for (int j = 0; j < size; j++) {
                        overlaps.add(Maps.newHashMap(overlaps.get(j)));
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

    public OverlapIterator overlapIterator() {
        checkState(overlaps != null);
        return new OverlapIterator(overlaps);
    }

    public void makeOverlapsUnique() {
        OverlapIterator it = overlapIterator();
        // position      index    swaps
        Map<Integer, Map<Integer, List<Integer>>> potentialSwaps = findPotentialSwaps();

        // position      indices
        Map<Integer, List<Integer>> permutations = Maps.newHashMap();
        for (int pos = 0; pos < size; pos++) {
            List<Integer> list = Lists.newArrayList(Collections.nCopies(overlapLength, 0));
            permutations.put(pos, list);
        }

        List<Overlap> visited = new ArrayList<>();
        int maxIndex = overlapLength - 1;
        while (it.hasNext()) {
            Overlap overlap = it.next();
            int pos = it.getCurrentPosition(); // must be after .next() call!!!

            // while current overlap is not unique (to all previous overlaps), start swapping
            while (matchesAnyVisited(overlap, visited)) {
                // ensure we haven't overflowed possible permutation indices
                checkState(permutations.get(pos).get(maxIndex) < potentialSwaps.get(pos).get(maxIndex).size());

                // swap all possible indices of current overlap based on permutation indices
                for (int overlapIndex = 0; overlapIndex < overlapLength; overlapIndex++) {
                    int permutationIndex = permutations.get(pos).get(overlapIndex);
                    List<Integer> swaps = potentialSwaps.get(pos).get(overlapIndex);
                    if (swaps == null || swaps.isEmpty()) continue;
                    int preAttachIndex = swaps.get(permutationIndex);
                    overlap.swapWithPreAttachments(overlapIndex, preAttachIndex);
                }

                // increment permutation indices
                int overlapPos = 0;
                List<Integer> perm = permutations.get(0);
                perm.set(overlapPos, perm.get(overlapPos) + 1);
                while (perm.get(overlapPos) >= potentialSwaps.get(pos).get(overlapPos).size()) {
                    perm.set(overlapPos, 0);
                    overlapPos++;
                    perm.set(overlapPos, perm.get(overlapPos) + 1);
                    if (overlapPos >= perm.size()) {
                        break;
                    }
                }
            }

            // add overlap to visted
            visited.add(overlap);
        }
    }

    private boolean matchesAnyVisited(Overlap overlap, List<Overlap> visited) {
        for (Overlap prev : visited) {
            if (Sequence.regionsMatch(overlap, prev)) {
                return true;
            }
        }
        return false;
    }

    private Map<Integer, Map<Integer, List<Integer>>> findPotentialSwaps() {
        Map<Integer, Map<Integer, List<Integer>>> map = Maps.newHashMap();
        for (int overlapPos = 0; overlapPos < size - 1; overlapPos++) {
            map.put(overlapPos, findPotentialSwaps(overlapPos));
        }
        return map;
    }

    private Map<Integer, List<Integer>> findPotentialSwaps(int pos) {
        Map<Integer, List<Integer>> overlapToCodons = Maps.newHashMap();
        for (int overlapIndex = 0; overlapIndex < overlapLength; overlapIndex++) {
            int start = pos * smalligo;
            for (int codonIndex = 0; codonIndex < smalligo; codonIndex++) {
                int globalCodonIndex = codonIndex + start;
                int globalOverlapIndex = start + smalligo + overlapIndex;
                Codon target = protein.get(globalCodonIndex);
                Codon current = protein.get(globalOverlapIndex);
                if (target != current && target.getAminoAcid() == current.getAminoAcid()) {
                    boolean important = false;
                    for (Overlap overlap : overlaps.get(pos)) {
                        Codon codonAtPos = overlap.get(overlapIndex);
                        if (codonsOfInterest.containsValue(codonAtPos)) {
                            important = true;
                            break;
                        }
                    }
                    for (Oligo oligo : oligos.get(pos)) {
                        Codon codonAtPos = oligo.get(codonIndex);
                        if (codonsOfInterest.containsValue(codonAtPos)) {
                            important = true;
                            break;
                        }
                    }
                    if (important) {
                        continue;
                    }
                    if (overlapToCodons.containsKey(overlapIndex)) {
                        overlapToCodons.get(overlapIndex).add(codonIndex);
                    } else {
                        List<Integer> list = Lists.newArrayList();
                        list.add(codonIndex);
                        overlapToCodons.put(overlapIndex, list);
                    }
                }

            }
        }
        return overlapToCodons;
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