package edu.tcnj.oligos.library;

import com.google.common.base.Strings;
import com.google.common.collect.EnumBiMap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
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
                    EnumBiMap<AminoAcid, Codon> codonsOfInterest, Map<Codon, Double> minFreq) {
        this.protein = protein;
        this.size = size;
        this.oligoLength = oligoLength;
        this.overlapLength = overlapLength;
        this.smalligo = oligoLength - overlapLength;
        this.designs = designs;
        this.codonsOfInterest = codonsOfInterest;

        this.codonFrequencies = calcFrequencies();
        setBaseFrequencies(this.protein, minFreq);
    }

    private Map<AminoAcid, Map<Codon, Double>> calcFrequencies() {
        Map<AminoAcid, Multiset<Codon>> counts = Maps.newHashMap();
        for (Codon codon : protein) {
            AminoAcid acid = codon.getAminoAcid();
            if (codonsOfInterest.containsKey(acid) && !codonsOfInterest.containsValue(codon)) {
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
                if (number == 0) continue;
                double freq = number / (double) total;
                frequencies.get(acid).put(possible, freq);
            }
        }
        return frequencies;
    }

    private Protein setBaseFrequencies(Protein protein, Map<Codon, Double> minFreq) {
        List<Codon> sequence = Lists.newArrayList(protein.getSequence());
        for (Map.Entry<Codon, Double> entry : minFreq.entrySet()) {
            Codon codonOfInterest = entry.getKey();
            AminoAcid acidOfInterest = codonOfInterest.getAminoAcid();
            Design design = designs.get(codonOfInterest);
            double baseFreq = entry.getValue();

            List<Integer> allCodonSpots = Lists.newArrayList();
            for (int i = 0; i < sequence.size(); i++) {
                if (sequence.get(i).getAminoAcid() == acidOfInterest) allCodonSpots.add(i);
            }

            List<Integer> potentialBaseCodonSpots = Lists.newArrayList();
            for (int pos = 0; pos < size; pos++) {
                if (design.getDeltasForPosition(pos) == null) {
                    int start = pos == 0 ? 0 : pos * smalligo + overlapLength;
                    for (int i = start; i < (pos + 1) * smalligo; i++) {
                        if (protein.get(i).getAminoAcid() == acidOfInterest) {
                            potentialBaseCodonSpots.add(i);
                        }
                    }
                }
                if (pos == size - 1 || design.getDeltasForRange(new Fragment.Range(pos, pos + 1)) == null) {
                    int start = (pos + 1) * smalligo;
                    for (int i = start; i < start + overlapLength; i++) {
                        if (protein.get(i).getAminoAcid() == acidOfInterest) {
                            potentialBaseCodonSpots.add(i);
                        }
                    }
                }
            }
            Collections.shuffle(potentialBaseCodonSpots);
            int baseOccurrences = (int) (baseFreq * allCodonSpots.size() + 0.5);
            potentialBaseCodonSpots = potentialBaseCodonSpots.subList(0, baseOccurrences);
            for (int i = 0; i < baseOccurrences; i++) {
                sequence.set(potentialBaseCodonSpots.get(i), codonOfInterest);
            }


            allCodonSpots.removeAll(potentialBaseCodonSpots);
            Collections.shuffle(allCodonSpots);
            Map<Codon, Integer> counts = findCodonCounts(codonFrequencies.get(acidOfInterest),
                    allCodonSpots.size());
            int index = 0;
            for (Map.Entry<Codon, Integer> count : counts.entrySet()) {
                for (int j = 0; j < count.getValue(); j++) {
                    sequence.set(allCodonSpots.get(index), count.getKey());
                    index++;
                }
            }
        }
        protein.setSequence(sequence);
        return protein;
    }

    public static Map<Codon, Integer> findCodonCounts(Map<Codon, Double> frequencies, int totalCount) {
        Map<Codon, Integer> counts = Maps.newHashMap();
        int totalSoFar = 0;
        for (Map.Entry<Codon, Double> entry : frequencies.entrySet()) {
            int thisCount = (int) (entry.getValue() * totalCount);
            totalSoFar += thisCount;
            counts.put(entry.getKey(), thisCount);
        }
        while (totalSoFar < totalCount) {
            double maxDelta = Double.NEGATIVE_INFINITY;
            Codon codonToIncrease = Codon.PAD;
            for (Map.Entry<Codon, Double> entry : frequencies.entrySet()) {
                double thisDelta = entry.getValue() * totalCount - counts.get(entry.getKey());
                if (thisDelta > maxDelta) {
                    maxDelta = thisDelta;
                    codonToIncrease = entry.getKey();
                }
            }
            counts.put(codonToIncrease, counts.get(codonToIncrease) + 1);
            totalSoFar++;
        }
        return counts;
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

    public void createOverlaps() {
        if (overlaps != null) return;
        overlaps = Maps.newHashMap();
        for (int i = 0; i < size - 1; i++) {
            createOverlapsForPosition(i);
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
                boolean matches = true;
                for (Map.Entry<Codon, Integer> entry : deltas.entrySet()) {
                    if (!Objects.equals(entry.getValue(), oligo.getDeltas().get(entry.getKey()))) {
                        matches = false;
                        break;
                    }
                }
                if (matches) {
                    overlap.linkPreAttachment(oligo);
                }
            }
            for (Oligo oligo : oligos.get(pos + 1)) {
                boolean matches = true;
                for (Map.Entry<Codon, Integer> entry : deltas.entrySet()) {
                    if (!Objects.equals(entry.getValue(), oligo.getDeltas().get(entry.getKey()))) {
                        matches = false;
                        break;
                    }
                }
                if (matches) {
                    overlap.linkPostAttachment(oligo);
                }
            }
            overlap.finalizeLink(overlapLength);
            retOverlaps.add(overlap);
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
        for (int pos = 0; pos < size - 1; pos++) {
            List<Integer> list = Lists.newArrayList(Collections.nCopies(overlapLength, -1));
            permutations.put(pos, list);
        }

        List<Overlap> visited = new ArrayList<>();
        int maxIndex = overlapLength - 1;
        while (it.hasNext()) {
            Overlap overlap = it.next();
            int pos = it.getCurrentPosition(); // must be after .next() call!!!
            boolean matches = matchesAnyVisited(overlap, visited);
            // while current overlap is not unique (to all previous overlaps), start swapping
            while (matches) {
                // ensure we haven't overflowed possible permutation indices
                checkState(permutations.get(pos).get(maxIndex) < potentialSwaps.get(pos).get(maxIndex).size());

                // swap all possible indices of current overlap based on permutation indices
                doOverlapPermutation(overlap, permutations.get(pos), potentialSwaps.get(pos));

                matches = matchesAnyVisited(overlap, visited);

                if (matches) {
                    //Undo the swap if it is not unique (by swapping again with the same permutation indices)
                    doOverlapPermutation(overlap, permutations.get(pos), potentialSwaps.get(pos));
                }

                // increment permutation indices
                int overlapPos = 0;
                List<Integer> perm = permutations.get(pos);
                perm.set(overlapPos, perm.get(overlapPos) + 1);
                while (perm.get(overlapPos) >= potentialSwaps.get(pos).get(overlapPos).size()
                        && overlapPos + 1 < perm.size()) {
                    perm.set(overlapPos, -1);
                    overlapPos++;
//                    if (overlapPos >= perm.size()) {
//                        break;
//                    }
                    perm.set(overlapPos, perm.get(overlapPos) + 1);
                }
            }

            // add overlap to visted
            visited.add(overlap);
        }
    }

    private void doOverlapPermutation(Overlap overlap, List<Integer> permutationIndices,
                                      Map<Integer, List<Integer>> potentialSwap) {
        for (int overlapIndex = 0; overlapIndex < overlapLength; overlapIndex++) {
            int permutationIndex = permutationIndices.get(overlapIndex);
            if (permutationIndex == -1) continue;
            List<Integer> swaps = potentialSwap.get(overlapIndex);
            if (swaps == null || swaps.isEmpty()) continue;
            int attachIndex = swaps.get(permutationIndex);

            if (attachIndex < oligoLength) {
                overlap.swapWithAttachments(overlapIndex, attachIndex, overlap.getPreAttachments());
            } else {
                overlap.swapWithAttachments(overlapIndex, attachIndex - oligoLength, overlap.getPostAttachments());
            }
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
            overlapToCodons.put(overlapIndex, Lists.<Integer>newArrayList());
            int codonIndexStart = pos == 0 ? 0 : overlapLength;
            for (int codonIndex = codonIndexStart; codonIndex < smalligo; codonIndex++) {
                int globalCodonIndex = codonIndex + start;
                int globalOverlapIndex = start + smalligo + overlapIndex;
                Codon target = protein.get(globalCodonIndex);
                Codon current = protein.get(globalOverlapIndex);
                if (target != current && target.getAminoAcid() == current.getAminoAcid()) {
                    boolean important = false;
                    for (Overlap overlap : overlaps.get(pos)) {
//                        Codon codonInOverlap = overlap.get(overlapIndex);
//                        if (codonsOfInterest.containsValue(codonInOverlap)) {
//                            important = true;
//                            break;
//                        }
                        for (Oligo oligo : overlap.getPreAttachments()) {
                            Codon codonAtPos = oligo.get(codonIndex);
                            if (codonsOfInterest.containsValue(codonAtPos) && !overlap.getDeltas().containsKey(codonAtPos)) {
                                important = true;
                                break;
                            }
                        }
                    }

                    if (important) {
                        continue;
                    }

                    overlapToCodons.get(overlapIndex).add(codonIndex);
                }

            }
            //Post attachments
            int codonIndexEnd = pos == size - 1 ? oligoLength : smalligo;
            for (int codonIndex = overlapLength; codonIndex < codonIndexEnd; codonIndex++) {
                int globalCodonIndex = codonIndex + start + smalligo;
                int globalOverlapIndex = start + smalligo + overlapIndex;
                Codon target = protein.get(globalCodonIndex);
                Codon current = protein.get(globalOverlapIndex);
                if (target != current && target.getAminoAcid() == current.getAminoAcid()) {
                    boolean important = false;
                    for (Overlap overlap : overlaps.get(pos)) {
//                        Codon codonInOverlap = overlap.get(overlapIndex);
//                        if (codonsOfInterest.containsValue(codonInOverlap)) {
//                            important = true;
//                            break;
//                        }
                        for (Oligo oligo : overlap.getPostAttachments()) {
                            Codon codonAtPos = oligo.get(codonIndex);
                            if (codonsOfInterest.containsValue(codonAtPos) && !overlap.getDeltas().containsKey(codonAtPos)) {
                                important = true;
                                break;
                            }
                        }
                    }

                    if (important) {
                        continue;
                    }

                    overlapToCodons.get(overlapIndex).add(codonIndex + oligoLength);
                }

            }
        }
        return overlapToCodons;
    }

    public Map<Integer, List<Oligo>> getOligos() {
        return oligos;
    }

    public static class Builder {
        private String proteinRNA = "";
        private int seqStart = -1;
        private int seqEnd;
        private int oligoLength = -1;
        private int overlapSize = -1;
        private Map<Codon, Design> designs;
        private EnumBiMap<AminoAcid, Codon> codonsOfInterest;
        private Map<Codon, Double> minFrequencies;

        public Builder withSequenceLength(int start, int end) {
            checkArgument(start < end,
                    "Invalid start and end positions for sequence: %s, %s", start, end);
            this.seqStart = start;
            this.seqEnd = end;
            return this;
        }

        public Builder withProteinFromRNA(String rna) {
            checkArgument(!Strings.isNullOrEmpty(rna), "Empty input RNA.");
            checkArgument(rna.length() % 3 == 0, "Input RNA length not multiple of 3.");
            this.proteinRNA = rna;
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

        public Builder withMinFrequencies(Map<Codon, Double> frequencies) {
//            for (Map.Entry<Codon, Double> entry : frequencies.entrySet()) {
//                checkArgument(codonsOfInterest.containsValue(entry.getKey()),
//                        "Frequencies must be for codons of interest");
//            }
            this.minFrequencies = frequencies;
            return this;
        }

        public Library build() {
            checkState(!proteinRNA.isEmpty());
            checkState(designs != null);
            checkState(codonsOfInterest != null);
            checkState(oligoLength != -1);
            checkState(overlapSize != -1);
            checkState(minFrequencies != null);
            Protein protein = (seqStart == -1)
                    ? new Protein(new Sequence(proteinRNA))
                    : new Protein(new Sequence(proteinRNA), seqStart, seqEnd);
            int size = ((seqEnd - seqStart) - overlapSize) / (oligoLength - overlapSize);

            return new Library(protein, size, oligoLength, overlapSize, designs, codonsOfInterest, minFrequencies);
        }
    }

}