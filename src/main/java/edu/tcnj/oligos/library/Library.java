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

/**
 * Stores all information needed for building the oligos needed to create all gene variants desired.
 */
public class Library {

    private Protein protein;

    private final int size; // number of positions
    private final int oligoLength; // size of an oligo
    private final int overlapLength; // size of an overlap
    private final int smalligo; // oligoLength - overlapLength; value needed in many calculations

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
        //Calculate the relative frequencies of all codons of the same acid as the codons of interest.
        //Used for filling in the rest of the codons once the codons of interest are set.

        //Use Multisets to count the occurrences of each acid
        Map<AminoAcid, Multiset<Codon>> counts = Maps.newHashMap();
        for (Codon codon : protein) {
            AminoAcid acid = codon.getAminoAcid();
            if (codonsOfInterest.containsKey(acid) && !codonsOfInterest.containsValue(codon)) {
                Multiset<Codon> codons = counts.get(acid);
                if (codons == null) {
                    codons = HashMultiset.create();
                    counts.put(acid, codons);
                }
                codons.add(codon);
            }
        }

        //Based on the multiset, calculate the relative percentage for each codon for each acid
        Map<AminoAcid, Map<Codon, Double>> frequencies = Maps.newHashMap();
        for (Map.Entry<AminoAcid, Multiset<Codon>> entry : counts.entrySet()) {
            AminoAcid acid = entry.getKey();
            frequencies.put(acid, Maps.<Codon, Double>newHashMap());
            Multiset<Codon> codons = entry.getValue();
            int total = codons.size();
            for (Codon possible : codons.elementSet()) {
                int number = codons.count(possible);
                double freq = number / (double) total;
                frequencies.get(acid).put(possible, freq);
            }
        }
        return frequencies;
    }

    //Set the codons of interest frequencies as desired for the design,
    //and fill in all other codons of the same acid based on their frequencies
    private void setBaseFrequencies(Protein protein, Map<Codon, Double> minFreq) {
        List<Codon> sequence = Lists.newArrayList(protein.getSequence());
        //For every codon of interest
        for (Map.Entry<Codon, Double> entry : minFreq.entrySet()) {
            Codon codonOfInterest = entry.getKey();
            AminoAcid acidOfInterest = codonOfInterest.getAminoAcid();
            Design design = designs.get(codonOfInterest);
            double baseFreq = entry.getValue();

            //Get every occurrence of a codon of this amino acid.
            //Used to figure out the number of codons needed to hit the base percentage,
            //as well as for filling in other codons after the codon of interest is set
            List<Integer> allCodonSpots = Lists.newArrayList();
            for (int i = 0; i < sequence.size(); i++) {
                if (sequence.get(i).getAminoAcid() == acidOfInterest) allCodonSpots.add(i);
            }

            //Find potential spots to put the codon of interest
            //(the codon can go anywhere not in a range where that codon is controlled)
            List<Integer> potentialCOIspots = Lists.newArrayList();
            for (int pos = 0; pos < size; pos++) {
                if (design.getDeltasForPosition(pos) == null) {
                    //If this position doesn't have a delta for it, than any spot in
                    //this position (excluding the overlaps) can be a potential spot
                    int start = pos == 0 ? 0 : pos * smalligo + overlapLength;
                    for (int i = start; i < (pos + 1) * smalligo; i++) {
                        if (protein.get(i).getAminoAcid() == acidOfInterest) {
                            potentialCOIspots.add(i);
                        }
                    }
                }
                if (pos == size - 1 || design.getDeltasForRange(new Fragment.Range(pos, pos + 1)) == null) {
                    //For the overlap following this position: if the overlap isn't part of a
                    //controlled range for this codon (i.e. when both the spot preceding and succeeding
                    //it are from the same fragment) then spots in it can be used for the base frequency
                    int start = (pos + 1) * smalligo;
                    for (int i = start; i < start + overlapLength; i++) {
                        if (protein.get(i).getAminoAcid() == acidOfInterest) {
                            potentialCOIspots.add(i);
                        }
                    }
                }
            }
            //Shuffle the positions to fill them in in random order
            Collections.shuffle(potentialCOIspots);
            //Calculate the base integer number of occurrences that is closest to the base percentage
            int baseOccurrences = (int) (baseFreq * allCodonSpots.size() + 0.5);
            //Trim the potential spots to only have as many spots as needed
            //(necessary for filtering out the used positions from allCodonSpots later)
            potentialCOIspots = potentialCOIspots.subList(0, baseOccurrences);
            for (Integer potentialSpot : potentialCOIspots) {
                sequence.set(potentialSpot, codonOfInterest);
            }

            //Fill in all the spots not used for the codon of interest based on their original frequencies
            allCodonSpots.removeAll(potentialCOIspots);
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
    }

    /**
     * get the best approximate counts to fill the spots with the given relative frequencies
     *
     * @param frequencies the relative frequencies of all codons to be used
     * @param totalCount  the number of spots to be filled
     * @return a map from codon to the number of times it should appear (counts will sum to totalCount)
     */
    public static Map<Codon, Integer> findCodonCounts(Map<Codon, Double> frequencies, int totalCount) {
        Map<Codon, Integer> counts = Maps.newHashMap();
        int totalSoFar = 0;
        //Get approximate counts for every codon, rounding down always
        for (Map.Entry<Codon, Double> entry : frequencies.entrySet()) {
            int thisCount = (int) (entry.getValue() * totalCount);
            totalSoFar += thisCount;
            counts.put(entry.getKey(), thisCount);
        }
        //Hit the total needed by continuously adding one to the count which
        //is furthest under its expected (exact, non-integral) count
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

    /**
     * Create oligos for every position, figuring out how many are needed based on the design.
     * Does not fill them based on the design; all oligos will be copies from the original rna sequence
     */
    public void createOligos() {
        if (oligos != null) return;
        oligos = Maps.newHashMap();
        for (int i = 0; i < size; i++) {
            createOligosForPosition(i);
        }
    }

    private void createOligosForPosition(int pos) {
        //for a given position (oligo # in overall sequence), create as many oligos as necessary to make the design
        checkArgument(0 <= pos && pos < size);
        List<Map<Codon, Integer>> oligos = Lists.newArrayList();
        //Add a blank map; to be cloned/filled later
        oligos.add(Maps.<Codon, Integer>newHashMap());
        for (Map.Entry<Codon, Design> entry : designs.entrySet()) {
            Design design = entry.getValue();
            List<Integer> deltasForPos = design.getDeltasForPosition(pos);
            //If a design covers this position, it needs oligos for each delta level needed.
            //The total number of oligos needed in a position is the product of the number of
            //levels needed for each design at this position
            if (deltasForPos != null) {
                int size = oligos.size();
                for (int i = 1; i < deltasForPos.size(); i++) {
                    for (int j = 0; j < size; j++) {
                        //Make as many copies as needed to add this design
                        oligos.add(Maps.newHashMap(oligos.get(j)));
                    }
                }
                //Fill in the originals and the copies with different levels
                for (int i = 0; i < deltasForPos.size(); i++) {
                    for (int j = i * size; j < (i + 1) * size; j++) {
                        Map<Codon, Integer> oligo = oligos.get(j);
                        oligo.put(entry.getKey(), deltasForPos.get(i));
                    }
                }
            }
        }
        //Given the maps of deltas, make oligos that have both
        //the proper section of the original RNA as well as these deltas
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

    /**
     * Create overlaps for every position, figuring out how many are needed based on the design.
     * Does not make the overlaps unique, but sets each to match with its corresponding oligos.
     */
    public void createOverlaps() {
        if (overlaps != null) return;
        overlaps = Maps.newHashMap();
        for (int i = 0; i < size - 1; i++) {
            createOverlapsForPosition(i);
        }
    }

    private void createOverlapsForPosition(int pos) {
        //Same logic as the createOligosForPosition method. Only change is that for an
        //overlap to be made, a design has to span both of its bordering positions.
        //Overlap positions are indexed by the oligo position before it
        //(i.e. overlap 2 would connect oligos from positions 2 and 3).
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
        //Link an ovelap to all pre- and post- attaching oligos that it corresponds to.
        //for an oligo to match an overlap, it must share all the same codon->delta pairs
        //that the overlap has (although the oligo can have more pairs than the overlap)
        List<Overlap> retOverlaps = Lists.newArrayList();
        for (Map<Codon, Integer> deltas : overlaps) {
            Overlap overlap = new Overlap(null, deltas);
            //For pre attachments: if an oligo in this position matches all key/value pairs, link it to the overlap
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
            //For post attachments
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
            //Set the overlap's content based on its first pre attachment, and make sure
            //that content matches in all of its attached oligos
            overlap.finalizeLink(overlapLength);
            retOverlaps.add(overlap);
        }

        this.overlaps.put(pos, retOverlaps);
    }

    public OverlapIterator overlapIterator() {
        checkState(overlaps != null);
        return new OverlapIterator(overlaps);
    }

    /**
     * Change overlaps by swapping codons of from the same acid so that no two overlaps are the same.
     * This allows the oligos to connect only when they are supposed to.
     */
    public void makeOverlapsUnique() {
        OverlapIterator it = overlapIterator();
        // position      index    swaps
        Map<Integer, Map<Integer, List<Integer>>> potentialSwaps = findPotentialSwaps();

        // position      indices
        Map<Integer, List<Integer>> permutations = Maps.newHashMap();
        for (int pos = 0; pos < size - 1; pos++) {
            //Initialize the permutation indices to -1.
            //A -1 in part of the overlap means don't swap anything.
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

                // increment permutation indices (to go through all possible permutations)
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

            // add overlap to visited
            visited.add(overlap);
        }
    }

    private void doOverlapPermutation(Overlap overlap, List<Integer> permutationIndices,
                                      Map<Integer, List<Integer>> potentialSwap) {
        //Based on the given permutation indices, swap each part of the oligo
        for (int overlapIndex = 0; overlapIndex < overlapLength; overlapIndex++) {
            int permutationIndex = permutationIndices.get(overlapIndex);

            //-1 means "do nothing" so skip swaps for the given codon
            if (permutationIndex == -1) continue;

            //From the potentialSwaps list, get the swap corresponding to this index
            List<Integer> swaps = potentialSwap.get(overlapIndex);
            if (swaps == null || swaps.isEmpty()) continue;
            int attachIndex = swaps.get(permutationIndex);

            //Swap with either the pre or post attachments based on the index
            if (attachIndex < oligoLength) {
                overlap.swapWithAttachments(overlapIndex, attachIndex, overlap.getPreAttachments());
            } else {
                overlap.swapWithAttachments(overlapIndex, attachIndex - oligoLength, overlap.getPostAttachments());
            }
        }
    }

    private boolean matchesAnyVisited(Overlap overlap, List<Overlap> visited) {
        //Make sure that this overlap is unique to all the others made so far
        for (Overlap prev : visited) {
            if (Sequence.regionsMatch(overlap, prev)) {
                return true;
            }
        }
        return false;
    }

    private Map<Integer, Map<Integer, List<Integer>>> findPotentialSwaps() {
        //Find all possible swaps for all overlap positions
        Map<Integer, Map<Integer, List<Integer>>> map = Maps.newHashMap();
        for (int overlapPos = 0; overlapPos < size - 1; overlapPos++) {
            map.put(overlapPos, findPotentialSwaps(overlapPos));
        }
        return map;
    }

    private Map<Integer, List<Integer>> findPotentialSwaps(int pos) {
        //Find potential swaps for a given overlap position. A swap can happen when the positions share the same amino acid,
        //and no oligo has a controlled codon in that position.
        Map<Integer, List<Integer>> overlapToCodons = Maps.newHashMap();
        //For every codon in the overlap region, find potential swaps
        for (int overlapIndex = 0; overlapIndex < overlapLength; overlapIndex++) {
            //Offest to check amino acid (start is the position that the preattachments start in the overall sequence)
            int start = pos * smalligo;
            overlapToCodons.put(overlapIndex, Lists.<Integer>newArrayList());
            int codonIndexStart = pos == 0 ? 0 : overlapLength;
            //For every codon in the preceding oligos
            for (int codonIndex = codonIndexStart; codonIndex < smalligo; codonIndex++) {
                int globalCodonIndex = codonIndex + start;
                int globalOverlapIndex = start + smalligo + overlapIndex;
                Codon target = protein.get(globalCodonIndex);
                Codon current = protein.get(globalOverlapIndex);
                //If the amino acids match, this could be a swap
                if (target != current && target.getAminoAcid() == current.getAminoAcid()) {
                    boolean important = false;
                    for (Overlap overlap : overlaps.get(pos)) {
//                        Codon codonInOverlap = overlap.get(overlapIndex);
//                        if (codonsOfInterest.containsValue(codonInOverlap)) {
//                            important = true;
//                            break;
//                        }
                        //If any oligo in this position has a codon of interest at this codon, it can't be swapped.
                        for (Oligo oligo : overlap.getPreAttachments()) {
                            Codon codonAtPos = oligo.get(codonIndex);
                            if (codonsOfInterest.containsValue(codonAtPos) && !overlap.getDeltas().containsKey(codonAtPos)) {
                                important = true;
                                break;
                            }
                        }
                    }

                    if (!important) {
                        overlapToCodons.get(overlapIndex).add(codonIndex);
                    }
                }
            }
            //Same logic as above, for post attachments. When a potential swap is found in post attachments, its
            //index is increased by oligoLength to differentiate it from pre attachments.
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

                    if (!important) {
                        overlapToCodons.get(overlapIndex).add(codonIndex + oligoLength);
                    }
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
        private int seqStart = Integer.MIN_VALUE;
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
            Protein protein = (seqStart == Integer.MIN_VALUE)
                    ? new Protein(new Sequence(proteinRNA))
                    : new Protein(new Sequence(proteinRNA), seqStart, seqEnd);
            int size = ((seqEnd - seqStart) - overlapSize) / (oligoLength - overlapSize);

            return new Library(protein, size, oligoLength, overlapSize, designs, codonsOfInterest, minFrequencies);
        }
    }

}