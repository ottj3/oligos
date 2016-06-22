package edu.tcnj.oligos.library;

import com.google.common.base.Strings;
import com.google.common.collect.EnumBiMap;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Codon;

import java.util.Map;

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

    private Library(Protein protein, int size, int oligoLength, int overlapLength, Map<Codon, Design> designs, EnumBiMap<AminoAcid, Codon> codonsOfInterest) {
        this.protein = protein;
        this.size = size;
        this.oligoLength = oligoLength;
        this.overlapLength = overlapLength;
        this.designs = designs;
        this.codonsOfInterest = codonsOfInterest;
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