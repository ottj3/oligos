package edu.tcnj.oligos.library;

import com.google.common.base.Joiner;
import com.google.common.base.MoreObjects;
import com.google.common.collect.Lists;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Codon;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * A special sequence to store both a list of codons, and the
 * corresponding list of amino acids that make up the protein
 */
public class Protein extends Sequence {
    private List<AminoAcid> aaSeq;

    public Protein(Sequence sequence) {
        this(sequence, 0, sequence.size());
    }

    /**
     * Construct a protein from a subsequence of codons. Pads the codon
     * sequence and the amino acid sequence according to start and end.
     *
     * @param sequence the codon sequence used to construct this protein
     * @param start    the index of the first codon (offset from the start of sequence), can be positive or negative
     * @param end      the index (exclusive) of the last codon from the sequence to be used. Can be greater than seq.size()
     */
    public Protein(Sequence sequence, int start, int end) {
        List<AminoAcid> tempSeq = new ArrayList<>(end - start);
        List<Codon> tempCodons = new ArrayList<>(end - start);
        while (start < 0) {
            tempSeq.add(AminoAcid.PAD);
            tempCodons.add(Codon.PAD);
            start++;
        }
        for (int i = 0; i < end; i++) {
            if (i < start || i >= sequence.size()) {
                tempSeq.add(AminoAcid.PAD);
                tempCodons.add(Codon.PAD);
            } else {
                Codon c = sequence.get(i - start);
                tempCodons.add(c);
                tempSeq.add(c.getAminoAcid());
            }
        }
        this.aaSeq = Collections.unmodifiableList(tempSeq);
        this.setSequence(Collections.unmodifiableList(tempCodons));
    }

    /**
     * Constructs a protein based on a string of bases
     *
     * @param string a string of RNA to be used to fill the codon sequence and then the amino acid sequence
     */
    public Protein(String string) {
        super(string);
        this.aaSeq = Lists.newArrayList();
        for (Codon codon : this.sequence) {
            this.aaSeq.add(codon.getAminoAcid());
        }
        this.aaSeq = Collections.unmodifiableList(aaSeq);
    }

    public List<AminoAcid> getAminoAcidSequence() {
        return aaSeq;
    }

    @Override
    public String toString() {
        return MoreObjects.toStringHelper(this)
                .add("aaSeq", Joiner.on("").join(aaSeq))
                .add("seq", Joiner.on("").join(sequence))
                .toString();
    }
}
