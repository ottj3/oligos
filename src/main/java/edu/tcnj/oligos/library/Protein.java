package edu.tcnj.oligos.library;

import com.google.common.base.MoreObjects;
import edu.tcnj.oligos.data.AminoAcid;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Protein extends Sequence {
    private List<AminoAcid> aaSeq;

    public Protein(Sequence sequence) {
        this(sequence, 0, sequence.size());
    }

    public Protein(Sequence sequence, int start, int end) {
        List<AminoAcid> tempSeq = new ArrayList<>(sequence.size());
        for (int i = 0; i < sequence.size(); i++) {
            if (i < start || i > end) {
                tempSeq.add(AminoAcid.PAD);
            } else {
                tempSeq.add(sequence.get(i).getAminoAcid());
            }
        }
        this.aaSeq = Collections.unmodifiableList(tempSeq);
        this.setSequence(Collections.unmodifiableList(sequence));
    }

    /*
    // Uncomment these constructors and method if, for some reason, we have to
    // get an input amino acid sequence without a corresponding RNA sequence.

    public Protein(String string) {
        this(string, 0, string.length());
    }
    public Protein(String string, int start, int end) {
        List<AminoAcid> tempSeq = new ArrayList<>();
        char[] chars = string.toCharArray();
        for (int i = 0; i < chars.length; i++) {
            AminoAcid acid;
            if (i < start || i > end) {
                acid = AminoAcid.PAD;
            } else {
                acid = AminoAcid.getAcidForSymbol(String.valueOf(chars[i]));
            }
            tempSeq.add(acid);
        }
        this.aaSeq = Collections.unmodifiableList(tempSeq);
        this.setSequence(acidsToWildcardCodons(tempSeq));
    }
    private static List<Codon> acidsToWildcardCodons(List<AminoAcid> seq) {
        List<Codon> list = new ArrayList<>(seq.size() * 3);
        for (AminoAcid aminoAcid : seq) {
            list.add(aminoAcid.getWildcard());
        }
        return list;
    }
    */

    public List<AminoAcid> getAminoAcidSequence() {
        return aaSeq;
    }

    @Override
    public String toString() {
        return MoreObjects.toStringHelper(this)
                .add("aaSeq", aaSeq.toArray())
                .add("seq", getSequence().toArray())
                .toString();
    }
}
