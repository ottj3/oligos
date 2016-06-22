package edu.tcnj.oligos.data;

import java.util.HashMap;
import java.util.Map;

public enum AminoAcid {
    ALA("A", "Alanine", Codon.ALA),
    ARG("R", "Arginine", Codon.ARG),
    ASN("N", "Asparagine", Codon.ASN),
    ASP("D", "Aspartic acid", Codon.ASP),
    CYS("C", "Cysteine", Codon.CYS),
    GLN("Q", "Glutamine", Codon.GLN),
    GLU("E", "Glutamic acid", Codon.GLU),
    GLY("G", "Glycine", Codon.GLY),
    HIS("H", "Histidine", Codon.HIS),
    LEU("L", "Leucine", Codon.LEU),
    LYS("K", "Lysine", Codon.LYS),
    ILE("I", "Isoleucine", Codon.ILE),
    MET("M", "Methionine", Codon.MET),
    PHE("F", "Phenylalanine", Codon.PHE),
    PRO("P", "Proline", Codon.PRO),
    SER("S", "Serine", Codon.SER),
    THR("T", "Threonine", Codon.THR),
    TRP("W", "Tryptophan", Codon.TRP),
    TYR("Y", "Tyrosine", Codon.TYR),
    VAL("V", "Valine", Codon.VAL),
    STOP(".", "Stop", Codon.STOP);

    AminoAcid(String ch, String name, Codon wildcard) {
        this.ch = ch;
        this.name = name;
        this.wildcard = wildcard;
    }

    private String ch;
    private String name;
    private Codon wildcard;
    private static Map<String, AminoAcid> seqMap = new HashMap<>();
    static {
        for (AminoAcid aa : AminoAcid.values()) {
            seqMap.put(aa.getCh(), aa);
        }
    }

    public String getCh() {
        return ch;
    }

    public String getName() {
        return name;
    }

    public Codon getWildcard() {
        return wildcard;
    }

    public static AminoAcid getAcidForSymbol(String character) {
        return seqMap.get(character);
    }
}
