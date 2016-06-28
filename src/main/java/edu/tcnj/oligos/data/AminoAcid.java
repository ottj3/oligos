package edu.tcnj.oligos.data;

import java.util.HashMap;
import java.util.Map;

public enum AminoAcid {
    ALA("A", "Alanine"),
    ARG("R", "Arginine"),
    ASN("N", "Asparagine"),
    ASP("D", "Aspartic acid"),
    CYS("C", "Cysteine"),
    GLN("Q", "Glutamine"),
    GLU("E", "Glutamic acid"),
    GLY("G", "Glycine"),
    HIS("H", "Histidine"),
    LEU("L", "Leucine"),
    LYS("K", "Lysine"),
    ILE("I", "Isoleucine"),
    MET("M", "Methionine"),
    PHE("F", "Phenylalanine"),
    PRO("P", "Proline"),
    SER("S", "Serine"),
    THR("T", "Threonine"),
    TRP("W", "Tryptophan"),
    TYR("Y", "Tyrosine"),
    VAL("V", "Valine"),
    STOP("*", "Stop"),
    PAD("X", "Padding");

    AminoAcid(String ch, String name) {
        this.ch = ch;
        this.name = name;
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

    void setWildcard(Codon wildcard) {
        this.wildcard = wildcard;
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

    @Override
    public String toString() {
        return ch;
    }
}
