package edu.tcnj.oligos.data;

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
    STOP(".", "Stop");

    public String getCh() {
        return ch;
    }

    public String getName() {
        return name;
    }

    private String ch;
    private String name;

    AminoAcid(String ch, String name) {
        this.ch = ch;
        this.name = name;
    }
}
