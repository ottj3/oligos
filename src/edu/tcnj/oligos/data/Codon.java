package edu.tcnj.oligos.data;

public enum Codon {
    AAA("AAA"),
    ACA("ACA"),
    ATA("ATA"),
    AGA("AGA"),
    CAA("CAA"),
    CCA("CCA"),
    CTA("CTA"),
    CGA("CGA"),
    TAA("TAA"),
    TCA("TCA"),
    TTA("TTA"),
    TGA("TGA"),
    GAA("GAA"),
    GCA("GCA"),
    GTA("GTA"),
    GGA("GGA"),
    AAC("AAC"),
    ACC("ACC"),
    ATC("ATC"),
    AGC("AGC"),
    CAC("CAC"),
    CCC("CCC"),
    CTC("CTC"),
    CGC("CGC"),
    TAC("TAC"),
    TCC("TCC"),
    TTC("TTC"),
    TGC("TGC"),
    GAC("GAC"),
    GCC("GCC"),
    GTC("GTC"),
    GGC("GGC"),
    AAT("AAT"),
    ACT("ACT"),
    ATT("ATT"),
    AGT("AGT"),
    CAT("CAT"),
    CCT("CCT"),
    CTT("CTT"),
    CGT("CGT"),
    TAT("TAT"),
    TCT("TCT"),
    TTT("TTT"),
    TGT("TGT"),
    GAT("GAT"),
    GCT("GCT"),
    GTT("GTT"),
    GGT("GGT"),
    AAG("AAG"),
    ACG("ACG"),
    ATG("ATG"),
    AGG("AGG"),
    CAG("CAG"),
    CCG("CCG"),
    CTG("CTG"),
    CGG("CGG"),
    TAG("TAG"),
    TCG("TCG"),
    TTG("TTG"),
    TGG("TGG"),
    GAG("GAG"),
    GCG("GCG"),
    GTG("GTG"),
    GGG("GGG");

    private String name;

    Codon(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }
}
