package edu.tcnj.oligos.data;

public enum Codon {
    ATG("ATG", AminoAcid.MET),
    ACA("ACA", AminoAcid.THR),
    ATA("ATA", AminoAcid.ILE),
    AGA("AGA", AminoAcid.ARG),
    CAA("CAA", AminoAcid.GLN),
    CCA("CCA", AminoAcid.PRO),
    CTA("CTA", AminoAcid.LEU),
    CGA("CGA", AminoAcid.ARG),
    TCA("TCA", AminoAcid.SER),
    TTA("TTA", AminoAcid.LEU),
    GAA("GAA", AminoAcid.GLU),
    GCA("GCA", AminoAcid.ALA),
    GGA("GGA", AminoAcid.GLY),
    AAC("AAC", AminoAcid.ASN),
    ACC("ACC", AminoAcid.THR),
    ATC("ATC", AminoAcid.ILE),
    AGC("AGC", AminoAcid.SER),
    CAC("CAC", AminoAcid.HIS),
    CCC("CCC", AminoAcid.PRO),
    CTC("CTC", AminoAcid.LEU),
    CGC("CGC", AminoAcid.ARG),
    TAC("TAC", AminoAcid.TYR),
    TCC("TCC", AminoAcid.SER),
    TTC("TTC", AminoAcid.PHE),
    TGC("TGC", AminoAcid.CYS),
    GAC("GAC", AminoAcid.ASP),
    GCC("GCC", AminoAcid.ALA),
    GGC("GGC", AminoAcid.GLY),
    AAT("AAT", AminoAcid.ASN),
    ACT("ACT", AminoAcid.THR),
    ATT("ATT", AminoAcid.ILE),
    AGT("AGT", AminoAcid.SER),
    CAT("CAT", AminoAcid.HIS),
    CCT("CCT", AminoAcid.PRO),
    CTT("CTT", AminoAcid.LEU),
    CGT("CGT", AminoAcid.ARG),
    TAT("TAT", AminoAcid.TYR),
    TCT("TCT", AminoAcid.SER),
    TTT("TTT", AminoAcid.PHE),
    TGT("TGT", AminoAcid.CYS),
    GAT("GAT", AminoAcid.ASP),
    GCT("GCT", AminoAcid.ALA),
    GGT("GGT", AminoAcid.GLY),
    ACG("ACG", AminoAcid.THR),
    AGG("AGG", AminoAcid.ARG),
    CAG("CAG", AminoAcid.GLN),
    CCG("CCG", AminoAcid.PRO),
    CTG("CTG", AminoAcid.LEU),
    CGG("CGG", AminoAcid.ARG),
    TCG("TCG", AminoAcid.SER),
    TTG("TTG", AminoAcid.LEU),
    TGG("TGG", AminoAcid.TRP),
    GAG("GAG", AminoAcid.GLU),
    AAA("AAA", AminoAcid.LYS),
    AAG("AAG", AminoAcid.LYS),
    GCG("GCG", AminoAcid.ALA),
    GGG("GGG", AminoAcid.GLY),
    GTA("GTA", AminoAcid.VAL),
    GTT("GTT", AminoAcid.VAL),
    GTC("GTC", AminoAcid.VAL),
    GTG("GTG", AminoAcid.VAL),
    TAA("TAA", AminoAcid.STOP),
    TGA("TGA", AminoAcid.STOP),
    TAG("TAG", AminoAcid.STOP);

    private String name;
    private AminoAcid aa;

    Codon(String name, AminoAcid aa) {
        this.name = name;
        this.aa = aa;
    }

    public String getName() {
        return name;
    }

    public AminoAcid getAminoAcid() {
        return aa;
    }
}
