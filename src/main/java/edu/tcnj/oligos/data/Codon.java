package edu.tcnj.oligos.data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
    TAG("TAG", AminoAcid.STOP),

    ALA("???", AminoAcid.ALA),
    ARG("???", AminoAcid.ARG),
    ASN("???", AminoAcid.ASN),
    ASP("???", AminoAcid.ASP),
    CYS("???", AminoAcid.CYS),
    GLN("???", AminoAcid.GLN),
    GLU("???", AminoAcid.GLU),
    GLY("???", AminoAcid.GLY),
    HIS("???", AminoAcid.HIS),
    LEU("???", AminoAcid.LEU),
    LYS("???", AminoAcid.LYS),
    ILE("???", AminoAcid.ILE),
    MET("???", AminoAcid.MET),
    PHE("???", AminoAcid.PHE),
    PRO("???", AminoAcid.PRO),
    SER("???", AminoAcid.SER),
    THR("???", AminoAcid.THR),
    TRP("???", AminoAcid.TRP),
    TYR("???", AminoAcid.TYR),
    VAL("???", AminoAcid.VAL),
    STOP("???", AminoAcid.STOP),

    PAD("???", AminoAcid.PAD);

    private String bases;
    private AminoAcid aa;
    private static Map<AminoAcid, List<Codon>> acidCodonMap = new HashMap<>();
    static {
        for (Codon c : Codon.values()) {
            if (acidCodonMap.containsKey(c.getAminoAcid())) {
                acidCodonMap.get(c.getAminoAcid()).add(c);
            } else {
                List<Codon> list = new ArrayList<>();
                list.add(c);
                acidCodonMap.put(c.getAminoAcid(), list);
            }
        }
    }

    Codon(String bases, AminoAcid aa) {
        this.bases = bases;
        this.aa = aa;
    }

    public String getBases() {
        return bases;
    }

    public AminoAcid getAminoAcid() {
        return aa;
    }

    public static List<Codon> getCodonsForAcid(AminoAcid acid) {
        return Collections.unmodifiableList(acidCodonMap.get(acid));
    }

    @Override
    public String toString() {
        return bases.equals("???") ? aa.getName() : getBases();
    }
}
