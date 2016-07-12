package edu.tcnj.oligos.data;

import com.google.common.collect.Collections2;
import com.google.common.collect.Lists;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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

    ALA("???", AminoAcid.ALA, true),
    ARG("???", AminoAcid.ARG, true),
    ASN("???", AminoAcid.ASN, true),
    ASP("???", AminoAcid.ASP, true),
    CYS("???", AminoAcid.CYS, true),
    GLN("???", AminoAcid.GLN, true),
    GLU("???", AminoAcid.GLU, true),
    GLY("???", AminoAcid.GLY, true),
    HIS("???", AminoAcid.HIS, true),
    LEU("???", AminoAcid.LEU, true),
    LYS("???", AminoAcid.LYS, true),
    ILE("???", AminoAcid.ILE, true),
    MET("???", AminoAcid.MET, true),
    PHE("???", AminoAcid.PHE, true),
    PRO("???", AminoAcid.PRO, true),
    SER("???", AminoAcid.SER, true),
    THR("???", AminoAcid.THR, true),
    TRP("???", AminoAcid.TRP, true),
    TYR("???", AminoAcid.TYR, true),
    VAL("???", AminoAcid.VAL, true),
    STOP("???", AminoAcid.STOP, true),

    PAD("???", AminoAcid.PAD, true);

    private String bases;
    private AminoAcid aa;
    private List<Base> baseSeq;
    private static Map<AminoAcid, List<Codon>> acidCodonMap = new HashMap<>();

    static {
        for (Codon c : Codon.values()) {
            if (c.bases.equals("???")) continue;
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
        this(bases, aa, false);
    }

    Codon(String bases, AminoAcid aa, boolean isWildcard) {
        this.bases = bases;
        this.aa = aa;
        if (isWildcard) {
            aa.setWildcard(this);
        }
        List<Base> baseList = new ArrayList<>();
        if (isWildcard) {
            if (this.aa == AminoAcid.PAD) {
                baseList.addAll(Arrays.asList(Base.Z, Base.Z, Base.Z));
            } else {
                baseList.addAll(Arrays.asList(Base.N, Base.N, Base.N));
            }
        } else {
            for (char c : bases.toCharArray()) {
                baseList.add(Base.valueOf(String.valueOf(c)));
            }
        }
        this.baseSeq = Collections.unmodifiableList(baseList);

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
        return bases.equals("???") ? aa.name() : getBases();
    }

    public List<Base> toBases() {
        return baseSeq;
    }
}
