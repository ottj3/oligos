package edu.tcnj.oligos.ui;

import com.google.common.base.Joiner;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Codon;
import edu.tcnj.oligos.ext.PythonHandler;
import edu.tcnj.oligos.library.Design;
import edu.tcnj.oligos.library.Library;
import edu.tcnj.oligos.library.Oligo;
import edu.tcnj.oligos.library.Protein;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

class Runner {
    private String seq;
    private int start;
    private int end;
    private int offset;
    private int oligoSize;
    private int overlapSize;
    private List<String> codons = new ArrayList<>();
    private List<Double> mins = new ArrayList<>();
    private List<Double> maxs = new ArrayList<>();
    private List<Integer> numLevels = new ArrayList<>();
    private String res;
    private Library lastLib;

    Runner(String seq, int start, int end, int offset, int oligoSize, int overlapSize,
                  List<String> codons, List<Double> mins, List<Double> maxs, List<Integer> numLevels) {
        this.seq = seq;
        this.start = start;
        this.end = end == 0 ? seq.length() : end;
        this.offset = offset;
        this.oligoSize = oligoSize;
        this.overlapSize = overlapSize;
        this.codons = codons;
        this.mins = mins;
        this.maxs = maxs;
        this.numLevels = numLevels;
    }

    void run() {
        String trimmedRNA = seq.substring(start, end);
        Protein gene = new Protein(trimmedRNA);
        String proteinString = Joiner.on("").join(gene.getAminoAcidSequence());
        for (int i = offset/3; i < 0; i++) {
            proteinString = "*" + proteinString;
        }
        String aoi = "";
        for (String s : codons) {
            aoi += Codon.valueOf(s).getAminoAcid().getCh();
        }

        PythonHandler pyth = new PythonHandler(proteinString, aoi, oligoSize / 3, overlapSize / 3,
                Doubles.toArray(mins), Doubles.toArray(maxs), Ints.toArray(numLevels));
        Map<AminoAcid, Design> designMap = pyth.run();


        // get amino acids from codons
        Map<AminoAcid, Codon> codonsOfInterest = Maps.newEnumMap(AminoAcid.class);
        for (String codonStr : codons) {
            Codon codon = Codon.valueOf(codonStr);
            AminoAcid acid = codon.getAminoAcid();
            codonsOfInterest.put(acid, codon);
        }
        Map<Codon, Design> codonDesignMap = Maps.newEnumMap(Codon.class);
        for (Map.Entry<AminoAcid, Design> entry : designMap.entrySet()) {
            codonDesignMap.put(codonsOfInterest.get(entry.getKey()), entry.getValue());
        }

        Map<Codon, Double> baseFrequencyMap = Maps.newEnumMap(Codon.class);
        for (int i = 0; i < codons.size(); i++) {
            Codon codon = Codon.valueOf(codons.get(i));
            double baseFreq = mins.get(i);
            baseFrequencyMap.put(codon, baseFreq);
        }

        Library.Builder builder = new Library.Builder();
        // calculate the number of extra characters we need to pad in order to get a full oligo at the last position
        int smalligo = oligoSize - overlapSize;
        int offsetToEndOfLastSmalligo = (((end - (start + offset)) - overlapSize));
        int rnaEnd = (oligoSize - (offsetToEndOfLastSmalligo % smalligo)) + offsetToEndOfLastSmalligo + offset;
        Library lib = builder
                .withProteinFromRNA(trimmedRNA)
                .withOligoSize(oligoSize / 3, overlapSize / 3)
                .withSequenceLength(offset / 3, rnaEnd / 3)
                .withCodonsOfInterest(codonsOfInterest)
                .withDesigns(codonDesignMap)
                .withMinFrequencies(baseFrequencyMap)
                .build();
        this.lastLib = lib;

        // run functions
        lib.createOligos();
        lib.fillFragments();

        lib.createOverlaps();
        lib.makeOverlapsUnique();

        // print output oligos
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<Integer, List<Oligo>> entry : lib.getOligos().entrySet()) {
            for (Oligo oligo : entry.getValue()) {
                StringBuilder oligoStr = new StringBuilder();
                for (Codon c : oligo.getSequence()) {
                    oligoStr.append(c.getBases());
                }
                sb.append(entry.getKey())
                        .append(" ")
                        .append(oligoStr.toString())
                        .append("\n");
            }
        }
        res = sb.toString();
    }

    public String getOligoOutput() {
        return res;
    }

    public Library getLastLib() {
        return lastLib;
    }
}
