package edu.tcnj.oligos.ui;

import com.google.common.base.Joiner;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Codon;
import edu.tcnj.oligos.ext.PythonHandler;
import edu.tcnj.oligos.library.BaseSequence;
import edu.tcnj.oligos.library.Design;
import edu.tcnj.oligos.library.Library;
import edu.tcnj.oligos.library.Protein;

import java.util.List;
import java.util.Map;

public class Runner {
    private String seq;
    private String scriptName;
    private int start;
    private int end;
    private int offset;
    private int oligoSize;
    private int overlapSize;
    private int minOverlapDiffs;
    private List<String> codons;
    private List<Double> mins;
    private List<Double> maxs;
    private List<Integer> numLevels;
    private List<BaseSequence> restrictions;
    private Library lastLib;

    public Runner(String scriptName, String seq, int start, int end, int offset, int oligoSize, int overlapSize,
           List<String> codons, List<Double> mins, List<Double> maxs, List<Integer> numLevels,
           List<BaseSequence> restrictions, int minOverlapDiffs) {
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
        this.restrictions = restrictions;
        this.scriptName = scriptName;
        this.minOverlapDiffs = minOverlapDiffs;
    }

    public void run() {
        this.lastLib = null;
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

        PythonHandler pyth = new PythonHandler(scriptName, proteinString, aoi, oligoSize / 3, overlapSize / 3,
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
                .withDifferencesNeeded(minOverlapDiffs)
                .withSequenceLength(offset / 3, rnaEnd / 3)
                .withCodonsOfInterest(codonsOfInterest)
                .withDesigns(codonDesignMap)
                .withRestrictions(restrictions)
                .build();
        this.lastLib = lib;

        // run functions
        lib.removeRestrictionEnzymes();
        lib.initBaseFrequencies(baseFrequencyMap);

        lib.createOligos();
        lib.fillFragments();

        lib.createOverlaps();
        lib.makeOverlapsUnique();
        lib.removePadding(start, offset, seq);
        lib.setExecutionPhase(Library.Phase.FINISHED);
    }

    public Library getLastLib() {
        return lastLib;
    }
}
