package edu.tcnj.oligos;

import com.google.common.base.Joiner;
import com.google.common.base.Strings;
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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class OligoDesigner {
    public static void main(String[] args) throws Exception {
        File file = new File("input.txt");
        FileInputStream fis = new FileInputStream(file);
        BufferedReader br = new BufferedReader(new InputStreamReader(fis));

        String seq;
        int start;
        int end;
        int offset;
        int oligoSize;
        int overlapSize;
        List<String> codons = new ArrayList<>();
        List<Double> mins = new ArrayList<>();
        List<Double> maxs = new ArrayList<>();
        List<Integer> numLevels = new ArrayList<>();

        String line;
        String[] split;

        // read sequence
        seq = br.readLine();

        // start and end
        line = br.readLine();
        if (!Strings.isNullOrEmpty(line)) {
            split = line.split(" ");
            start = Integer.parseInt(split[0]);
            end = Integer.parseInt(split[1]);
        } else {
            start = 0;
            end = seq.length() - 1;
        }

        // offset
        line = br.readLine();
        if (!Strings.isNullOrEmpty(line)) {
            offset = Integer.parseInt(line);
        } else {
            offset = 0;
        }

        // size and overlap
        line = br.readLine();
        if (!Strings.isNullOrEmpty(line)) {
            split = line.split(" ");
            oligoSize = Integer.parseInt(split[0]);
            overlapSize = Integer.parseInt(split[1]);
        } else {
            throw new IllegalArgumentException("Input had no oligo size or overlap length at expected line.");
        }

        // codons
        while ((line = br.readLine()) != null) {
            split = line.split(" ");
            codons.add(split[0]);
            mins.add(Double.parseDouble(split[1]));
            maxs.add(Double.parseDouble(split[2]));
            numLevels.add(Integer.parseInt(split[3]));
        }
        br.close();

        String trimmedRNA = seq.substring(start, end - 2);
        Protein gene = new Protein(trimmedRNA);
        String proteinString = Joiner.on("").join(gene.getAminoAcidSequence());
        String aoi = "";
        for (String s : codons) {
            aoi += Codon.valueOf(s).getAminoAcid().getCh();
        }

        PythonHandler pyth = new PythonHandler(proteinString, aoi, oligoSize / 3, overlapSize / 3,
                Doubles.toArray(mins), Doubles.toArray(maxs), Ints.toArray(numLevels));
        Map<AminoAcid, Design> designMap = pyth.run();


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

        lib.createOligos();
        lib.fillFragments();

        lib.createOverlaps();
        lib.makeOverlapsUnique();
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
        System.out.println(sb.toString());
    }
    // rna
    // start end
    // offset
    // size overlap
    // codon min max num
    // codon ...
    //
}
