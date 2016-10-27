package edu.tcnj.oligos.ui;

import com.google.common.base.Strings;
import edu.tcnj.oligos.data.Codon;
import edu.tcnj.oligos.library.BaseSequence;
import edu.tcnj.oligos.library.Library;
import edu.tcnj.oligos.library.Oligo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class FileInputOligoDesigner {
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
        int minOverlapDiffs;
        List<String> codons = new ArrayList<>();
        List<Double> mins = new ArrayList<>();
        List<Double> maxs = new ArrayList<>();
        List<Integer> numLevels = new ArrayList<>();
        List<BaseSequence> restrictions;

        String line;
        String[] split;

        // read sequence
        seq = br.readLine();

        // start end offset
        line = br.readLine();
        if (!Strings.isNullOrEmpty(line)) {
            split = line.split(" ");
            start = Integer.parseInt(split[0]);
            end = Integer.parseInt(split[1]);
            offset = Integer.parseInt(split[2]);
        } else {
            start = 0;
            end = 0;
            offset = 0;
        }

        // size overlap minOverlapDiffs
        line = br.readLine();
        if (!Strings.isNullOrEmpty(line)) {
            split = line.split(" ");
            oligoSize = Integer.parseInt(split[0]);
            overlapSize = Integer.parseInt(split[1]);
            minOverlapDiffs = Integer.parseInt(split[2]);
        } else {
            throw new IllegalArgumentException("Input had no oligo size or overlap length at expected line.");
        }

        // restriction sites
        line = br.readLine();
        if (!Strings.isNullOrEmpty(line)) {
            restrictions = Util.getRestrictionSites(line.split(","));
        } else {
            restrictions = null;
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

        Runner runner = new Runner("design.py", seq, start, end, offset, oligoSize, overlapSize, codons,
                mins, maxs, numLevels, restrictions, minOverlapDiffs);

        runner.run();
        Library lib = runner.getLastLib();

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
        System.out.println(sb.toString());
    }

}
