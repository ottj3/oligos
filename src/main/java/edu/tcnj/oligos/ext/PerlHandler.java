package edu.tcnj.oligos.ext;

import java.io.*;

public final class PerlHandler {

    private PerlHandler() {}

    private static boolean hasPerl;
    static {
        ProcessBuilder pb = new ProcessBuilder("perl", "-v");
        try {
            pb.start();
            hasPerl = true;
        } catch (IOException ignored) {
            System.out.println("Perl not found. Codon pair scores will be unknown.");
            hasPerl = false;
        }
    }

    public static Double getCodonPairScore(String gene) {
        if (!hasPerl) return null;
        File cpsScript = new File("cps.pl");
        if (!cpsScript.exists() || !cpsScript.canExecute()) {
            return null;
        }
        File temp = new File("tmp.fa");
        StringBuilder content = new StringBuilder();
        content.append(">\n").append(gene);

        try {
            FileWriter fw = new FileWriter(temp);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write(content.toString());
            bw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
            temp.delete();
            return null;
        }
        ProcessBuilder pb = new ProcessBuilder("perl","cps.pl", "-p e.coli.k12.codon_pair_info.expected", "tmp.fa");
        StringBuilder plOutSB = null;
        try {
            Process proc = pb.start();
            int ev = proc.waitFor();
            if (ev == 0) {
                plOutSB = new StringBuilder();
                BufferedReader br = new BufferedReader(new InputStreamReader(proc.getInputStream()));
                String line;
                while ((line = br.readLine()) != null) {
                    plOutSB.append(line);
                }
                br.close();
            }
            temp.delete();
        } catch (IOException e) {
            e.printStackTrace();
            temp.delete();
            return null;
        } catch (InterruptedException e) {
            e.printStackTrace();
            temp.delete();
            return null;
        }

        if (plOutSB == null) {
            temp.delete();
            return null;
        }
        String plOut = plOutSB.toString();
        String[] outSplit = plOut.split("\\|");
        return outSplit.length == 3 ? Double.valueOf(outSplit[2]) : null;
    }
}
