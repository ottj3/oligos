package edu.tcnj.oligos.ext;

import jep.Jep;
import jep.JepException;

import java.util.Collection;

public class PythonHandler {

    public Object getDesigns(String protein, String acids, int segmentLength, int overlapLength, double[] mins, double[] maxs, int[] numLevels) {
        return runScript(protein, acids, segmentLength, overlapLength, mins, maxs, numLevels);
    }

    public Object getExampleDesigns() {
        String P = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGVQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK";
        String AOI = "TEV";
        int SEGSIZE = 30;
        int OVERLAPSIZE = 6;
        double[] MINS = {.1, .1, .1};
        double[] MAXS = {.85, .85, .85};
        int[] NUMPTS = {4, 4, 6};
        return runScript(P, AOI, SEGSIZE, OVERLAPSIZE, MINS, MAXS, NUMPTS);
    }

    private Object runScript(String protein, String acids, int segmentLength, int overlapLength, double[] mins, double[] maxs, int[] numLevels) {
        try {
            Jep jep = new Jep(false);
            jep.runScript("design.py");
            Object res = jep.invoke("compute_best_design", protein, acids, segmentLength, overlapLength, mins, maxs, numLevels);
            return res;
        } catch (JepException e) {
            e.printStackTrace();
        }
        return null;
    }

    // test
    public static void main(String[] args) {
        Object res = (new PythonHandler()).getExampleDesigns();
        if (res instanceof Collection) {
            for (Object o : ((Collection) res)) {
                System.out.println(o);
            }
        }
    }
}