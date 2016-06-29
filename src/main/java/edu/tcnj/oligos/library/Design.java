package edu.tcnj.oligos.library;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Stores a map of fragment ranges to delta levels for one amino acid
 */
public class Design {

    private Map<Fragment.Range, List<Integer>> deltasForRange = new HashMap<>();

    public Design(Map<Fragment.Range, List<Integer>> deltasForRange) {
        this.deltasForRange = deltasForRange;
    }

    List<Integer> getDeltasForPosition(int pos) {
        for (Map.Entry<Fragment.Range, List<Integer>> range : deltasForRange.entrySet()) {
            if (range.getKey().contains(pos)) return range.getValue();
        }
        return null;
    }

    List<Integer> getDeltasForRange(Fragment.Range range) {
        for (Map.Entry<Fragment.Range, List<Integer>> entry : deltasForRange.entrySet()) {
            if (entry.getKey().contains(range)) {
                return entry.getValue();
            }
        }
        return null;
    }

    Iterator<Map.Entry<Fragment.Range, List<Integer>>> iterator() {
        return deltasForRange.entrySet().iterator();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("{");
        Iterator<Map.Entry<Fragment.Range, List<Integer>>> it = deltasForRange.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry<Fragment.Range, List<Integer>> entry = it.next();
            Fragment.Range range = entry.getKey();
            List<Integer> deltas = entry.getValue();
            sb.append("(");
            sb.append(range.getStartPosition());
            sb.append(",");
            sb.append(range.getEndPosition());
            sb.append(")");
            sb.append(" -> ");
            sb.append(Arrays.toString(deltas.toArray()));
            if (it.hasNext()) sb.append(",\t");
        }
        sb.append("}");
        return sb.toString();
    }
}
