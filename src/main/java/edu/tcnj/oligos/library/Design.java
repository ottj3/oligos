package edu.tcnj.oligos.library;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

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
}
