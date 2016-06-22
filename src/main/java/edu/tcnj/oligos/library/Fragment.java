package edu.tcnj.oligos.library;

import com.google.common.base.MoreObjects;

public class Fragment {

    private int startPosition;
    private int endPosition;

    public int getStartPosition() {
        return startPosition;
    }

    public int getEndPosition() {
        return endPosition;
    }

    @Override
    public String toString() {
        return MoreObjects.toStringHelper(this).add("start", startPosition).add("end", endPosition).toString();
    }

}
