package edu.tcnj.oligos.library;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;

public class SwapIterator {
    private Map<Integer, List<Integer>> potentialSwaps;
    private Map<Integer, Integer> swapIndices;
    private int baseValue;

    /**
     * Construct a SwapIterator for the given map of potentialSwaps,
     * using a base value of -1 (indicating a do-nothing state)
     *
     * @param potentialSwaps a map from an int (position) to a list of ints
     *                       (positions of codons with the same acid)
     */
    public SwapIterator(Map<Integer, List<Integer>> potentialSwaps) {
        this(potentialSwaps, -1);
    }

    /**
     * Construct a SwapIterator for the given map of PotentialSwaps and given baseValue
     *
     * @param potentialSwaps a map from an int (position) to a list of ints
     *                       (positions of codons with the same acid)
     * @param baseValue      the starting value for each iterator. Values &lt; 0 mean do nothing
     *                       (more negative means more do nothing states),
     *                       other values mean skip indices less than baseValue
     */
    public SwapIterator(Map<Integer, List<Integer>> potentialSwaps, int baseValue) {
        this.potentialSwaps = potentialSwaps;
        this.baseValue = baseValue;
        initializeIndices();
    }

    private void initializeIndices() {
        //Set the indices to baseValue for every position
        swapIndices = new HashMap<>();
        for (Map.Entry<Integer, List<Integer>> entry : potentialSwaps.entrySet()) {
            swapIndices.put(entry.getKey(), baseValue);
        }
    }

    /**
     * Get the current value without advancing the list
     *
     * @return A map from int (position) to int (a possible swap position)
     */
    public Map<Integer, Integer> peek() {
        if (!this.hasNext()) throw new NoSuchElementException("Out of potential swaps.");
        Map<Integer, Integer> currentSwap = new HashMap<>();
        //For every position
        for (Map.Entry<Integer, Integer> entry : swapIndices.entrySet()) {
            //If the value is negative (a do-nothing state) put that value in the map
            if (entry.getValue() <= -1) {
                currentSwap.put(entry.getKey(), entry.getValue());
            } else {
                //Otherwise, get the potentialSwap that this entry's indices refer to, and put it in the map
                currentSwap.put(entry.getKey(), potentialSwaps.get(entry.getKey()).get(entry.getValue()));
            }
        }
        return currentSwap;
    }

    /**
     * Get the current value and advance to the next permutation
     *
     * @return A map from int (position) to int (a possible swap position)
     */
    public Map<Integer, Integer> next() {
        //Get the current element to return it
        Map<Integer, Integer> nextVal = this.peek();
        //Step to the next permutation
        Iterator<Map.Entry<Integer, Integer>> it = swapIndices.entrySet().iterator();
        Map.Entry<Integer, Integer> next = it.next();
        int key = next.getKey();
        int value = next.getValue();
        swapIndices.put(key, value + 1);
        //If the current index overflows, reset it to baseValue and try incrementing the next index instead
        //(Essentially counting with the least significant bits on the left)
        //If the last index overflows, there are no permutations left so hasNext will later return false
        while (swapIndices.get(key) >= potentialSwaps.get(key).size()
                && it.hasNext()) {
            swapIndices.put(key, baseValue);
            next = it.next();
            key = next.getKey();
            value = next.getValue();
            swapIndices.put(key, value + 1);
        }
        return nextVal;
    }

    /**
     * Determine whether the iterator has run out of permutations
     *
     * @return true if there are still possible permutations, false otherwise
     */
    public boolean hasNext() {
        for (Map.Entry<Integer, Integer> entry : swapIndices.entrySet()) {
            if (entry.getValue() >= potentialSwaps.get(entry.getKey()).size()) return false;
        }
        return true;
    }
}
