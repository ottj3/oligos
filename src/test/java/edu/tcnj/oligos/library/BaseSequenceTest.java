package edu.tcnj.oligos.library;

import com.google.common.collect.Lists;
import edu.tcnj.oligos.data.Base;
import org.junit.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class BaseSequenceTest {
    @Test
    public void testContains() {
        BaseSequence one = new BaseSequence(Lists.asList(Base.A, Base.C, new Base[]{Base.N, Base.N, Base.N}));
        BaseSequence two = new BaseSequence(Lists.asList(Base.A, Base.C, new Base[]{Base.T, Base.T, Base.T}));
        BaseSequence three = new BaseSequence(Lists.asList(Base.T, Base.T, new Base[]{Base.T, Base.A}));
        BaseSequence four = new BaseSequence(Lists.asList(Base.T, Base.T, new Base[]{Base.A}));

        assertTrue(one.contains(two));
        assertTrue(two.contains(one));
        assertTrue(three.contains(four));
        assertFalse(four.contains(three));
    }
}
