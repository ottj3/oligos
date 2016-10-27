package edu.tcnj.oligos.ui;

import com.google.common.collect.Lists;
import edu.tcnj.oligos.data.Base;
import edu.tcnj.oligos.library.BaseSequence;

import java.util.List;

public class Util {
    static List<BaseSequence> getRestrictionSites(String[] sites) {
        List<BaseSequence> restrictions = Lists.newArrayList();
        for (String site : sites) {
            site = site.trim().toUpperCase();
            if (site.isEmpty()) continue;
            List<Base> restriction = Lists.newArrayList();
            for (char c : site.toCharArray()) {
                try {
                    Base b = Base.valueOf(String.valueOf(c));
                    restriction.add(b);
                } catch (IllegalArgumentException e) {
                    throw new IllegalArgumentException(
                            "Unrecognized or unsupported base " + c + " in restriction sites.", e);
                }
            }
            restrictions.add(new BaseSequence(restriction));
        }
        return restrictions;
    }
}