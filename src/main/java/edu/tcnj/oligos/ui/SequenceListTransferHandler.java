package edu.tcnj.oligos.ui;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import edu.tcnj.oligos.library.Oligo;
import edu.tcnj.oligos.library.Sequence;

import javax.swing.*;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.util.List;

class SequenceListTransferHandler extends TransferHandler {
    protected Transferable createTransferable(JComponent c) {
        JList list = (JList) c;
        int[] indices = list.getSelectedIndices();
        List<String> seqs = Lists.newArrayListWithCapacity(indices.length);
        for (int index : indices) {
            Sequence listElem = ((SequenceListModel) list.getModel()).getActualAt(index);
            seqs.add(listElem.toString());
        }
        String value = Joiner.on("\n").join(seqs);
        return new StringSelection(value);
    }

    @Override
    public int getSourceActions(JComponent c) {
        return TransferHandler.COPY;
    }
}
