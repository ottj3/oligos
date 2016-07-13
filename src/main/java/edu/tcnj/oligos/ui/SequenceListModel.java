package edu.tcnj.oligos.ui;

import edu.tcnj.oligos.library.Sequence;

import javax.swing.*;

abstract class SequenceListModel extends DefaultListModel<String> {
    abstract Sequence getActualAt(int i);
}
