package edu.tcnj.oligos.ui;

import com.google.common.primitives.Ints;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableColumn;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.*;
import java.util.List;
import java.util.regex.Pattern;

public class OligoDesigner {
    private ResourceBundle res = ResourceBundle.getBundle("edu.tcnj.oligos.ui.oligoDesigner");
    private static final Pattern rnaPattern = Pattern.compile("^[ACTG]*");

    private JTextField rnaInputField;
    private JSpinner oligoLengthSpinner;
    private JSpinner overlapSizeSpinner;
    private JTable codonTable;
    private JTextArea outputArea;
    private JButton resetAllInputButton;
    private JButton calculateOligosButton;
    private JLabel rnaSequenceLabel;
    private JLabel rnaSequenceLengthLabel;
    private JLabel oligoLengthLabel;
    private JLabel overlapSizeLabel;
    private JPanel mainPanel;
    private JPanel rnaSeqPanel;
    private JPanel oligoSizePanel;
    private JPanel oligoLengthPanel;
    private JPanel overlapSizePanel;
    private JPanel offsetsPanel;
    private JPanel codonsPanel;
    private JSpinner seqStartSpinner;
    private JSpinner seqEndSpinner;
    private JSpinner seqOffsetSpinner;
    private JButton remCodonButton;
    private JButton addCodonButton;

    public static void main(String[] args) {
        JFrame frame = new JFrame("OligoDesigner");
        frame.setContentPane(new OligoDesigner().mainPanel);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setVisible(true);
    }

    public OligoDesigner() {
        createUIComponents();
    }

    private void createUIComponents() {
        rnaSequenceLengthLabel.setText(String.format(res.getString("label.sequenceLength"), 0, 0));
        rnaInputField.setInputVerifier(new InputVerifier() {
            @Override
            public boolean verify(JComponent input) {
                JTextField field = ((JTextField) input);
                String text = field.getText();
                if (text.length() % 3 == 0 && rnaPattern.matcher(text).matches()) {
                    rnaSequenceLengthLabel.setText(String.format(
                            res.getString("label.sequenceLength"), text.length() / 3, text.length()));
                    return true;
                } else {
                    rnaSequenceLengthLabel.setText("Invalid input RNA!");
                    return false;
                }
            }
        });
        rnaInputField.addKeyListener(new KeyAdapter() {
            @Override
            public void keyPressed(KeyEvent e) {
                super.keyPressed(e);
                if (e.getKeyCode() == KeyEvent.VK_ENTER || e.getKeyCode() == KeyEvent.VK_TAB) {
                    rnaInputField.getInputVerifier().verify(rnaInputField);
                }
            }
        });

        oligoLengthSpinner.setModel(new SpinnerNumberModel());
        ((JSpinner.DefaultEditor) oligoLengthSpinner.getEditor()).getTextField().setInputVerifier(new InputVerifier() {
            @Override
            public boolean verify(JComponent jComponent) {
                int val = ((SpinnerNumberModel) oligoLengthSpinner.getModel()).getNumber().intValue();
                if (val <= rnaInputField.getText().length() && val >= 0) {
                    return true;
                } else {
                    oligoLengthSpinner.setValue(val < 0 ? 0 : rnaInputField.getText().length());
                    return false;
                }
            }
        });
        overlapSizeSpinner.setModel(new SpinnerNumberModel());
        ((JSpinner.DefaultEditor) overlapSizeSpinner.getEditor()).getTextField().setInputVerifier(new InputVerifier() {
            @Override
            public boolean verify(JComponent jComponent) {
                int val = ((SpinnerNumberModel) overlapSizeSpinner.getModel()).getNumber().intValue();
                if (val <= ((Integer) oligoLengthSpinner.getValue()) && val >= 0) {
                    return true;
                } else {
                    overlapSizeSpinner.setValue(val < 0 ? 0 : oligoLengthSpinner.getValue());
                    return false;
                }
            }
        });

        // start, end, offset spinners

        // table
        DefaultTableModel codonTableModel = ((DefaultTableModel) codonTable.getModel());
        TableColumn col1 = new TableColumn(0);
        col1.setHeaderValue("Codon");
        codonTableModel.addColumn(col1);

        TableColumn col2 = new TableColumn(1);
        col2.setHeaderValue("Min Freq");
        codonTableModel.addColumn(col2);

        TableColumn col3 = new TableColumn(2);
        col3.setHeaderValue("Max Freq");
        codonTableModel.addColumn(col3);

        TableColumn col4 = new TableColumn(3);
        col4.setHeaderValue("# Levels");
        codonTableModel.addColumn(col4);

        codonTableModel.setColumnIdentifiers(new Object[]{"Codon", "Min Freq", "Max Freq", "# Levels"});

        addCodonButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                ((DefaultTableModel) codonTable.getModel()).addRow(new Object[]{"", 0.0, 0.0, 0});
            }
        });
        remCodonButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                List<Integer> removed = Ints.asList(codonTable.getSelectedRows());
                Collections.sort(removed);
                Collections.reverse(removed);
                for (int row : removed) {
                    ((DefaultTableModel) codonTable.getModel()).removeRow(row);
                }
            }
        });
    }

}
