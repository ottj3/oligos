package edu.tcnj.oligos.ui;

import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;
import edu.tcnj.oligos.data.Codon;
import edu.tcnj.oligos.library.Design;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableColumn;
import java.awt.*;
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
        if ("true".equals(System.getProperty("oligoDesigner.fillInTestData"))) {
            setupTestData();
        }
    }

    private void setupTestData() {
        rnaInputField.setText("ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGGTTATGGTGTTCAATGCTTTGCGAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCTGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAA");
        oligoLengthSpinner.setValue(90);
        overlapSizeSpinner.setValue(18);
        ((DefaultTableModel) codonTable.getModel()).addRow(new Object[]{"ACT", 0.1, 0.85, 4});
        ((DefaultTableModel) codonTable.getModel()).addRow(new Object[]{"GAA", 0.1, 0.85, 6});
        ((DefaultTableModel) codonTable.getModel()).addRow(new Object[]{"GTT", 0.1, 0.85, 4});
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
        Font font = Font.decode("Courier 16");
        outputArea.setFont(font);

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
        final DefaultTableModel codonTableModel = ((DefaultTableModel) codonTable.getModel());
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

//        final InputVerifier verifier = new InputVerifier() {
//            @Override
//            public boolean verify(JComponent jComponent) {
//                JTextField field = ((JTextField) jComponent);
//                String text = field.getText();
//            }
//        };
//        codonTable.getColumn("Min Freq").setCellEditor(new DefaultCellEditor(new JTextField()) {
//            @Override
//            public boolean stopCellEditing() {
//                return verifier.verify(editorComponent) && super.stopCellEditing();
//            }
//        });
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

        calculateOligosButton.addActionListener(new ActionListener() {
            private int val(Object obj) {
                return (Integer) obj;
            }

            private List<String> codons() {
                List<String> codons = Lists.newArrayList();
                for (int i = 0; i < codonTable.getRowCount(); i++) {
                    codons.add(String.valueOf(codonTable.getValueAt(i, 0)));
                }
                return codons;
            }

            private List<Double> freqs(int col) {
                List<Double> freqs = Lists.newArrayList();
                for (int i = 0; i < codonTable.getRowCount(); i++) {
                    freqs.add(Double.parseDouble(String.valueOf(codonTable.getValueAt(i, col))));
                }
                return freqs;
            }

            private List<Integer> levels() {
                List<Integer> levels = Lists.newArrayList();
                for (int i = 0; i < codonTable.getRowCount(); i++) {
                    levels.add(Integer.parseInt(String.valueOf(codonTable.getValueAt(i, 3))));
                }
                return levels;
            }

            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                try {
                    Runner runner = new Runner(rnaInputField.getText(),
                            val(seqStartSpinner.getValue()),
                            val(seqEndSpinner.getValue()),
                            val(seqOffsetSpinner.getValue()),
                            val(oligoLengthSpinner.getValue()),
                            val(overlapSizeSpinner.getValue()),
                            codons(), freqs(1), freqs(2), levels()
                    );
                    runner.run();
                    outputArea.setText("Design\n");
                    outputArea.append("Codon\tRange->Deltas...\n");
                    for (Map.Entry<Codon, Design> entry : runner.getLastLib().getDesigns().entrySet()) {
                        outputArea.append(entry.getKey() + "\t" + entry.getValue() + "\n");
                    }
                    outputArea.append("\nOligos\n");
                    outputArea.append(runner.getOligoOutput());
                } catch (Exception ex) {
                    outputArea.setText("There was an error running the program.\n"
                            + ex.getMessage()
                            + "\nPlease see console for details.");
                    ex.printStackTrace();
                }
            }
        });
    }

}
