package edu.tcnj.oligos.ui;

import com.google.common.base.Joiner;
import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;
import com.intellij.uiDesigner.core.GridConstraints;
import com.intellij.uiDesigner.core.GridLayoutManager;
import com.intellij.uiDesigner.core.Spacer;
import edu.tcnj.oligos.data.Base;
import edu.tcnj.oligos.data.Codon;
import edu.tcnj.oligos.library.BaseSequence;
import edu.tcnj.oligos.library.Design;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.ResourceBundle;
import java.util.regex.Pattern;

public class OligoDesigner {
    private ResourceBundle res = ResourceBundle.getBundle("edu.tcnj.oligos.ui.oligoDesigner");
    private static final Pattern rnaPattern = Pattern.compile("^[ACTG]*");

    private JTextField rnaInputField;
    private JSpinner oligoLengthSpinner;
    private JSpinner overlapSizeSpinner;
    private JTable codonTable;
    private JTextArea outputArea;
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
    private JPanel extraOptionsPanel;
    private JPanel restrictionSitesPanel;
    private JCheckBox restrictionSitesCheckBox;
    private JTextArea restrictionSitesText;

    public static void main(String[] args) {
        JFrame frame = new JFrame("OligoDesigner");
        frame.setContentPane(new OligoDesigner().mainPanel);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setSize(800, frame.getHeight());
        frame.setVisible(true);
    }

    private OligoDesigner() {
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

        InputMap im = codonTable.getInputMap(JTable.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
        KeyStroke tab = KeyStroke.getKeyStroke(KeyEvent.VK_TAB, 0);
        KeyStroke enter = KeyStroke.getKeyStroke(KeyEvent.VK_ENTER, 0);
        Action tabAction = new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                JTable t = (JTable) e.getSource();
                int column = t.getSelectedColumn();
                int row = t.getSelectedRow();
                do {
                    if (row == -1) row = 0;
                    if (column == -1) column = 0;
                    else column++;

                    if (column == t.getColumnCount()) {
                        column = 0;
                        row++;
                        if (row == t.getRowCount()) row = 0;
                    }
                } while (!t.isCellEditable(row, column));
                t.changeSelection(row, column, false, false);
                t.editCellAt(row, column);
                t.getEditorComponent().requestFocusInWindow();
                ((JTextField) t.getEditorComponent()).selectAll();
            }
        };
        codonTable.getActionMap().put(im.get(tab), tabAction);
        codonTable.getActionMap().put(im.get(enter), tabAction);

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

        restrictionSitesText.setEnabled(false);
        restrictionSitesCheckBox.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                AbstractButton abstractButton = (AbstractButton) actionEvent.getSource();
                restrictionSitesText.setEnabled(abstractButton.getModel().isSelected());
                if (restrictionSitesText.isEnabled()) {
                    restrictionSitesText.requestFocusInWindow();
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

            private List<BaseSequence> restrictionSites() {
                String text = restrictionSitesText.getText();
                if (!restrictionSitesCheckBox.getModel().isSelected() || Strings.isNullOrEmpty(text)) {
                    return null;
                }
                String[] sites = text.split(",");
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
                            throw new IllegalArgumentException("Unrecognized or unsupported base " + c
                                    + " in restriction sites.", e);
                        }
                    }
                    restrictions.add(new BaseSequence(restriction));
                }
                return restrictions;
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
                            codons(), freqs(1), freqs(2), levels(),
                            restrictionSites()
                    );
                    runner.run();
                    outputArea.setText("Design\n");
                    outputArea.append("Codon\tRange->Deltas...\n");
                    for (Map.Entry<Codon, Design> entry : runner.getLastLib().getDesigns().entrySet()) {
                        outputArea.append(entry.getKey() + "\t" + entry.getValue() + "\n");
                    }
                    List<BaseSequence> restrictions = runner.getLastLib().getRestrictions();
                    if (restrictions != null && !restrictions.isEmpty()) {
                        outputArea.append("\nRestriction Enzymes\n");
                        outputArea.append(Joiner.on(",").join(runner.getLastLib().getRestrictions()) + "\n");
                    }
                    outputArea.append("\nOligos\n");
                    outputArea.append(runner.getOligoOutput());
                } catch (Exception ex) {
                    outputArea.setText("There was an error running the program.\n"
                            + (ex.getMessage() == null ? "" : ex.getMessage() + "\n")
                            + "Please see console for details.");
                    ex.printStackTrace();
                }
            }
        });
    }

    {
// GUI initializer generated by IntelliJ IDEA GUI Designer
// >>> IMPORTANT!! <<<
// DO NOT EDIT OR ADD ANY CODE HERE!
        $$$setupUI$$$();
    }

    /**
     * Method generated by IntelliJ IDEA GUI Designer
     * >>> IMPORTANT!! <<<
     * DO NOT edit this method OR call it in your code!
     *
     * @noinspection ALL
     */
    private void $$$setupUI$$$() {
        mainPanel = new JPanel();
        mainPanel.setLayout(new GridLayoutManager(1, 2, new Insets(0, 0, 0, 0), -1, -1));
        final JPanel panel1 = new JPanel();
        panel1.setLayout(new GridLayoutManager(3, 1, new Insets(0, 0, 0, 0), -1, -1));
        mainPanel.add(panel1, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null, 0, false));
        final JPanel panel2 = new JPanel();
        panel2.setLayout(new GridLayoutManager(3, 1, new Insets(0, 0, 0, 0), -1, -1));
        panel1.add(panel2, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        rnaSeqPanel = new JPanel();
        rnaSeqPanel.setLayout(new GridLayoutManager(3, 1, new Insets(0, 0, 0, 0), -1, -1));
        panel2.add(rnaSeqPanel, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        rnaSequenceLabel = new JLabel();
        this.$$$loadLabelText$$$(rnaSequenceLabel, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("label.sequence"));
        rnaSeqPanel.add(rnaSequenceLabel, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        rnaInputField = new JTextField();
        rnaSeqPanel.add(rnaInputField, new GridConstraints(1, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_FIXED, null, new Dimension(150, -1), null, 0, false));
        rnaSequenceLengthLabel = new JLabel();
        this.$$$loadLabelText$$$(rnaSequenceLengthLabel, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("label.sequenceLength"));
        rnaSeqPanel.add(rnaSequenceLengthLabel, new GridConstraints(2, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        oligoSizePanel = new JPanel();
        oligoSizePanel.setLayout(new GridLayoutManager(1, 2, new Insets(0, 0, 0, 0), -1, -1));
        panel2.add(oligoSizePanel, new GridConstraints(1, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        oligoLengthPanel = new JPanel();
        oligoLengthPanel.setLayout(new GridLayoutManager(2, 1, new Insets(0, 0, 0, 0), -1, -1));
        oligoSizePanel.add(oligoLengthPanel, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        oligoLengthLabel = new JLabel();
        this.$$$loadLabelText$$$(oligoLengthLabel, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("label.oligoLength"));
        oligoLengthPanel.add(oligoLengthLabel, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        oligoLengthSpinner = new JSpinner();
        oligoLengthPanel.add(oligoLengthSpinner, new GridConstraints(1, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        overlapSizePanel = new JPanel();
        overlapSizePanel.setLayout(new GridLayoutManager(2, 2, new Insets(0, 0, 0, 0), -1, -1));
        oligoSizePanel.add(overlapSizePanel, new GridConstraints(0, 1, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        overlapSizeLabel = new JLabel();
        this.$$$loadLabelText$$$(overlapSizeLabel, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("label.overlapSize"));
        overlapSizePanel.add(overlapSizeLabel, new GridConstraints(0, 0, 1, 2, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        overlapSizeSpinner = new JSpinner();
        overlapSizePanel.add(overlapSizeSpinner, new GridConstraints(1, 0, 1, 2, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        offsetsPanel = new JPanel();
        offsetsPanel.setLayout(new GridLayoutManager(1, 3, new Insets(0, 0, 0, 0), -1, -1));
        panel2.add(offsetsPanel, new GridConstraints(2, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        final JPanel panel3 = new JPanel();
        panel3.setLayout(new GridLayoutManager(2, 1, new Insets(0, 0, 0, 0), -1, -1));
        offsetsPanel.add(panel3, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        final JLabel label1 = new JLabel();
        this.$$$loadLabelText$$$(label1, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("label.seqStart"));
        panel3.add(label1, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        seqStartSpinner = new JSpinner();
        panel3.add(seqStartSpinner, new GridConstraints(1, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        final JPanel panel4 = new JPanel();
        panel4.setLayout(new GridLayoutManager(2, 1, new Insets(0, 0, 0, 0), -1, -1));
        offsetsPanel.add(panel4, new GridConstraints(0, 1, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        final JLabel label2 = new JLabel();
        this.$$$loadLabelText$$$(label2, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("label.seqEnd"));
        panel4.add(label2, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        seqEndSpinner = new JSpinner();
        panel4.add(seqEndSpinner, new GridConstraints(1, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        final JPanel panel5 = new JPanel();
        panel5.setLayout(new GridLayoutManager(2, 1, new Insets(0, 0, 0, 0), -1, -1));
        offsetsPanel.add(panel5, new GridConstraints(0, 2, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        final JLabel label3 = new JLabel();
        this.$$$loadLabelText$$$(label3, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("label.seqOffset"));
        panel5.add(label3, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        seqOffsetSpinner = new JSpinner();
        panel5.add(seqOffsetSpinner, new GridConstraints(1, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        codonsPanel = new JPanel();
        codonsPanel.setLayout(new GridLayoutManager(2, 4, new Insets(0, 0, 0, 0), -1, -1));
        panel1.add(codonsPanel, new GridConstraints(1, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null, 0, false));
        final JLabel label4 = new JLabel();
        this.$$$loadLabelText$$$(label4, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("label.codonsOfInterest"));
        codonsPanel.add(label4, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        final JScrollPane scrollPane1 = new JScrollPane();
        codonsPanel.add(scrollPane1, new GridConstraints(1, 0, 1, 4, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null, 0, false));
        codonTable = new JTable();
        codonTable.setFillsViewportHeight(true);
        codonTable.putClientProperty("JTable.autoStartsEdit", Boolean.FALSE);
        scrollPane1.setViewportView(codonTable);
        remCodonButton = new JButton();
        this.$$$loadButtonText$$$(remCodonButton, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("button.remCodon"));
        codonsPanel.add(remCodonButton, new GridConstraints(0, 3, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        addCodonButton = new JButton();
        this.$$$loadButtonText$$$(addCodonButton, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("button.addCodon"));
        codonsPanel.add(addCodonButton, new GridConstraints(0, 2, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        final Spacer spacer1 = new Spacer();
        codonsPanel.add(spacer1, new GridConstraints(0, 1, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, 1, null, null, null, 0, false));
        extraOptionsPanel = new JPanel();
        extraOptionsPanel.setLayout(new GridLayoutManager(1, 2, new Insets(0, 0, 0, 0), -1, -1));
        panel1.add(extraOptionsPanel, new GridConstraints(2, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        restrictionSitesPanel = new JPanel();
        restrictionSitesPanel.setLayout(new GridLayoutManager(2, 1, new Insets(0, 0, 0, 0), -1, -1));
        extraOptionsPanel.add(restrictionSitesPanel, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        restrictionSitesCheckBox = new JCheckBox();
        this.$$$loadButtonText$$$(restrictionSitesCheckBox, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("button.enableRestrictions"));
        restrictionSitesPanel.add(restrictionSitesCheckBox, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        restrictionSitesText = new JTextArea();
        restrictionSitesText.setLineWrap(true);
        restrictionSitesText.setWrapStyleWord(true);
        restrictionSitesPanel.add(restrictionSitesText, new GridConstraints(1, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_WANT_GROW, null, new Dimension(150, 50), null, 0, false));
        final Spacer spacer2 = new Spacer();
        extraOptionsPanel.add(spacer2, new GridConstraints(0, 1, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, 1, null, null, null, 0, false));
        final JPanel panel6 = new JPanel();
        panel6.setLayout(new GridLayoutManager(2, 1, new Insets(0, 0, 0, 0), -1, -1));
        mainPanel.add(panel6, new GridConstraints(0, 1, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        final JScrollPane scrollPane2 = new JScrollPane();
        panel6.add(scrollPane2, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null, 0, false));
        outputArea = new JTextArea();
        outputArea.setEditable(false);
        outputArea.setFont(new Font("Courier", outputArea.getFont().getStyle(), 16));
        scrollPane2.setViewportView(outputArea);
        calculateOligosButton = new JButton();
        this.$$$loadButtonText$$$(calculateOligosButton, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("button.calculate"));
        panel6.add(calculateOligosButton, new GridConstraints(1, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
    }

    /**
     * @noinspection ALL
     */
    private void $$$loadLabelText$$$(JLabel component, String text) {
        StringBuffer result = new StringBuffer();
        boolean haveMnemonic = false;
        char mnemonic = '\0';
        int mnemonicIndex = -1;
        for (int i = 0; i < text.length(); i++) {
            if (text.charAt(i) == '&') {
                i++;
                if (i == text.length()) break;
                if (!haveMnemonic && text.charAt(i) != '&') {
                    haveMnemonic = true;
                    mnemonic = text.charAt(i);
                    mnemonicIndex = result.length();
                }
            }
            result.append(text.charAt(i));
        }
        component.setText(result.toString());
        if (haveMnemonic) {
            component.setDisplayedMnemonic(mnemonic);
            component.setDisplayedMnemonicIndex(mnemonicIndex);
        }
    }

    /**
     * @noinspection ALL
     */
    private void $$$loadButtonText$$$(AbstractButton component, String text) {
        StringBuffer result = new StringBuffer();
        boolean haveMnemonic = false;
        char mnemonic = '\0';
        int mnemonicIndex = -1;
        for (int i = 0; i < text.length(); i++) {
            if (text.charAt(i) == '&') {
                i++;
                if (i == text.length()) break;
                if (!haveMnemonic && text.charAt(i) != '&') {
                    haveMnemonic = true;
                    mnemonic = text.charAt(i);
                    mnemonicIndex = result.length();
                }
            }
            result.append(text.charAt(i));
        }
        component.setText(result.toString());
        if (haveMnemonic) {
            component.setMnemonic(mnemonic);
            component.setDisplayedMnemonicIndex(mnemonicIndex);
        }
    }

    /**
     * @noinspection ALL
     */
    public JComponent $$$getRootComponent$$$() {
        return mainPanel;
    }
}
