package edu.tcnj.oligos.ui;

import com.google.common.base.Joiner;
import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;
import com.intellij.uiDesigner.core.GridConstraints;
import com.intellij.uiDesigner.core.GridLayoutManager;
import com.intellij.uiDesigner.core.Spacer;
import edu.tcnj.oligos.data.AminoAcid;
import edu.tcnj.oligos.data.Base;
import edu.tcnj.oligos.data.Codon;
import edu.tcnj.oligos.library.*;
import sun.awt.AppContext;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.event.*;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.ResourceBundle;
import java.util.concurrent.Executors;
import java.util.regex.Pattern;

public class OligoDesigner {
    private ResourceBundle res = ResourceBundle.getBundle("edu.tcnj.oligos.ui.oligoDesigner");
    private static final Pattern rnaPattern = Pattern.compile("^[ACTG]*");
    private SwingWorker runnerThread = null;

    private JTextField rnaInputField;
    private JSpinner oligoLengthSpinner;
    private JSpinner overlapSizeSpinner;
    private JTable codonTable;
    private JTextArea outputInfoArea;
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
    private JTabbedPane outputTabs;
    private JList outputOligoList;
    private JTextArea outputOligoInfo;
    private JList outputGeneList;
    private JTextArea outputGeneInfo;
    private JButton cancelCalculateButton;

    public static void main(String[] args) {
        JFrame frame = new JFrame("Oligo Designer");
        frame.setContentPane(new OligoDesigner().mainPanel);
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        frame.pack();
        frame.setSize(1200, 600);
        frame.setVisible(true);
        AppContext.getAppContext().put(SwingWorker.class, Executors.newSingleThreadExecutor());
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
        rnaSequenceLengthLabel.setText(String.format(
                res.getString("label.sequenceLength"),
                rnaInputField.getText().length() / 3,
                rnaInputField.getText().length()));
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
                if (val <= rnaInputField.getText().length() && val > 0) {
                    return true;
                } else {
                    oligoLengthSpinner.setValue(val <= 0 ? 1 : rnaInputField.getText().length());
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
//                    // unfinished
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

        outputOligoList.setModel(new OligoListModel());
        outputOligoList.setTransferHandler(new SequenceListTransferHandler());

        outputGeneList.setModel(new GeneListModel());
        outputGeneList.setTransferHandler(new SequenceListTransferHandler());

        outputOligoList.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                super.mouseClicked(e);
                int index = outputOligoList.locationToIndex(e.getPoint());
                if (index == -1) return;
                Oligo oligo = ((OligoListModel) outputOligoList.getModel()).getActualAt(index);
                outputOligoInfo.setText("Sequence: ");
                outputOligoInfo.append(oligo.toString() + "\n\n");
                outputOligoInfo.append(String.format(res.getString("info.oligo.position"),
                        Integer.valueOf(((String) outputOligoList.getModel().getElementAt(index)).split(" ")[0]))
                        + "\n\n");
                String deltas = "Codon -> Delta\n";
                for (Map.Entry<Codon, Integer> entry : oligo.getDeltas().entrySet()) {
                    deltas += entry.getKey() + " (" + entry.getKey().getAminoAcid() + ") -> " + entry.getValue() + "\n";
                }
                outputOligoInfo.append(deltas);
            }
        });
        outputGeneList.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                super.mouseClicked(e);
                int index = outputGeneList.locationToIndex(e.getPoint());
                if (index == -1) return;
                Gene gene = ((GeneListModel) outputGeneList.getModel()).getActualAt(index);
                outputGeneInfo.setText("Sequence: ");
                outputGeneInfo.append(gene.toString() + "\n\n");
                String deltas = "Codon Frequencies\n\n";
                for (Map.Entry<AminoAcid, Map<Codon, Double>> acidEntry : gene.getFreqs().entrySet()) {
                    deltas += acidEntry.getKey().getName() + ":\n";
                    List<String> codonInfo = Lists.newArrayList();
                    for (Map.Entry<Codon, Double> codonEntry : acidEntry.getValue().entrySet()) {
                        codonInfo.add(codonEntry.getKey() + " -> " + ((int) (codonEntry.getValue() * 100 + 0.5)) + "%");
                    }
                    deltas += Joiner.on(", ").join(codonInfo) + "\n";
                }
                outputGeneInfo.append(deltas);
            }
        });
        DefaultListCellRenderer alternatingListRenderer = new DefaultListCellRenderer() {
            @Override
            public Component getListCellRendererComponent(JList jList, Object o, int i, boolean b, boolean b1) {
                super.getListCellRendererComponent(jList, o, i, b, b1);
                String[] split = String.valueOf(o).split(" ");
                if (Integer.valueOf(split[0]) % 2 == 0) {
                    this.setBackground(this.getBackground().darker());
                }
                return this;
            }
        };
        outputOligoList.setCellRenderer(alternatingListRenderer);
        //outputGeneList.setCellRenderer(alternatingListRenderer);

        cancelCalculateButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                if (runnerThread != null) {
                    runnerThread.cancel(true);
                }
                runnerThread = null;
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
                calculateOligosButton.setEnabled(false);
                calculateOligosButton.setVisible(false);
                final Runner runner = new Runner(rnaInputField.getText(),
                        val(seqStartSpinner.getValue()),
                        val(seqEndSpinner.getValue()),
                        val(seqOffsetSpinner.getValue()),
                        val(oligoLengthSpinner.getValue()),
                        val(overlapSizeSpinner.getValue()),
                        codons(), freqs(1), freqs(2), levels(),
                        restrictionSites()
                );
                final Runnable run = new Runnable() {
                    @Override
                    public void run() {
                        try {
                            runner.run();
                        } catch (final Exception ex) {
                            SwingUtilities.invokeLater(new Runnable() {
                                @Override
                                public void run() {
                                    handleException(ex);
                                }
                            });
                            return;
                        }
                        SwingUtilities.invokeLater(new Runnable() {
                            @Override
                            public void run() {
                                outputInfoArea.setText("Design\n");
                                outputInfoArea.append("Codon\tRange->Deltas...\n");
                                for (Map.Entry<Codon, Design> entry : runner.getLastLib().getDesigns().entrySet()) {
                                    outputInfoArea.append(entry.getKey() + "\t" + entry.getValue() + "\n");
                                }
                                List<BaseSequence> restrictions = runner.getLastLib().getRestrictions();
                                if (restrictions != null && !restrictions.isEmpty()) {
                                    outputInfoArea.append("\nRestriction Enzymes\n");
                                    outputInfoArea.append(Joiner.on(",").join(runner.getLastLib().getRestrictions()) + "\n");
                                }
                                ((OligoListModel) outputOligoList.getModel()).setOligos(runner.getLastLib().getOligos());
                                ((GeneListModel) outputGeneList.getModel()).addGenes(runner.getLastLib());
                                cancelCalculateButton.setVisible(false);
                                cancelCalculateButton.setEnabled(false);
                                calculateOligosButton.setVisible(true);
                                calculateOligosButton.setEnabled(true);
                            }
                        });
                    }
                };
                runnerThread = new SwingWorker<Void, Void>() {
                    @Override
                    protected Void doInBackground() {
                        run.run();
                        return null;
                    }
                };
                runnerThread.execute();

                cancelCalculateButton.setVisible(true);
                cancelCalculateButton.setEnabled(true);
            }
        });
    }

    private void handleException(Exception ex) {
        if (ex instanceof RuntimeException) {
            if (ex.getCause() instanceof InterruptedException) {
                outputInfoArea.setText(res.getString("exception.interrupt"));
            } else if (ex.getCause() instanceof OutOfSwapsException) {
                outputInfoArea.setText(res.getString("exception.noSwaps").replaceAll("\\\\n", "\n"));
            }
        } else {
            outputInfoArea.setText(res.getString("exception.generic")
                    + (ex.getMessage() == null ? "." : ":\n" + ex.getMessage()));
        }
        outputInfoArea.append("\nSee console for details.");
        ex.printStackTrace();

        cancelCalculateButton.setVisible(false);
        cancelCalculateButton.setEnabled(false);
        calculateOligosButton.setVisible(true);
        calculateOligosButton.setEnabled(true);
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
        mainPanel.setLayout(new GridLayoutManager(1, 1, new Insets(0, 0, 0, 0), -1, -1));
        final JSplitPane splitPane1 = new JSplitPane();
        mainPanel.add(splitPane1, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, new Dimension(200, 200), null, 0, false));
        final JPanel panel1 = new JPanel();
        panel1.setLayout(new GridLayoutManager(3, 1, new Insets(0, 0, 0, 0), -1, -1));
        splitPane1.setLeftComponent(panel1);
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
        rnaSeqPanel.add(rnaSequenceLengthLabel, new GridConstraints(2, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
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
        codonsPanel.add(scrollPane1, new GridConstraints(1, 0, 1, 4, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
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
        panel1.add(extraOptionsPanel, new GridConstraints(2, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null, 0, false));
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
        panel6.setLayout(new GridLayoutManager(2, 2, new Insets(0, 0, 0, 0), -1, -1));
        splitPane1.setRightComponent(panel6);
        final JSplitPane splitPane2 = new JSplitPane();
        splitPane2.setDividerLocation(200);
        splitPane2.setOrientation(0);
        panel6.add(splitPane2, new GridConstraints(0, 0, 1, 2, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, new Dimension(200, 200), null, 0, false));
        final JScrollPane scrollPane2 = new JScrollPane();
        splitPane2.setLeftComponent(scrollPane2);
        outputInfoArea = new JTextArea();
        outputInfoArea.setEditable(false);
        outputInfoArea.setFont(new Font("Courier", outputInfoArea.getFont().getStyle(), 16));
        scrollPane2.setViewportView(outputInfoArea);
        outputTabs = new JTabbedPane();
        splitPane2.setRightComponent(outputTabs);
        final JSplitPane splitPane3 = new JSplitPane();
        splitPane3.setDividerLocation(260);
        splitPane3.setOneTouchExpandable(true);
        outputTabs.addTab(ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("tab.oligoView"), splitPane3);
        final JPanel panel7 = new JPanel();
        panel7.setLayout(new GridLayoutManager(1, 1, new Insets(0, 0, 0, 0), -1, -1));
        splitPane3.setLeftComponent(panel7);
        final JScrollPane scrollPane3 = new JScrollPane();
        panel7.add(scrollPane3, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null, 0, false));
        scrollPane3.setBorder(BorderFactory.createTitledBorder(BorderFactory.createLineBorder(Color.black), null));
        outputOligoList = new JList();
        scrollPane3.setViewportView(outputOligoList);
        final JPanel panel8 = new JPanel();
        panel8.setLayout(new GridLayoutManager(1, 1, new Insets(0, 0, 0, 0), -1, -1));
        splitPane3.setRightComponent(panel8);
        final JScrollPane scrollPane4 = new JScrollPane();
        panel8.add(scrollPane4, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null, 0, false));
        outputOligoInfo = new JTextArea();
        scrollPane4.setViewportView(outputOligoInfo);
        final JSplitPane splitPane4 = new JSplitPane();
        splitPane4.setDividerLocation(260);
        splitPane4.setOneTouchExpandable(true);
        outputTabs.addTab(ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("tab.geneView"), splitPane4);
        final JPanel panel9 = new JPanel();
        panel9.setLayout(new GridLayoutManager(1, 1, new Insets(0, 0, 0, 0), -1, -1));
        splitPane4.setLeftComponent(panel9);
        final JScrollPane scrollPane5 = new JScrollPane();
        panel9.add(scrollPane5, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null, 0, false));
        scrollPane5.setBorder(BorderFactory.createTitledBorder(BorderFactory.createLineBorder(Color.black), null));
        outputGeneList = new JList();
        scrollPane5.setViewportView(outputGeneList);
        final JPanel panel10 = new JPanel();
        panel10.setLayout(new GridLayoutManager(1, 1, new Insets(0, 0, 0, 0), -1, -1));
        splitPane4.setRightComponent(panel10);
        final JScrollPane scrollPane6 = new JScrollPane();
        panel10.add(scrollPane6, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null, 0, false));
        outputGeneInfo = new JTextArea();
        scrollPane6.setViewportView(outputGeneInfo);
        calculateOligosButton = new JButton();
        this.$$$loadButtonText$$$(calculateOligosButton, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("button.calculate"));
        panel6.add(calculateOligosButton, new GridConstraints(1, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        cancelCalculateButton = new JButton();
        cancelCalculateButton.setEnabled(false);
        this.$$$loadButtonText$$$(cancelCalculateButton, ResourceBundle.getBundle("edu/tcnj/oligos/ui/oligoDesigner").getString("button.cancelCalc"));
        cancelCalculateButton.setVisible(false);
        panel6.add(cancelCalculateButton, new GridConstraints(1, 1, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
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
