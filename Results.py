import sys


from PyQt5 import QtWidgets, uic, QtCore, QtGui, Qt
from bs4 import BeautifulSoup
import requests
import webbrowser
from Scoring import OnTargetScore
import Algorithms
from Algorithms import SeqTranslate
from CSPRparser import CSPRparser
import GlobalSettings
import os
from APIs import Kegg
import OffTarget
# =========================================================================================
# CLASS NAME: Results
# Inputs: Takes information from the main application window and displays the gRNA results
# Outputs: The display of the gRNA results search
# =========================================================================================


class Results(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        super(Results, self).__init__(parent)
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'resultsWindow.ui'), self)

        self.setWindowTitle('Results')
        self.geneViewer.setReadOnly(True)
        # Scoring Class object #
        self.onscore = OnTargetScore()
        self.S = SeqTranslate()
        self.dbpath = ""
        # Main Data container
        # Keys: Gene names
        # Values: #
        self.AllData = {}
        self.highlighted = {}

        self.co_target_endo_list = list()
        self.startpos = 0
        self.endpos = 0
        self.directory = ""
        self.geneDict = dict() # dictionary passed into transfer_data
        self.geneNTDict = dict() #dictionary passed into transfer_data, same key as geneDict, but hols the NTSEQ

        self.switcher = [1,1,1,1,1,1,1,1,1]  # for keeping track of where we are in the sorting clicking for each column

        # Target Table settings #
        self.targetTable.setColumnCount(9)  # hardcoded because there will always be nine columns
        self.targetTable.setShowGrid(False)
        self.targetTable.setHorizontalHeaderLabels("Location;Sequence;Strand;PAM;Score;Off-Target;Off-Target;Details;Endonuclease(s)".split(";"))
        self.targetTable.horizontalHeader().setSectionsClickable(True)
        self.targetTable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.targetTable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.targetTable.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

        self.back_button.clicked.connect(self.goBack)
        self.targetTable.horizontalHeader().sectionClicked.connect(self.table_sorting)
        self.actionSave.triggered.connect(self.save_data)
        self.actionOpen.triggered.connect(self.open_data)
        self.actionOff_Target_Analysis.triggered.connect(self.Off_Target_Analysis)
        self.actionCoTargeting.triggered.connect(self.open_cotarget)
        self.changeEndoButton.clicked.connect(self.changeEndonuclease)
        self.displayGeneViewer.stateChanged.connect(self.checkGeneViewer)
        self.cotarget_checkbox.stateChanged.connect(self.prep_cotarget_checkbox)
        self.highlight_gene_viewer_button.clicked.connect(self.highlight_gene_viewer)
        self.checkBoxSelectAll.stateChanged.connect(self.selectAll)
        self.pushButton_Deselect_All.clicked.connect(self.deselectAll)
        self.gene_viewer_settings_button.clicked.connect(self.changeGeneViewerSettings)
        self.change_start_end_button.clicked.connect(self.change_indicies)
        self.actionTo_CSV.triggered.connect(self.open_export_to_csv)


        #self.targetTable.itemSelectionChanged.connect(self.item_select)
        self.minScoreLine.setText("0")

        # Connecting the filters to the displayGeneData function
        self.fivegseqCheckBox.stateChanged.connect(self.displayGeneData)
        self.minScoreLine.textChanged.connect(self.displayGeneData)

        # Setting up the score filter:
        self.scoreSlider.setMinimum(0)
        self.scoreSlider.setMaximum(100)
        self.scoreSlider.setTracking(False)

        self.scoreSlider.valueChanged.connect(self.update_score_filter)

        #bool used to make sure only 1 instance of the OffTarget window is created
        self.first_boot = True
        #OTA is used to hold the row numbers of the items selected by user for OffTargetAnalysis
        #using this helps speed up updating the chart
        self.OTA = []



        self.detail_output_list = []
        self.rows_and_seq_list = []
        self.seq_and_avg_list = []
        self.files_list = []
        self.seq_finder_cspr_file = ''

        self.mwfg = self.frameGeometry() ##Center window
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center() ##Center window



    # this function opens the export_to_csv window
    # first it makes sure that the user actually has some highlighted targets that they want exported
    def open_export_to_csv(self):
        select_items = self.targetTable.selectedItems()
        if len(select_items) <= 0:
            QtWidgets.QMessageBox.question(self, "Nothing Selected",
                                           "No targets were highlighted."
                                           "Please highlight the targets you want to be exported to a CSV File!",
                                           QtWidgets.QMessageBox.Ok)
            return
        # now launch the window
        GlobalSettings.mainWindow.export_csv_window.launch(select_items)

    def change_indicies(self):

        # make sure the gene viewer is on
        if not self.displayGeneViewer.isChecked():
            QtWidgets.QMessageBox.question(self, "Gene Viewer Error",
                                           "Gene Viewer display is off! "
                                           "Please turn the Gene Viewer on in order to highlight the sequences selected",
                                           QtWidgets.QMessageBox.Ok)
            return


        # change the start and end values
        prevTuple = self.geneDict[self.comboBoxGene.currentText()]
        tempTuple = (self.geneDict[self.comboBoxGene.currentText()][0], int(self.lineEditStart.displayText()), int(self.lineEditEnd.displayText()))

        # make sure that the difference between indicies is not too large
        if abs(tempTuple[1] - tempTuple[2]) > 50000:
            #print(abs(tempTuple[1] - tempTuple[2]))
            QtWidgets.QMessageBox.question(self, "Sequence Too Long",
                                           "The sequence is too long! "
                                           "Please choose indicies that will make the sequence less than 50,000!",
                                           QtWidgets.QMessageBox.Ok)
            self.lineEditStart.setText(str(self.geneDict[self.comboBoxGene.currentText()][1]))
            self.lineEditEnd.setText(str(self.geneDict[self.comboBoxGene.currentText()][2]))
            return


        # if the user is using kegg
        if GlobalSettings.mainWindow.gene_viewer_settings.file_type == "kegg":
            tempTuple = (self.geneDict[self.comboBoxGene.currentText()][0], int(self.lineEditStart.displayText()), int(self.lineEditEnd.displayText()), self.geneDict[self.comboBoxGene.currentText()][3])
            self.geneDict[self.comboBoxGene.currentText()] = tempTuple
            # build the URL string
            # build the URL one way if the chrom from KEGG is -1
            if tempTuple[3] == -1:
                url = "https://www.genome.jp/dbget-bin/cut_sequence_genes.pl?FROM="
                url = url + str(self.geneDict[self.comboBoxGene.currentText()][1])
                url = url + "&TO="
                url = url + str(self.geneDict[self.comboBoxGene.currentText()][2])
                url = url + "&VECTOR=1&ORG="
                url = url + GlobalSettings.mainWindow.Annotations_Organism.currentText().split(" ")[1]
            # build it a different way if the chrom is not -1
            else:
                url = "https://www.genome.jp/dbget-bin/cut_sequence_genes.pl?FROM="
                url = url + str(self.geneDict[self.comboBoxGene.currentText()][1])
                url = url + "&TO="
                url = url + str(self.geneDict[self.comboBoxGene.currentText()][2])
                url = url + "&VECTOR=1&ORG="
                url = url + GlobalSettings.mainWindow.Annotations_Organism.currentText().split(" ")[1]
                url = url + '&CHR=' + tempTuple[3]


            # soup time
            source = requests.get(url)
            plain_text = source.text
            soup = BeautifulSoup(plain_text, "html.parser")


            buffer = ""

            # use the try in case user gives bad input
            try:
                for item in soup.find('pre'):
                    # get the first line
                    if ">" in item.string:
                        temp = item.string.replace("\n", "")
                        temp = temp.replace(" ", "")
                        newLineIndex = temp.find(")") + 1

                        if buffer == "":
                            buffer = temp[newLineIndex:]
                        else:
                            buffer = buffer + temp[newLineIndex:]
                    # get every other lie
                    else:
                        if buffer == "":
                            buffer = item.string
                        else:
                            buffer = buffer + item.string

                buffer = buffer.replace("\n", "")
                self.geneNTDict[self.comboBoxGene.currentText()] = buffer
            except:
                QtWidgets.QMessageBox.question(self, "Error",
                                               'An error has occured. It is possible that KEGG does not allow for NT Sequence changes for this organism.',
                                               QtWidgets.QMessageBox.Ok)
                self.lineEditStart.setText(str(self.geneDict[self.comboBoxGene.currentText()][1]))
                self.lineEditEnd.setText(str(self.geneDict[self.comboBoxGene.currentText()][2]))
                return
        # if the user is using gbff
        elif GlobalSettings.mainWindow.gene_viewer_settings.file_type == "gbff":
            self.geneDict[self.comboBoxGene.currentText()] = tempTuple
            sequence = GlobalSettings.mainWindow.gene_viewer_settings.gbff_sequence_finder(self.geneDict[self.comboBoxGene.currentText()])
            self.geneNTDict[self.comboBoxGene.currentText()] = sequence
        # if the user is using fna
        elif GlobalSettings.mainWindow.gene_viewer_settings.file_type == "fna":
            self.geneDict[self.comboBoxGene.currentText()] = tempTuple
            sequence = GlobalSettings.mainWindow.gene_viewer_settings.fna_sequence_finder(self.geneDict[self.comboBoxGene.currentText()])
            self.geneNTDict[self.comboBoxGene.currentText()] = sequence

        # check and see if we need to add lowercase letters
        changeInStart = tempTuple[1] - prevTuple[1]
        changeInEnd = tempTuple[2] - prevTuple[2]

        # check and see if the sequence is extended at all
        # if it is, make the extended part lower-case as opposed to upper case
        if changeInStart != 0 and changeInStart < 0:
            tempString = self.geneNTDict[self.comboBoxGene.currentText()][:abs(changeInStart)].lower()
            tempString = tempString + self.geneNTDict[self.comboBoxGene.currentText()][abs(changeInStart):]
            self.geneNTDict[self.comboBoxGene.currentText()] = tempString
        if changeInEnd != 0 and changeInEnd > 0:
            tempString = self.geneNTDict[self.comboBoxGene.currentText()][len(self.geneNTDict[self.comboBoxGene.currentText()]) - abs(changeInEnd):].lower()
            tempString = self.geneNTDict[self.comboBoxGene.currentText()][:len(self.geneNTDict[self.comboBoxGene.currentText()]) - abs(changeInEnd)] + tempString
            self.geneNTDict[self.comboBoxGene.currentText()] = tempString

        # update the gene viewer
        self.checkGeneViewer()

    # this function goes through and deselects all of the Off-Target checkboxes
    # also sets the selectAllShown check box to unchecked as well
    def deselectAll(self):
        for i in range(self.targetTable.rowCount()):
            self.targetTable.cellWidget(i, 6).setCheckState(0)
        self.checkBoxSelectAll.setCheckState(0)

    # this function listens for a stateChange in selectAllShown
    # if it is checked, it selects all shown
    # if it is unchecked, it deselects all shown
    # Note: it is a little buggy, possibly because when you change the minimum score it resets it all
    def selectAll(self):
        if self.checkBoxSelectAll.isChecked():
            for i in range(self.targetTable.rowCount()):
                self.targetTable.cellWidget(i, 6).setCheckState(2)
        elif not self.checkBoxSelectAll.isChecked():
            for i in range(self.targetTable.rowCount()):
                self.targetTable.cellWidget(i, 6).setCheckState(0)

    # hightlights the sequences found in the gene viewer
    # highlighting should stay the exact same with fasta and genbank files, as this function only edits what
    #   is currently in the gene viewer text table anyways
    def highlight_gene_viewer(self):
        # make sure gene viewer is enabled
        if not self.displayGeneViewer.isChecked():
            QtWidgets.QMessageBox.question(self, "Gene Viewer Error",
                                           "Gene Viewer display is off! "
                                           "Please turn the Gene Viewer on in order to highlight the sequences selected",
                                           QtWidgets.QMessageBox.Ok)
            return

        # variables needed
        k = Kegg()
        cursor = self.geneViewer.textCursor()
        format = QtGui.QTextCharFormat()
        noMatchString = ""

        # reset the gene viewer text
        self.geneViewer.setText(self.geneNTDict[self.comboBoxGene.currentText()])

        # check and make sure still is actually highlighted!
        selectedList = self.targetTable.selectedItems()
        if len(selectedList) <= 0:
            QtWidgets.QMessageBox.question(self, "Nothing Selected",
                                           "No targets were highlighted "
                                           "Please highlight the targets you want to be highlighted in the gene viewer!",
                                           QtWidgets.QMessageBox.Ok)
            return
        # this is the loop that actually goes in and highlights all the things
        for i in range(self.targetTable.rowCount()):
            if self.targetTable.item(i, 0).isSelected():
                # get the strand and sequence strings
                locationString = self.targetTable.item(i,0).text()
                strandString = self.targetTable.item(i, 2).text()
                sequenceString = self.targetTable.item(i, 1).text()
                printSequence = ""
                movementIndex = 0
                left_right = ""
                #print("Length of geneViewer: ", len(self.geneViewer.toPlainText()))


                # get the location
                location = int(locationString) - self.geneDict[self.comboBoxGene.currentText()][1]

                # get which way it's moving, and the real location. This is for checking edge cases
                if "Cas12" in self.endonucleaseBox.currentText():
                    # movement is always 24
                    movementIndex = 24

                    # if the strand is positive, it moves to the right, if the strand is negative, it moves to the left
                    if strandString == "-":
                        left_right = "-"
                        location = (location - len(self.targetTable.item(i, 3).text())) + 1
                    elif strandString == "+":
                        location = (location + len(self.targetTable.item(i,3).text())) + 1
                        left_right = "+"
                # if the endo is Cas9
                elif "Cas9" in self.endonucleaseBox.currentText():
                    # movement is always 20
                    movementIndex = 20
                    location = location + 2

                    # if the strand is negative, it moves to the right if the strand is positive it moves to the left
                    if strandString == "-":
                        left_right = "+"
                        #location = location + len(self.targetTable.item(i, 3).text())
                    elif strandString == "+":
                        left_right = "-"
                        #location = location - len(self.targetTable.item(i, 3).text())

                # get the right color and the revcom
                if strandString == "+":
                    format.setBackground(QtGui.QBrush(QtGui.QColor("green")))
                elif strandString == "-":
                    format.setBackground(QtGui.QBrush(QtGui.QColor("red")))
                    sequenceString = k.revcom(sequenceString)

                testSequence = sequenceString.lower()
                testGeneViewer = self.geneViewer.toPlainText().lower()

                #print("Location is: ", location)
                #print("Length of geneome viewer: ", len(self.geneViewer.toPlainText()))

                # check to see if the sequence is in the gene viewer to behind with
                if testSequence in testGeneViewer:
                    #print("In the if testSequence in testGeneViewer")
                    indexInViewer = testGeneViewer.find(testSequence)
                    cursor.setPosition(indexInViewer)
                    for i in range(len(sequenceString)):
                        cursor.movePosition(QtGui.QTextCursor.NextCharacter, 1)
                    cursor.mergeCharFormat(format)

                # below elif's are edge cases, in case the sequence is not fully in the gene viewer
                # if the start is too far to the left, but part of the sequence is in gene viewer
                # and it's being built right-to-left
                elif left_right == "-" and location - movementIndex < 0 and location > 0:
                    #print("in the first elif")
                    cursor.setPosition(0)
                    for i in range(location + 1):
                        cursor.movePosition(QtGui.QTextCursor.NextCharacter, 1)
                    cursor.mergeCharFormat(format)
                # being built left-to-right
                # if start is too far to the left, but start + total movement is in the geneviewer
                elif left_right == "+" and location < 0 and location + movementIndex > 0:
                    #print("In the second elif")
                    cursor.setPosition(0)
                    for i in range(location + movementIndex):
                        cursor.movePosition(QtGui.QTextCursor.NextCharacter, 1)
                    cursor.mergeCharFormat(format)
                # being build right-to-left
                # if the location is too far to the right, but location-movement is in the gene viewer
                elif left_right == "-" and location > len(self.geneViewer.toPlainText()) and location - movementIndex < len(self.geneViewer.toPlainText()):
                    #print("In the third elif statement")
                    cursor.setPosition((location - movementIndex) + 1)
                    cursor.movePosition(QtGui.QTextCursor.End, 1)
                    cursor.mergeCharFormat(format)
                # being built left-to-right
                # if the location + movement is too far to the right, but location is in the geneviewer
                elif left_right == "+" and location + movementIndex > len(self.geneViewer.toPlainText()) and location < len(self.geneViewer.toPlainText()):
                    #print("In the fourth elif")
                    cursor.setPosition(location)
                    cursor.movePosition(QtGui.QTextCursor.End, 1)
                    cursor.mergeCharFormat(format)

                # else, it is not able to be found
                else:
                    if noMatchString == "":
                        noMatchString = sequenceString
                    else:
                        noMatchString = noMatchString + ";;" + sequenceString


        # if any of the sequences return 0 matches, show the user which ones were not found
        if len(noMatchString) >= 5:
            QtWidgets.QMessageBox.question(self, "Warning", "The following sequence(s) were not found in the Gene Viewer"
                                                                "text:\n\t" + noMatchString, QtWidgets.QMessageBox.Ok)

    # this function updates the gene viewer based on the user clicking 'display on'
    # if it is check marked, it displays the correct data
    # if it is un-marked, it hides the data
    def checkGeneViewer(self):
        if self.displayGeneViewer.isChecked():
            self.lineEditStart.setText(str(self.geneDict[self.comboBoxGene.currentText()][1]))
            self.lineEditEnd.setText(str(self.geneDict[self.comboBoxGene.currentText()][2]))
            self.geneViewer.setText(self.geneNTDict[self.comboBoxGene.currentText()])
        elif not self.displayGeneViewer.isChecked():
            self.lineEditStart.clear()
            self.lineEditEnd.clear()
            self.geneViewer.clear()

    # this function just opens ing when the user clicks the CoTargeting button
    # opened the same way that main opens it
    def open_cotarget(self):
        endo_list = list()

        if self.endonucleaseBox.count() <= 1:
            QtWidgets.QMessageBox.question(self, "Not Enough Endonucleases",
                                           "There are not enough endonucleases with this organism. "
                                           "At least 2 endonucleases are required for this function. "
                                           "Use Analyze New Genome to create CSPR files with other endonucleases.",
                                           QtWidgets.QMessageBox.Ok)
            return

        for i in range(self.endonucleaseBox.count()):
            endo_list.append(self.endonucleaseBox.itemText(i))

        GlobalSettings.mainWindow.CoTargeting.launch(endo_list, GlobalSettings.mainWindow.orgChoice.currentText())

    # this function goes through and calls transfer_data again.
    # Uses data from the mainWindow in Globalsettings, but that's because that info should not change
    # unless the user closes out of the Results window
    def changeEndonuclease(self):
        full_org = str(GlobalSettings.mainWindow.orgChoice.currentText())
        organism = GlobalSettings.mainWindow.shortHand[full_org]

        endoChoice = self.endonucleaseBox.currentText().split("|")

        # make sure the user actually selects a new endonuclease
        if self.endo == endoChoice:
            QtWidgets.QMessageBox.question(self, "Select a different Endonuclease",
                                           "Please be sure to select a different endonuclease!",
                                           QtWidgets.QMessageBox.Ok)
            return

        # enable the cotarget checkbox if needed
        if len(endoChoice) > 1:
            self.cotarget_checkbox.setEnabled(True)
            self.cotarget_checkbox.setChecked(0)
        else:
            self.cotarget_checkbox.setEnabled(False)
            self.cotarget_checkbox.setChecked(0)

        self.transfer_data(organism, endoChoice, GlobalSettings.CSPR_DB, self.geneDict,
                           self.geneNTDict, "")

    # Function that is used to set up the results page.
    # it calls get_targets, which in turn calls display data
    def transfer_data(self, org, endo, path, geneposdict, geneNTSeqDict, fasta):
        # set all of the classes variables
        self.org = org
        self.endo = endo
        self.directory = path
        self.fasta_ref = fasta
        self.comboBoxGene.clear()
        self.AllData.clear()
        self.geneDict = geneposdict
        self.geneNTDict = geneNTSeqDict

        self.highlighted.clear()
        self.detail_output_list.clear()
        self.seq_and_avg_list.clear()
        self.rows_and_seq_list.clear()
        self.OTA.clear()

        for gene in geneposdict:
            detail_output1 = {}
            rows_and_seq2 = {}
            seq_and_avg3 = {}
            self.detail_output_list.append(detail_output1)
            self.seq_and_avg_list.append(seq_and_avg3)
            self.rows_and_seq_list.append(rows_and_seq2)
            self.comboBoxGene.addItem(gene)
            self.get_targets(gene, geneposdict[gene])



        # Enable the combobox to be toggled now that the data is in AllData
        self.comboBoxGene.currentTextChanged.connect(self.displayGeneData)
        self.first_boot = True

    def goBack(self):
        # check and see if they searched on Sequence. If so, delete the temp CSPR file
        if len(self.seq_finder_cspr_file) > 0:
            os.remove(self.seq_finder_cspr_file)
            self.displayGeneViewer.setEnabled(False)
            self.gene_viewer_settings_button.setEnabled(True)
            GlobalSettings.mainWindow.pushButton_ViewTargets.setEnabled(False)
            self.seq_finder_cspr_file = ''
        GlobalSettings.mainWindow.show()
        self.cotarget_checkbox.setChecked(0)
        self.hide()

    # called when the user hits 'gene viewer settings'
    def changeGeneViewerSettings(self):
        GlobalSettings.mainWindow.gene_viewer_settings.show()

    # this is the function that sets up the co-targetting.
    # it is called from the coTargetting class, when the user hits submit
    # myBool is whether or not to change the endoChoice comboBox
    def populate_cotarget_table(self, myBool = True):
        # make a string of the combinitation, separated by commas's
        endoBoxString = ""
        for i in range(len(self.co_target_endo_list)):
            if endoBoxString == "":
                endoBoxString = self.co_target_endo_list[i]
            else:
                endoBoxString = endoBoxString + '|' + self.co_target_endo_list[i]

        # put the new endoChoice at the beginning. THis is the only way i could find to do it
        # get a list of all endo choices, and put the newest at the front
        endoBoxList = list()
        endoBoxList.append(endoBoxString)
        for i in range(self.endonucleaseBox.count()):
            endoBoxList.append(self.endonucleaseBox.itemText(i))

        # clear the current endo choices, and append the new order
        if myBool:
            self.endonucleaseBox.clear()
            for i in range(len(endoBoxList)):
                self.endonucleaseBox.addItem(endoBoxList[i])

        # enable the cotarget checkbox
        self.cotarget_checkbox.setEnabled(True)
        self.cotarget_checkbox.setChecked(0)

        # add it to the endoBox choices, and then call transfer_data
        self.transfer_data(self.org, self.co_target_endo_list, GlobalSettings.CSPR_DB, self.geneDict, self.geneNTDict, "")

    # prep function for the checkbox for cotargeting
    # if the checkbox is checked, just go ahead and displayGeneData
    # if not, call populate_cotarget_table, as a reset to get all of the data there
    def prep_cotarget_checkbox(self):
        if self.cotarget_checkbox.isChecked():
            self.displayGeneData()
        elif not self.cotarget_checkbox.isChecked():
            self.populate_cotarget_table(myBool=False)


    # Function grabs the information from the .cspr file and adds them to the AllData dictionary
    #changed to now call CSPRparser's function. Same function essentially, just cleaned up here
    def get_targets(self, genename, pos_tuple):
        #get the right files
        for endo in self.endo:
            if self.directory.find("/") != -1:
                file = (self.directory+"/" + self.org + "_" + endo + ".cspr")
            else:
                file = (self.directory + "/" + self.org + "_" + endo + ".cspr")

            #create the parser, read the targets store it. then display the GeneData screen
            parser = CSPRparser(file)

            # if genename is not in the dict, make that spot into a list
            if genename not in self.AllData:
                self.AllData[genename] = list()
            # now append parser's data to it
            self.AllData[genename].append(parser.read_targets(genename, pos_tuple, endo))

            # for each list item
            for item in self.AllData[genename]:
                # for each tuple item
                for i in range(len(item)):
                    self.highlighted[item[i][1]] = False
        self.displayGeneData()

        # if the endo choice is greater than 1, call the combine function
        if len(self.endo) > 1:
            self.combine_coTargets(genename)


    ###############################################################################################################
    # Main Function for updating the Table.  Connected to all filter buttons and the Gene toggling of the combobox.
    ###############################################################################################################
    def displayGeneData(self):

        curgene = str(self.comboBoxGene.currentText())  # Gets the current gene
        # Creates the set object from the list of the current gene:
        if curgene=='' or len(self.AllData)<1:
            return

        subset_display = set()

        # set the start and end numbers, as well as set the geneViewer text, if the displayGeneViewer is checked
        if self.displayGeneViewer.isChecked():
            self.lineEditStart.setText(str(self.geneDict[self.comboBoxGene.currentText()][1]))
            self.lineEditEnd.setText(str(self.geneDict[self.comboBoxGene.currentText()][2]))
            self.geneViewer.setText(self.geneNTDict[self.comboBoxGene.currentText()])

        # if this checkBox is checked, remove the single endo
        if self.cotarget_checkbox.isChecked():
            self.remove_single_endo(curgene)

        # Removing all sequences below minimum score and creating the set:
        # for each list item
        for item in self.AllData[curgene]:
            # for each tuple item
            for i in range(len(item)):
                if int(item[i][3]) > int(self.minScoreLine.text()):
                    # Removing all non 5' G sequences:
                    if self.fivegseqCheckBox.isChecked():
                        if item[i][1].startswith("G"):
                            subset_display.add(item[i])
                    else:
                        subset_display.add(item[i])

        self.targetTable.setRowCount(len(subset_display))

        index = 0
        #changed the number items to use setData so that sorting will work correctly
        #because before the numbers were interpretted as strings and not numbers
        for item in subset_display:
            num = int(item[0])
            loc = QtWidgets.QTableWidgetItem()
            loc.setData(QtCore.Qt.EditRole, num)
            seq = QtWidgets.QTableWidgetItem(item[1])
            strand = QtWidgets.QTableWidgetItem(str(item[4]))
            PAM = QtWidgets.QTableWidgetItem(item[2])
            num1 = int(item[3])
            endonuclease = QtWidgets.QTableWidgetItem(item[5])
            score = QtWidgets.QTableWidgetItem()
            score.setData(QtCore.Qt.EditRole, num1)
            self.targetTable.setItem(index, 0, loc)
            self.targetTable.setItem(index, 1, seq)
            self.targetTable.setItem(index, 2, strand)
            self.targetTable.setItem(index, 3, PAM)
            self.targetTable.setItem(index, 4, score)
            self.targetTable.setItem(index, 5, QtWidgets.QTableWidgetItem("--.--"))
            self.targetTable.removeCellWidget(index, 7)
            self.targetTable.setItem(index, 8, endonuclease)
            ckbox = QtWidgets.QCheckBox()
            ckbox.clicked.connect(self.search_gene)
            self.targetTable.setCellWidget(index,6,ckbox)
            if(self.highlighted[str(item[1])] == True):
                ckbox.click()
            if (item[1] in self.seq_and_avg_list[self.comboBoxGene.currentIndex()].keys()):
                OT = QtWidgets.QTableWidgetItem()
                OT.setData(QtCore.Qt.EditRole, self.seq_and_avg_list[self.comboBoxGene.currentIndex()][item[1]])
                self.targetTable.setItem(index, 5, OT)
            if (item[1] in self.detail_output_list[self.comboBoxGene.currentIndex()].keys()):
                details = QtWidgets.QPushButton()
                details.setText("Details")
                details.clicked.connect(self.show_details)
                self.targetTable.setCellWidget(index, 7, details)
            index += 1

        self.targetTable.resizeColumnsToContents()

    # this function is only entered if the user checks the show only cotargeted sequence checkbox
    def remove_single_endo(self, curgene):
        removalDict = dict()
        # go through and figure out which ones need to be shown
        for i in range(len(self.AllData[curgene])):
            for j in range(len(self.AllData[curgene][i])):
                endoData = self.AllData[curgene][i][j][5].split("|")
                if len(endoData) == 1:
                    if i not in removalDict:
                        removalDict[i] = list()

                    removalDict[i].append(j)

        # now go through and delete them. But you have to go in reverse
        for item in removalDict:
            # for the reverse of that list. This is to keep the program from crashing
            # easier than building a new list honestly
            for index in reversed(removalDict[item]):
                self.AllData[curgene][item].pop(index)

    # this function goes through and combines table rows that have the same location and PAM dir
    # it edits the dictionary data itself
    # currently it does not take the PAM direction into account
    # parameter curgene:  the key to which part in the dictionary to look at
    def combine_coTargets(self, curgene):
        deletingDict = dict()
        endoList = list()

        #print(self.AllData[curgene])
        # get the endo list data
        for i in range(len(self.AllData[curgene])):
            # if one of them is empty, just return because co-targeting is useless for that one
            if len(self.AllData[curgene][i]) == 0:
                return
            endoList.append(self.AllData[curgene][i][0][5])


        # for each endoNuclease in the curGene block
        for i in range(len(self.AllData[curgene])):
            endoData1 = endoList[i]
            # for each target in that gene
            for j in range(len(self.AllData[curgene][i])):
                # get first locations endo
                locationData1 = self.AllData[curgene][i][j][0]
                sequenceData1 = self.AllData[curgene][i][j][1]
                pamData1 = self.AllData[curgene][i][j][2]
                scoreData1 = self.AllData[curgene][i][j][3]
                strandData1 = self.AllData[curgene][i][j][4]

                # for each endoNuclease in the curGene block
                for k in range(len(self.AllData[curgene])):
                    # if k == i then we are on the same endo target list, so break out because there can't be any combinations
                    if k == i:
                        break
                    endoData2 = endoList[k]
                    # for each target in that gene
                    for l in range(len(self.AllData[curgene][k])):
                        locationData2 = self.AllData[curgene][k][l][0]
                        strandData2 = self.AllData[curgene][k][l][4]
                        pamData2 = self.AllData[curgene][k][l][2]

                        # check which PAM is longer, and store the longer one. Otherwise, just store the first one
                        if len(pamData1) > len(pamData2):
                            storePam = pamData1
                        elif len(pamData1) < len(pamData2):
                            storePam = pamData2
                        else:
                            storePam = pamData1

                        # get the directions
                        dir1 = self.S.endo_info[endoData1][3]
                        dir2 = self.S.endo_info[endoData2][3]

                        # check if can be combined
                        if locationData1 == locationData2 and endoData1 != endoData2 and endoData2 not in endoData1:
                            if dir1 == dir2 and strandData1 == strandData2:
                                storeEndo = self.AllData[curgene][i][j][5]
                                if endoData2 not in storeEndo:
                                    storeEndo = storeEndo + "|" + endoData2
                                # combine the endo data
                                    self.AllData[curgene][i][j] = (locationData1, sequenceData1, storePam, scoreData1, strandData1, storeEndo)

                                # store which ones to delete
                                if k not in deletingDict:
                                    deletingDict[k] = list()
                                deletingDict[k].append(l)


        # delete the ones that need to be deleted
        for item in deletingDict:
            # go in the reverse of that list. This is to keep the program from crashing.
            # it is easier than building a new list honestly
            for index in reversed(deletingDict[item]):
                self.AllData[curgene][item].pop(index)
    ########################################## END UPDATING FUNCTION #############################################

    def search_gene(self):
        search_trms = []
        checkBox = self.sender()
        index = self.targetTable.indexAt(checkBox.pos())
        # print(index.column(), index.row(), checkBox.isChecked())
        seq = self.targetTable.item(index.row(),1).text()
        self.highlighted[str(seq)] = checkBox.isChecked()

        x=1

    def item_select(self):
        print(self.targetTable.selectedItems())

    def table_sorting(self, logicalIndex):
        self.switcher[logicalIndex] *= -1
        if self.switcher[logicalIndex] == -1:
            self.targetTable.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
        else:
            self.targetTable.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)

        #update the OTA list in case user already has stuff selected, because their rows
        #will most likely change after the sort
        self.OTA.clear()
        for row in range(self.targetTable.rowCount()):
            if (self.targetTable.cellWidget(row, 6).isChecked()):
                self.OTA.append(row)

    #linked to when the user pushes tools->off target analysis
    def Off_Target_Analysis(self):
        #build temp file for offtarget to read from
        f = open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'OffTargetFolder') + '/temp.txt','w+')
        self.OTA.clear()
        for row in range(self.targetTable.rowCount()):
            if(self.targetTable.cellWidget(row,6).isChecked()):
                self.OTA.append(row)
                loc = self.targetTable.item(row, 0).text()
                seq = self.targetTable.item(row,1).text()
                strand = self.targetTable.item(row,2).text()
                pam = self.targetTable.item(row,3).text()
                score = self.targetTable.item(row,4).text()
                loc = self.S.compress(int(loc), 64)
                self.rows_and_seq_list[self.comboBoxGene.currentIndex()][seq] = row
                seq = self.S.compress(seq, 64)
                pam = self.S.compress(pam, 64)
                score = self.S.compress(int(score), 64)
                output = str(loc) + ',' + str(seq) + str(strand) + str(pam) + ',' + score
                f.write(output + '\n')
        f.close()
        #only make off target object if first time, otherwise just
        #reshow the object
        if(self.first_boot == True):
            self.first_boot = False
            self.off_tar_win = OffTarget.OffTarget()
            self.off_tar_win.submitButton.clicked.connect(self.refresh_data)
        self.off_tar_win.show()
        f.close()

    #refresh data is linked to the submit button on the off target analysis UI
    def refresh_data(self):
        #setup filename based on output name given in OffTarget
        filename = self.off_tar_win.output_path


        # if the user hits submit without running thr program, do nothing
        if filename == '':
            return


        for rows in range(0,self.targetTable.rowCount()):
            self.targetTable.removeCellWidget(rows,7)
            self.targetTable.removeCellWidget(rows,5)
            self.targetTable.setItem(rows, 5, QtWidgets.QTableWidgetItem("--.--"))
        self.off_tar_win.hide()
        filename = filename[:len(filename)-1]
        filename = filename[1:]
        filename = filename.replace(r'\\', '/')
        filename = filename.replace('"', '')
        self.files_list.append(filename)
        out_file = open(filename, "r")
        #read the first line : either AVG or DETAILED OUTPUT
        output_type = out_file.readline()
        output_type = output_type.strip('\r\n')
        #length = len(self.OTA)
        #parse based on whether avg or detailed output
        if(output_type == "AVG OUTPUT"):
            for line in out_file:
                #line = out_file.readline()
                line = line.strip('\n')
                if (line != ''):
                    values = line.split(":")
                    seq = str(values[0])
                    row = self.rows_and_seq_list[self.comboBoxGene.currentIndex()][values[0]]
                    OT = QtWidgets.QTableWidgetItem()
                    OT.setData(QtCore.Qt.EditRole, values[1])
                    self.targetTable.setItem(row, 5, OT)
                    self.seq_and_avg_list[self.comboBoxGene.currentIndex()][values[0]] = values[1]
        else:
            details_bool = False
            temp_list = []
            values = []
            for line in out_file:
                #line = out_file.readline()
                line = line.strip('\n')
                if(line.find(':') != -1):
                    if(details_bool == True):
                        self.detail_output_list[self.comboBoxGene.currentIndex()][values[0]] = temp_list
                        details_bool = False
                        temp_list = []
                    values = line.split(":")
                    self.seq_and_avg_list[self.comboBoxGene.currentIndex()][values[0]] = values[1]
                    row = self.rows_and_seq_list[self.comboBoxGene.currentIndex()][values[0]]
                    OT = QtWidgets.QTableWidgetItem()
                    OT.setData(QtCore.Qt.EditRole, values[1])
                    self.targetTable.setItem(row, 5, OT)
                else:
                    details_bool = True
                    temp_list.append(line)
                    details = QtWidgets.QPushButton()
                    details.setText("Details")
                    details.clicked.connect(self.show_details)
                    self.targetTable.setCellWidget(row, 7, details)
            if(details_bool == True):
                self.detail_output_list[self.comboBoxGene.currentIndex()][values[0]] = temp_list

            #print(self.detail_output_list)
            self.targetTable.resizeColumnsToContents()
            out_file.close()

    def show_details(self):
        #create msg box popup for the details | used html to make it easier to style with bold
        button = self.sender()
        index = self.targetTable.indexAt(button.pos())
        msg = QtWidgets.QMessageBox()
        msg.setWindowTitle("Details")
        msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
        key = str(self.targetTable.item(index.row(),1).text())
        temp_str = ''
        for items in self.detail_output_list[self.comboBoxGene.currentIndex()][key]:
            temp_str += items + "<br>"

        chromo_str = "<html><b>Chromsome: Location, Sequence, Strand, PAM, Score:<br></b></html>"
        input_str = self.targetTable.item(index.row(),0).text() + ' , ' + key + ' , ' + \
                    self.targetTable.item(index.row(),2).text() + ' , ' + self.targetTable.item(index.row(),3).text() + \
                    ' , ' + self.targetTable.item(index.row(),4).text() + "<br><br>"
        detail_str = "<html><b>Deatailed Output: Score, Chromsome, Location, Sequence:<br></b></html>"
        msg.setText(chromo_str + input_str + detail_str + temp_str)
        msg.exec()


    # Function for displaying the target in the gene viewer
    """def displayGene(self,fastafile=None, Kegg=False, NCBI=False):
        organism_genome = list()  # list of chromosomes/scaffolds
        if fastafile:
            f = open(fastafile)
            chr_string = str()
            for line in f:
                if not line.startswith(">"):
                    chr_string += line[:-1]
                else:
                    organism_genome.append(chr_string)
                    chr_string = ""
            return organism_genome
        elif Kegg:
            # Get the gene from the Kegg database
        elif NCBI:
            # Get the gene from NCBI database (RefSeq)
        else:
            return "Error: Cannot find reference sequence.  Search Kegg, NCBI, or download a FASTA file to create a genome reference."""""

    # -----------------------------------------------------------------------------------------------------#
    # ---- All Filter functions below ---------------------------------------------------------------------#
    # -----------------------------------------------------------------------------------------------------#
    def update_score_filter(self):
        self.minScoreLine.setText(str(self.scoreSlider.value()))

    #allows user to save what is currently in the table
    def save_data(self):
        filename = QtWidgets.QFileDialog.getSaveFileName(self,
                                      "Enter Text File Name", ".txt",
                                      "Text Document (*.txt)" )

        if (str(filename[0]) != ''):
            f = open(str(filename[0]), "w+")
            for genomes in self.AllData:
                f.write("***"+genomes + "\n")
                for items in self.AllData[genomes]:
                    for i in items:
                        for j in i:
                            f.write(str(j) + '|')

                        f.write(str(self.highlighted[items[1][1]]))
                        f.write("\n")
            f.close()

    #open any saved .txt of previous tables opened
    def open_data(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File')
        if (os.path.isfile(str(filename[0]))):
            f = open(str(filename[0]), "r+")
            temp_str = f.readline()
            if(temp_str.startswith("***") == False):
                msg = QtWidgets.QMessageBox()
                msg.setWindowTitle("Error")
                msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.setText("<font size=4>" + "Not correct file format" + "</font>")
                msg.exec()
                f.close()
            else:
                f.close()
                self.AllData.clear()
                self.highlighted.clear()
                first = 1
                list1 = []
                list2 = []
                s = ""
                self.comboBoxGene.clear()
                if (os.path.isfile(str(filename[0]))):
                    f = open(str(filename[0]), "r+")
                    for line in f:
                        if(line.startswith("***")):
                            if(first == 0):
                                list2.append(list1)
                                self.AllData[s] = list2
                            else:
                                first = 0
                            s = line[3:]
                            s = s.strip('\n')
                            list1 = []
                            self.comboBoxGene.addItem(s)
                        else:
                            temp = line.split("|")
                            h = temp.pop()
                            h = h.strip('\n')
                            self.highlighted[temp[1]] = eval(h)
                            tup = tuple(temp)
                            list1.append(tup)
                    list2.append(list1)
                    self.AllData[s] = list2
                    f.close()
                    self.displayGeneData()

    # this function calls the closingWindow class.
    def closeEvent(self, event):
        # check and see if they searched on Sequence. If so, delete the temp CSPR file
        if len(self.seq_finder_cspr_file) > 0:
            os.remove(self.seq_finder_cspr_file)
            GlobalSettings.mainWindow.pushButton_ViewTargets.setEnabled(False)
            self.seq_finder_cspr_file = ''

        GlobalSettings.mainWindow.closeFunction()
        event.accept()

############################################
# CLASS NAME: geneViewerSettings
# It's essentially a little window where the user can tell the program which file to use for geneviewer
############################################
class geneViewerSettings(QtWidgets.QDialog):
    def __init__(self):
        # Qt init stuff
        super(geneViewerSettings, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "geneViewerSettings.ui"), self)
        self.setWindowTitle("Change Gene Viewer Settings")
        self.setWindowIcon(Qt.QIcon(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "cas9image.png")))

        # button connections
        self.kegg_radio_button.clicked.connect(self.change_file_type)
        self.gbff_radio_button.clicked.connect(self.change_file_type)
        self.fna_radio_button.clicked.connect(self.change_file_type)
        self.browse_button.clicked.connect(self.browseForFile)
        self.cancel_button.clicked.connect(self.cancelFunction)
        self.submit_button.clicked.connect(self.submitFunction)

        # class variables
        self.file_type = ""

    # this function is called when the user changes the file type
    # it just sets a class variable to the type of file selected
    def change_file_type(self):
        if self.kegg_radio_button.isChecked():
            self.file_type = "kegg"
        elif self.gbff_radio_button.isChecked():
            self.file_type = "gbff"
        elif self.fna_radio_button.isChecked():
            self.file_type = "fna"


    # this function is only called when the user selects browse for a file option
    # it opens a window such that the user can search for a file to use for the gene viewer sequence
    def browseForFile(self):
        # return out if the user has selected Kegg
        if self.kegg_radio_button.isChecked():
            return

        # make sure that either GBFF or FNA is checked
        if not self.gbff_radio_button.isChecked() and not self.fna_radio_button.isChecked():
            QtWidgets.QMessageBox.question(self, "Nothing Selected",
                                           "Please select either the GBFF or the FNA radio button.",
                                           QtWidgets.QMessageBox.Ok)
            return

        # open a window so that the user can select a file
        filed = QtWidgets.QFileDialog()
        myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose an Annotation File")

        # make sure they choose the correct type of file
        if self.file_type not in myFile[0]:
            QtWidgets.QMessageBox.question(self, "Wrong type of file selected",
                                           "Please select the same type of file selected in the radio buttons.",
                                           QtWidgets.QMessageBox.Ok)
            self.file_name_edit.setText("")
            return

        # if the file is not empty, then set it
        if (myFile[0] != ""):
            self.file_name_edit.setText(myFile[0])

    # this function is only called when the user clicks on cancel.
    # it resets the text for the file chosen, while also unchecking all of the radio buttons
    # then it hides the window
    def cancelFunction(self):
        self.file_name_edit.setText("Please choose a file!")

        # setting them to False does not seem to work, no idea why
        self.kegg_radio_button.setChecked(False)
        self.gbff_radio_button.setChecked(False)
        self.fna_radio_button.setChecked(False)
        self.hide()

    # this function is only called when the user clicks on submit
    # it will find all of the sequences for all of the genes in the comboGeneBox in the results window
    # It will go based on the lengths stored in the comboGeneBox dictionary
    def submitFunction(self):
        if not self.kegg_radio_button.isChecked() and (self.file_name_edit.displayText() == "" or self.file_name_edit.displayText() == "Please choose a file!"):
            #print("No file chosen")
            return

        sequence = ""

        # for each gene selected from the results window
        for item in GlobalSettings.mainWindow.Results.geneDict:
            if self.file_type == "kegg":
                k = Kegg()
                organism = GlobalSettings.mainWindow.Annotations_Organism.currentText().split(" ")[1]
                nameFull = item.split(" ")
                name = nameFull[len(nameFull) - 1]
                # get kegg's ntsequence and store it
                nt_sequence = k.get_nt_sequence(organism + ":" + name)
                GlobalSettings.mainWindow.Results.geneNTDict[item] = nt_sequence
            if self.file_type == "fna":
                sequence = self.fna_sequence_finder(GlobalSettings.mainWindow.Results.geneDict[item])
                GlobalSettings.mainWindow.Results.geneNTDict[item] = sequence
            if self.file_type == "gbff":
                sequence = self.gbff_sequence_finder(GlobalSettings.mainWindow.Results.geneDict[item])
                GlobalSettings.mainWindow.Results.geneNTDict[item] = sequence


        GlobalSettings.mainWindow.Results.displayGeneViewer.setEnabled(True)
        GlobalSettings.mainWindow.Results.lineEditStart.setEnabled(True)
        GlobalSettings.mainWindow.Results.lineEditEnd.setEnabled(True)
        GlobalSettings.mainWindow.Results.change_start_end_button.setEnabled(True)
        GlobalSettings.mainWindow.Results.displayGeneViewer.setChecked(1)
        GlobalSettings.mainWindow.Results.checkGeneViewer()
        self.hide()

    # this function gets the sequence out of the GBFF file
    # may have indexing issues
    def gbff_sequence_finder(self, location_data):
        # start up the function
        fileStream = open(self.file_name_edit.displayText())
        buffer = fileStream.readline()
        index = 1
        pre_sequence = ""

        # get to the first chromosome's origin
        while True:
            if "ORIGIN" in buffer:
                buffer = fileStream.readline()
                break
            buffer = fileStream.readline()

        # skip all of the data until we are at the chromosome we care about
        while index != location_data[0]:
            if "ORIGIN" in buffer:
                index += 1
            buffer = fileStream.readline()

        # get the entire chromesome into a string
        while "//" not in buffer:
            if "LOCUS" in buffer or buffer == "":
                break
            # replace digits with spaces (if i replace them with nothing the program will crash)
            for i in range(len(buffer)):
                if buffer[i].isdigit():
                    buffer = buffer.replace(buffer[i], " ")

            # replace all of the spaces
            buffer = buffer.replace(" ", "")

            # append it to the stored version of the entire string
            if pre_sequence == "":
                pre_sequence = buffer
            else:
                pre_sequence = pre_sequence + buffer

            buffer = fileStream.readline()

        # take out the endlines and uppercase the string
        pre_sequence = pre_sequence.replace("\n", "")
        pre_sequence = pre_sequence.upper()
        #print("Length of the pre-sequence: ", len(pre_sequence))

        # we are having an issue here. Sometimes the length of the pre_sequence string is not large enough
        # need to talk to brian to see what could be causing that
        # get the correct location and return
        ret_sequence = pre_sequence[location_data[1]:location_data[2]]
        return ret_sequence

    # this function is the function that actually finds the sequence
    # May have indexing issues here
    def fna_sequence_finder(self, location_data):
        # Open the file and set the index to 0
        fileStream = open(self.file_name_edit.displayText())
        index = 1

        # skip the file until we get to the chromesome we want
        buffer = fileStream.readline()
        while index != location_data[0]:
            buffer = fileStream.readline()
            if buffer.startswith(">"):
                index += 1

        buffer = fileStream.readline()

        # now go through and get the whole chromesome
        sequence = buffer
        while ">" not in buffer:
            buffer = fileStream.readline()
            if not buffer.startswith(">"):
                sequence = sequence + buffer
                # make sure to break out if the end of file is reached
            if buffer == "":
                break

        # uppercase that chromesome, and change all endlines with spaces
        sequence = sequence.upper()
        sequence = sequence.replace("\n", "")
        #print("Length of the pre-sequence: ", len(sequence))
        # now set the return sequence to the substring that we want
        NTSequence = sequence[location_data[1]:location_data[2]]

        return NTSequence






# Window opening and GUI launching code for debugging #
# ----------------------------------------------------------------------------------------------------- #
"""
app = Qt.QApplication(sys.argv)
app.setOrganizationName("TrinhLab-UTK")
app.setApplicationName("CASPER")
window = Results()
window.transfer_data("yli", "spCas9", "/Users/brianmendoza/Dropbox/CrisprDB/", {"phos.carboxylase":(2,3030460,3032157)}, "/Volumes/Seagate_Drive/FASTAs/yli.fna")
sys.exit(app.exec_())"""
