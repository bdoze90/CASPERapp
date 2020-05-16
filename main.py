import sys
import os, platform
import io
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
from APIs import Kegg, SeqFromFasta
from bioservices import KEGG
from Bio import Entrez
from CoTargeting import CoTargeting
from closingWin import closingWindow

from Results import Results, geneViewerSettings
from NewGenome import NewGenome, NCBI_Search_File
from NewEndonuclease import NewEndonuclease
from kishanAnnotation import createGraph
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QDir, QUrl

import webbrowser
import requests
import GlobalSettings
from bs4 import BeautifulSoup
import multitargeting
from AnnotationParser import Annotation_Parser
from NCBI_API import Assembly
from export_to_csv import export_csv_window
from cspr_chromesome_selection import cspr_chromesome_selection
from generateLib import genLibrary
from functools import partial
############################## MT Libraries #####################
import operator
from PyQt5.QtChart import (QBarCategoryAxis,QBarSet, QChartView, QBarSeries,QChart,QLineSeries)
from Algorithms import SeqTranslate
from CSPRparser import CSPRparser
import populationAnalysis
############################## MT Libraries #####################




# =========================================================================================
# CLASS NAME: AnnotationsWindow
# Inputs: Greg: fill in this information
# Outputs: Greg: fill in this information
# =========================================================================================

class AnnotationsWindow(QtWidgets.QMainWindow):
    def __init__(self, info_path):
        super(AnnotationsWindow, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'Annotation Details.ui'), self)
        self.setWindowIcon(QtGui.QIcon(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "cas9image.png")))
        self.Submit_button.clicked.connect(self.submit)
        self.Go_Back_Button.clicked.connect(self.go_Back)
        self.select_all_checkbox.stateChanged.connect(self.select_all_genes)
        self.mainWindow=""
        self.type = ""

    def submit(self):
        if self.type == "kegg":
            self.mainWindow.collect_table_data()
            self.hide()
            self.mainWindow.show()
        elif self.type == "nonkegg":
            self.mainWindow.collect_table_data_nonkegg()
            self.hide()
            self.mainWindow.show()

    def go_Back(self):
        self.tableWidget.clear()
        self.mainWindow.checkBoxes.clear()
        self.mainWindow.searches.clear()
        self.tableWidget.setColumnCount(0)
        self.mainWindow.show()
        self.mainWindow.progressBar.setValue(0)
        self.hide()

    # this function is very similar to the other fill_table, it just works with the other types of annotation files
    def fill_table_nonKegg(self, mainWindow):
        self.tableWidget.clearContents()
        self.mainWindow = mainWindow
        index = 0
        self.tableWidget.setColumnCount(5)
        self.mainWindow.progressBar.setValue(85)
        self.tableWidget.setHorizontalHeaderLabels("Description;Chromosome #;Type;Gene ID;Select".split(";"))
        mainWindow.checkBoxes = []
        self.type = "nonkegg"

        # below chain of loops goes through and figures out how many rows are needed
        #for searchValue in mainWindow.searches:
        #    for definition in mainWindow.searches[searchValue]:
        #        for gene in mainWindow.searches[searchValue][definition]:
        #            index += 1
        #self.tableWidget.setRowCount(index)

        index = 0
        for searchValue in mainWindow.searches:
            for definition in mainWindow.searches[searchValue]:
                for gene in mainWindow.searches[searchValue][definition]:
                    if(gene[2] == 'gene' or gene[2] == 'tRNA'):
                        self.tableWidget.setRowCount(index + 1)
                        # set the checkbox
                        ckbox = QtWidgets.QCheckBox()
                        self.tableWidget.setCellWidget(index, 4, ckbox)

                        # set the description part of the window as well as set the correct data for the checkbox
                        if definition != gene[0]:
                            defin_obj = QtWidgets.QTableWidgetItem(definition)
                            self.tableWidget.setItem(index, 0, defin_obj)
                            mainWindow.checkBoxes.append([definition])
                        else:
                            checkValue = searchValue.upper()
                            defin_obj = QtWidgets.QTableWidgetItem(checkValue)
                            self.tableWidget.setItem(index, 0 ,defin_obj)
                            mainWindow.checkBoxes.append([checkValue])
                        mainWindow.checkBoxes[len(mainWindow.checkBoxes) - 1].append(gene)
                        mainWindow.checkBoxes[len(mainWindow.checkBoxes) - 1].append(ckbox)

                        # set the type in the window
                        type_obj = QtWidgets.QTableWidgetItem(gene[2])
                        self.tableWidget.setItem(index, 2, type_obj)

                        # set the gene id in the window
                        gene_id_obj = QtWidgets.QTableWidgetItem(gene[0])
                        self.tableWidget.setItem(index, 3, gene_id_obj)

                        chrom_number = QtWidgets.QTableWidgetItem(str(gene[1]))
                        self.tableWidget.setItem(index, 1, chrom_number)

                        index += 1
        index = 0
        self.tableWidget.resizeColumnsToContents()

        # if show all is checked, show the window so the user can select the genes they want
        if mainWindow.Show_All_Results.isChecked():
            mainWindow.hide()
            self.show()
        else: # show all not checked
            if (len(mainWindow.checkBoxes) > 15):  # check the size, throw an error if it is too large
                error = QtWidgets.QMessageBox.question(self, "Large File Found",
                                                       "This Annotation file and search parameter yieled many matches, and could cause a slow down.\n\n"
                                                       "Do you wish to continue?",
                                                       QtWidgets.QMessageBox.Yes |
                                                       QtWidgets.QMessageBox.No,
                                                       QtWidgets.QMessageBox.No)
                if (error == QtWidgets.QMessageBox.No):
                    return -2
            for obj in mainWindow.checkBoxes:  # check every match
                obj[2].setChecked(True)
            self.mainWindow.collect_table_data_nonkegg()
        return 0

    # this is the connection for the select all checkbox
    # selects/deselects all the genes in the table
    def select_all_genes(self):
        # check to see if we're selecting all of them or not
        if self.select_all_checkbox.isChecked():
            select_all = True
        else:
            select_all = False

        # go through and do the selection
        if self.type == 'kegg':
            for i in range(self.tableWidget.rowCount()):
                self.tableWidget.cellWidget(i, 2).setChecked(select_all)
        elif self.type == 'nonkegg':
            for i in range(self.tableWidget.rowCount()):
                self.tableWidget.cellWidget(i,4).setChecked(select_all)


    def fill_Table(self,mainWindow):
        self.tableWidget.clearContents()
        self.mainWindow = mainWindow
        index = 0
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setColumnCount(3)
        self.mainWindow.progressBar.setValue(25)
        self.tableWidget.setHorizontalHeaderLabels("Description;Gene ID;Select".split(";"))
        self.type = "kegg"

        mainWindow.checkBoxes = []
        for sValues in mainWindow.searches:
            for definition in mainWindow.searches[sValues]:
                for gene in mainWindow.searches[sValues][definition]:
                    index = index+1
        self.tableWidget.setRowCount(index)
        if index == 0:
            return -1
        index= 0

        for sValues in mainWindow.searches:
            for definition in mainWindow.searches[sValues]:
                defin_obj = QtWidgets.QTableWidgetItem(definition)
                self.tableWidget.setItem(index,0,defin_obj)
                for gene in mainWindow.searches[sValues][definition]:
                    ckbox = QtWidgets.QCheckBox()
                    mainWindow.checkBoxes.append([definition + " " + gene])
                    mainWindow.checkBoxes[len(mainWindow.checkBoxes)-1].append(ckbox)
                    gene_obj = QtWidgets.QTableWidgetItem(gene)
                    self.tableWidget.setItem(index,1 , gene_obj)
                    self.tableWidget.setCellWidget(index, 2,ckbox)
                    index = index+1
        self.tableWidget.resizeColumnsToContents()
        self.mainWindow.progressBar.setValue(50)
        if mainWindow.Show_All_Results.isChecked():
            mainWindow.hide()
            self.show()
        else: #Show all not checked
            if(len(mainWindow.checkBoxes) > 15):#check the size, throw an error if it is too large
                error = QtWidgets.QMessageBox.question(self, "Large File Found",
                                                       "This Annotation file and search parameter yieled many matches, and could cause a slow down.\n\n"
                                                       "Do you wish to continue?",
                                                       QtWidgets.QMessageBox.Yes |
                                                       QtWidgets.QMessageBox.No,
                                                       QtWidgets.QMessageBox.No)
                if(error == QtWidgets.QMessageBox.No):
                    return -2
            self.mainWindow.progressBar.setValue(65)
            for obj in mainWindow.checkBoxes: #check every match
                obj[1].setChecked(True)
            self.mainWindow.collect_table_data() #collect the data
        return 0

    # this function calls the closingWindow class.
    def closeEvent(self, event):
        GlobalSettings.mainWindow.closeFunction()
        event.accept()


# =========================================================================================
# CLASS NAME: CMainWindow
# Inputs: Takes in the path information from the startup window and also all input parameters
# that define the search for targets e.g. endonuclease, organism genome, gene target etc.
# Outputs: The results of the target search process by generating a new Results window
# =========================================================================================


class CMainWindow(QtWidgets.QMainWindow):

    def __init__(self,info_path):

        super(CMainWindow, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'CASPER_main.ui'), self)
        self.dbpath = ""
        self.info_path = info_path
        self.data = {} # each org genome name and the endonucleases along with it
        self.TNumbers = {} # the T numbers from a kegg search
        self.shortHand ={} #  each org's short name IE bacillus subtillis is bsu
        self.orgcodes = {}  # Stores the Kegg organism code by the format {full name : organism code}
        self.gene_list = {} # list of genes (no ides what they pertain to
        self.searches = {}
        self.checkBoxes = []
        self.add_orgo = []
        self.checked_info = {}
        self.check_ntseq_info = {} # the ntsequences that go along with the checked_info
        self.annotation_parser = Annotation_Parser()
        self.ncbi_searcher = Assembly()
        self.link_list = list() # the list of the downloadable links from the NCBI search
        self.organismDict = dict() # the dictionary for the links to download. Key is the description of the organism, value is the ID that can be found in link_list

        # --- Button Modifications --- #
        self.setWindowIcon(QtGui.QIcon(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "cas9image.png")))
        self.pushButton_FindTargets.clicked.connect(self.gather_settings)
        self.pushButton_ViewTargets.clicked.connect(self.view_results)
        self.pushButton_ViewTargets.setEnabled(False)
        self.radioButton_Gene.clicked.connect(self.toggle_annotation)
        self.radioButton_Position.clicked.connect(self.toggle_annotation)
        self.radioButton_Sequence.clicked.connect(self.toggle_annotation)
        self.Search_Button.clicked.connect(self.search_kegg_ncbi_browse_own)
        self.Annotation_Kegg.clicked.connect(self.change_annotation)
        self.Annotation_Ownfile.clicked.connect(self.change_annotation)
        self.NCBI_Select.clicked.connect(self.change_annotation)
        self.actionUpload_New_Genome.triggered.connect(self.launch_newGenome)
        self.actionUpload_New_Endonuclease.triggered.connect(self.launch_newEndonuclease)
        self.actionOpen_Genome_Browser.triggered.connect(self.launch_newGenomeBrowser)
        self.Add_Orgo_Button.clicked.connect(self.add_Orgo)
        self.Remove_Organism_Button.clicked.connect(self.remove_Orgo)
        self.endoChoice.currentIndexChanged.connect(self.endo_Changed)
        self.GenerateLibrary.clicked.connect(self.prep_genlib)
        self.actionExit.triggered.connect(self.close_app)
        self.visit_repo.triggered.connect(self.visit_repo_func)


        self.Search_Input.setEnabled(False)

        self.change_annotation()
        #self.test_button.clicked.connect(self.Ann_Window)

        self.progressBar.setMinimum(0)
        self.progressBar.setMaximum(100)
        self.progressBar.reset()
        self.Annotation_Window = AnnotationsWindow(info_path)

        #Hide Added orgo boxes
        self.Added_Org_Combo.hide()
        self.Remove_Organism_Button.hide()
        self.Added_Org_Label.hide()
        # --- Menubar commands --- #
        self.actionChange_Directory.triggered.connect(self.change_directory)
        self.actionMultitargeting.triggered.connect(self.changeto_multitargeting)
        self.actionPopulation_Analysis.triggered.connect(self.changeto_population_Analysis)
        self.actionNCBI.triggered.connect(self.open_ncbi_web_page)
        self.actionCasper2.triggered.connect(self.open_casper2_web_page)
        self.actionNCBI_BLAST.triggered.connect(self.open_ncbi_blast_web_page)

        self.Question_Button_add_org.clicked.connect(self.add_org_popup)
        self.Confused_Button.clicked.connect(self.annotation_information)

	# --- Setup for Gene Entry Field --- #
        self.geneEntryField.setPlainText("Example Inputs: \n"
                                               "Gene (LocusID): YOL086C  *for Saccharomyces Cerevisiae ADH1 gene* \n"
                                               "Position: chromosome,start,stop\n chromosome,start,stop...\n"
                                               "Sequence: *Pure sequence. CASPER will search for targets and report off"
                                               "targets based on the genome selected, if any*")

        #self.Kegg_Search_Imput.setPlainText("test")
        #show functionalities on window
        ############################self.view_my_results = Results()
        self.newGenome = NewGenome(info_path)
        self.newEndonuclease = NewEndonuclease()
        self.ncbi_search_dialog = NCBI_Search_File()
        self.CoTargeting = CoTargeting(info_path)
        self.Results = Results()
        self.gene_viewer_settings = geneViewerSettings()
        self.export_csv_window = export_csv_window()
        self.cspr_selector = cspr_chromesome_selection()
        self.genLib = genLibrary()
        self.myClosingWindow = closingWindow()


    def endo_Changed(self):
        i=3
        self.add_orgo.clear()
        self.Add_Orgo_Combo.clear()
        self.Added_Org_Combo.clear()
        self.addOrgoCombo()
        self.Added_Org_Combo.hide()
        self.Added_Org_Label.hide()
        self.Remove_Organism_Button.hide()

    ####---FUNCTIONS TO RUN EACH BUTTON---####
    def remove_Orgo(self):
        self.add_orgo.remove(self.Added_Org_Combo.currentText())
        self.Add_Orgo_Combo.addItem(self.Added_Org_Combo.currentText())
        self.Added_Org_Combo.removeItem(self.Added_Org_Combo.currentIndex())
        if len(self.add_orgo)==0:
            self.Added_Org_Combo.hide()
            self.Added_Org_Label.hide()
            self.Remove_Organism_Button.hide()

    def add_Orgo(self):
        if self.Add_Orgo_Combo.currentText() == "Select Organism":
            QtWidgets.QMessageBox.question(self, "Must Select Organism",
                                           "You must select an Organism to add",
                                           QtWidgets.QMessageBox.Ok)
            return

        self.add_orgo.append(self.Add_Orgo_Combo.currentText())
        self.Added_Org_Combo.addItem(self.Add_Orgo_Combo.currentText())
        self.Add_Orgo_Combo.removeItem(self.Add_Orgo_Combo.currentIndex())
        self.Added_Org_Combo.show()
        self.Added_Org_Label.show()
        self.Remove_Organism_Button.show()


    # this function prepares everything for the generate library function
    # it is very similar to the gather settings, how ever it stores the data instead of calling the Annotation Window class
    # it moves the data onto the generateLib function, and then opens that window
    def prep_genlib(self):
        # make sure the user actually inputs something
        inputstring = str(self.geneEntryField.toPlainText())
        if (inputstring.startswith("Example Inputs:") or inputstring == ""):
            QtWidgets.QMessageBox.question(self, "Error",
                                           "No Gene has been entered. Please enter a gene.",
                                           QtWidgets.QMessageBox.Ok)
        else:
            # standardize the input
            inputstring = inputstring.lower()
            found_matches_bool = True

            # call the respective function
            self.progressBar.setValue(10)
            if self.radioButton_Gene.isChecked():
                ginput = inputstring.split(',')
                found_matches_bool = self.run_results("gene", ginput, openAnnoWindow=False)
            elif self.radioButton_Position.isChecked() or self.radioButton_Sequence.isChecked():
                QtWidgets.QMessageBox.question(self, "Error", "Generate Library can only work with Gene (Locus ID)", QtWidgets.QMessageBox.Ok)
                return
            """
            elif self.radioButton_Position.isChecked():
                pinput = inputstring.split(';')
                found_matches_bool = self.run_results("position", pinput,openAnnoWindow=False)
            elif self.radioButton_Sequence.isChecked():
                sinput = inputstring
                found_matches_bool = self.run_results("sequence", sinput, openAnnoWindow=False)
            """
            # if matches are found
            if found_matches_bool == True:
                # get the cspr file name
                fileName = self.shortHand[self.orgChoice.currentText()]
                fileName = GlobalSettings.CSPR_DB + '/' + fileName + '_' + self.endoChoice.currentText() + '.cspr'

                # get whether its Kegg or not
                kegg_non = ''
                if self.Annotation_Kegg.isChecked():
                    kegg_non = 'kegg'
                else:
                    kegg_non = 'non_kegg'


                # launch generateLib
                self.progressBar.setValue(100)

                # calculate the total number of matches found
                tempSum = 0
                for item in self.searches:
                    tempSum += len(self.searches[item])

                # warn the user if the number is greater than 50
                if tempSum > 50:
                    error = QtWidgets.QMessageBox.question(self, "Many matches Found",
                                                       "More than 50 matches have been found. Continuing could cause a slow down. .\n\n"
                                                       "Do you wish to continue?",
                                                       QtWidgets.QMessageBox.Yes |
                                                       QtWidgets.QMessageBox.No,
                                                       QtWidgets.QMessageBox.No)
                    if(error == QtWidgets.QMessageBox.No):
                        self.searches.clear()
                        self.progressBar.setValue(0)
                        return -2

                self.genLib.launch(self.searches, fileName, kegg_non)
            else:
                self.progressBar.setValue(0)


    # Function for collecting the settings from the input field and transferring them to run_results
    def gather_settings(self):
        inputstring = str(self.geneEntryField.toPlainText())

        #Error check: make sure the user actually inputs something
        if(inputstring.startswith("Example Inputs:") or inputstring == ""):
            QtWidgets.QMessageBox.question(self, "Error",
                "No Gene has been entered. Please enter a gene." ,
                QtWidgets.QMessageBox.Ok)
        else:
            #standardize the input
            inputstring = inputstring.lower()

            self.progressBar.setValue(10)
            if self.radioButton_Gene.isChecked():
                ginput = inputstring.split(',')
                self.run_results("gene", ginput)
            elif self.radioButton_Position.isChecked():
                self.run_results("position", inputstring)
            elif self.radioButton_Sequence.isChecked():
                sinput = inputstring
                self.run_results("sequence", sinput)


# ---- Following functions are for running the auxillary algorithms and windows ---- #
    # this function is parses the annotation file given, and then goes through and goes onto results
    # it will call other versions of collect_table_data and fill_table that work with these file types
    # this function should work with the any type of annotation file, besides kegg.
    # this assumes that the parsers all store the data the same way, which gff and feature table do
    # please make sure the gbff parser stores the data in the same way
    # so far the gff files seems to all be different. Need to think about how we want to parse it
    def run_results_own_ncbi_file(self, inputstring, fileName, openAnnoWindow=True):
        self.annotation_parser = Annotation_Parser()
        self.annotation_parser.annotationFileName = fileName
        fileType = self.annotation_parser.find_which_file_version()

        # if the parser retuns the 'wrong file type' error
        if fileType == -1:
            QtWidgets.QMessageBox.question(self, "Error:",
                                           "We cannot parse the file type given. Please make sure to choose a GBFF, GFF, or Feature Table file."
                                           , QtWidgets.QMessageBox.Ok)
            self.progressBar.setValue(0)
            return


        self.progressBar.setValue(60)

        # this bit may not be needed here. Just a quick error check to make sure the chromesome numbers match
        full_org = str(self.orgChoice.currentText())
        cspr_file = (GlobalSettings.CSPR_DB + "/" + self.shortHand[full_org] + "_" + str(self.endoChoice.currentText()) + ".cspr")
        own_cspr_parser = CSPRparser(cspr_file)
        own_cspr_parser.read_first_lines()

        if len(own_cspr_parser.karystatsList) != self.annotation_parser.max_chrom:
            QtWidgets.QMessageBox.question(self, "Warning:",
                                           "The number of chromesomes do not match. This could cause errors."
                                           , QtWidgets.QMessageBox.Ok)

        ##may need to open a message window here
        ###################

        # now go through and search for the actual locus tag, in the case the user input that
        searchValues = self.separate_line(inputstring[0])
        self.searches.clear()
        for search in searchValues:
            search = self.removeWhiteSpace(search)
            if len(search) == 0:
                continue

            # set the searche's dict of search equal to a new dictionary
            self.searches[search] = {}

            # upper case it, because the files seem to be all uppercase for this part
            checkNormalDict = search.upper()
            if checkNormalDict in self.annotation_parser.reg_dict:  # if it is in the normal dictionary
                for item in self.annotation_parser.reg_dict[checkNormalDict]:  # for each list item in that position
                    if item[0] not in self.searches[search]:  # if its not in the search's position yet
                        self.searches[search][item[0]] = self.annotation_parser.reg_dict[checkNormalDict]
                    elif item not in self.searches[search][
                        item[0]]:  # assume it is in the searches position, but do not store duplicates
                        self.searches[search][item[0]].append(self.annotation_parser.reg_dict[checkNormalDict])
        if len(self.searches[searchValues[0]]) >= 1:  # if the previous search yielded results, do not continue
            if openAnnoWindow:
                self.Annotation_Window.fill_table_nonKegg(self)
                return
            else:
                return True

        self.progressBar.setValue(75)
        # reset, and search the parallel dictionary now
        self.searches = {}
        for search in searchValues:
            search = self.removeWhiteSpace(search)
            if len(search) == 0:
                continue

            self.searches[search] = {}
            for item in self.annotation_parser.para_dict:
                checkingItem = item.lower()  # lowercase now, to match the user's input

                if search in checkingItem:  # if what they are searching for is somewhere in that key
                    # loop through each list item now
                    for i in range(len(self.annotation_parser.para_dict[item])):
                        # now loop through the regular dict's position and store the ones we want
                        if self.annotation_parser.para_dict[item][i] != '':
                            for match in self.annotation_parser.reg_dict[self.annotation_parser.para_dict[item][i]]:
                                if item not in self.searches[search]:
                                    self.searches[search][item] = [match]
                                elif item not in self.searches[search][item]:
                                    self.searches[search][item].append(match)
        # if the search returns nothing, throw an error
        if len(self.searches[searchValues[0]]) <= 0:
            QtWidgets.QMessageBox.question(self, "No Matches Found",
                                           "No matches found with that search, please try again",
                                           QtWidgets.QMessageBox.Ok)
            self.progressBar.setValue(0)
            if openAnnoWindow:
                return
            else:
                return False


        # jsut testing as of now
        #for i in self.searches:
         #  print(i)
          # for j in self.searches[i]:
           #     print("\t", j)
            #    for k in self.searches[i][j]:
             #      print("\t\t", k)
        # if we get to this point, that means that the search yieleded results, so fill the table
        self.progressBar.setValue(80)
        # check whether this function call is for Annotation Window, or for generate Lib
        if openAnnoWindow:
            self.Annotation_Window.fill_table_nonKegg(self)
        else:
            return True
    def run_results(self, inputtype, inputstring, openAnnoWindow=True):
        kegginfo = Kegg()
        org = str(self.orgChoice.currentText())
        endo = str(self.endoChoice.currentText())
        progvalue = 15
        self.searches = {}
        self.gene_list = {}
        self.progressBar.setValue(progvalue)

        # set the gene viewer stuff to false so that the user can't open the gene viewer unless they have a
        # FNA/GBFF file
        # Note: if they are using kegg, the enabled's are set to True because it uses the kegg database by default
        self.Results.displayGeneViewer.setEnabled(False)
        self.Results.lineEditStart.setEnabled(False)
        self.Results.lineEditEnd.setEnabled(False)
        self.gene_viewer_settings.kegg_radio_button.setEnabled(False)
        self.Results.change_start_end_button.setEnabled(False)
        self.Results.displayGeneViewer.setChecked(0)


        if inputtype == "gene":
            # make sure an annotation file has been selected
            if self.Annotations_Organism.currentText() == "":
                error = QtWidgets.QMessageBox.question(self, "No Annotation",
                                                   "Please select an Annotation from either KEGG, NCBI, or provide you own Annotation File"
                                                   , QtWidgets.QMessageBox.Ok)
                self.progressBar.setValue(0)
                return

            # ncbi file search code
            if self.NCBI_Select.isChecked():
                type_of_annotation_file = ""
                compressed_file = ""

                # make sure they select a type of file
                if not self.gbff_button.isChecked() and not self.gff_button.isChecked() and not self.feature_table_button.isChecked():
                    QtWidgets.QMessageBox.question(self, "Error",
                                                   "Please select a type of annotation file to download. (Ex: Feature_Table, GFF, GBFF)"
                                                   , QtWidgets.QMessageBox.Ok)
                    self.progressBar.setValue(0)
                    return

                # check to see which file type they selected, and then make sure the program
                # downloads that file
                if self.gbff_button.isChecked():
                    type_of_annotation_file = "_genomic.gbff.gz"
                elif self.gff_button.isChecked():
                    type_of_annotation_file = "_genomic.gff.gz"
                elif self.feature_table_button.isChecked():
                    type_of_annotation_file = "_feature_table.txt.gz"

                # go through and find the link that works, and download that compressed file
                for i in range(len(self.link_list)):
                    if self.organismDict[self.Annotations_Organism.currentText()] in self.link_list[i]:
                        compressed_file = self.ncbi_searcher.download_compressed_annotation_file(self.link_list[i], type_of_annotation_file)
                        break

                #decompress that file, and then delete the compressed version
                if compressed_file:
                    self.progressBar.setValue(25)
                    storeFileName = self.ncbi_searcher.decompress_annotation_file(compressed_file, type_of_annotation_file)
                    file_names = [f for f in os.listdir(GlobalSettings.CSPR_DB) if os.path.isfile(os.path.join(GlobalSettings.CSPR_DB, f))]
                    for file in file_names:
                        if ".gz" in file:
                            #print("Deleting: ", file)
                            os.remove(file)

                    # now run results
                    self.progressBar.setValue(35)
                    myBool = self.run_results_own_ncbi_file(inputstring, storeFileName, openAnnoWindow=openAnnoWindow)
                    if not openAnnoWindow:
                        return myBool
                else:
                    QtWidgets.QMessageBox.question(self, "Error",
                                                   "The database does not have the type of file you have requested. Please try another type of file"
                                                   , QtWidgets.QMessageBox.Ok)
                    self.progressBar.setValue(0)
                    return

            # own annotation file code
            if self.Annotation_Ownfile.isChecked():
                # this now just goes onto the other version of run_results
                myBool = self.run_results_own_ncbi_file(inputstring, self.Annotations_Organism.currentText(), openAnnoWindow=openAnnoWindow)
                if not openAnnoWindow:
                    return myBool
            # KEGG's code
            elif self.Annotation_Kegg.isChecked():
                #check to make sure that both the annotation file and the cspr files have the same version
                #IE both have bsu or otherwise. Just warn the user that the program could crash, or that targets may not be found
                checkList = self.Annotations_Organism.currentText().split(" ")
                if(checkList[1] != self.shortHand[self.orgChoice.currentText()]):
                    error = QtWidgets.QMessageBox.question(self, "Mismatched Annotation File", "The annotation file you have selected does not match the .cspr file selected. Continuing could result in the program crashing. "
                                                                                                 "Targets may not be found as well.\n\n"
                                                                 "Do you wish to continue?",
                                                           QtWidgets.QMessageBox.Yes |
                                                           QtWidgets.QMessageBox.No,
                                                           QtWidgets.QMessageBox.No)

                    if(error == QtWidgets.QMessageBox.No):
                        self.progressBar.setValue(0)
                        return

                self.make_dictonary()
                self.Results.displayGeneViewer.setEnabled(True)
                self.Results.lineEditStart.setEnabled(True)
                self.Results.lineEditEnd.setEnabled(True)
                self.gene_viewer_settings.kegg_radio_button.setEnabled(True)
                self.Results.change_start_end_button.setEnabled(True)
                self.gene_viewer_settings.file_type = "kegg"
                list_sVal = self.separate_line(inputstring[0])
                for sValue in list_sVal:
                    sValue = self.removeWhiteSpace(sValue)
                    if len(sValue) == 0:
                        continue
                    self.searches[sValue] = {}
                    for defin in self.gene_list:
                        #set a string equal to a string version of defin lowercased. That way the case of user input does not matter
                        checkString = str(defin).lower()
                        if sValue in checkString or (sValue in self.gene_list[defin]):
                            if defin in self.searches[sValue]:
                                if self.gene_list[defin] not in self.searches[sValue][defin]:
                                    self.searches[sValue][defin].append(self.gene_list[defin])
                            else:
                                self.searches[sValue][defin] = self.gene_list[defin]
                            #print(self.searches[sValue][defin][0])
                if openAnnoWindow:
                    did_work = self.Annotation_Window.fill_Table(self)
                    if did_work == -1:
                        QtWidgets.QMessageBox.question(self, "Gene Database Error",
                                                       "The Gene you entered could not be found in the Kegg database. "
                                                       "Please make sure you entered everything correctly and try again.",
                                                       QtWidgets.QMessageBox.Ok)
                        self.progressBar.setValue(0)
                    elif did_work == -2: #if the user selects 'no' from the warning of a large file
                        self.progressBar.setValue(0)
                    else:
                        self.progressBar.setValue(100)
                else:
                    # if this function call is for generateLib and not for Annotation Window

                    # get the start/end/chrom data
                    self.searches = self.getKeggDataForGenLib()
                    # see if there's matches
                    if len(self.searches) == 0:
                        QtWidgets.QMessageBox.question(self, "Gene Database Error",
                                                   "The Gene you entered could not be found in the Kegg database. "
                                                   "Please make sure you entered everything correctly and try again.",
                                                   QtWidgets.QMessageBox.Ok)
                        self.progressBar.setValue(0)
                        return False
                    else:
                        return True
                        self.progressBar.setValue(100)
            else:
                self.progressBar.setValue(0)
                return
        # position code below
        if inputtype == "position":
            inputstring = inputstring.replace(' ', '')
            searchInput = inputstring.split('\n')
            full_org = str(self.orgChoice.currentText())
            self.checked_info.clear()
            self.check_ntseq_info.clear()

            for item in searchInput:
                searchIndicies = item.split(',')

                # make sure the right amount of arguments were passed
                if len(searchIndicies) != 3:
                    QtWidgets.QMessageBox.question(self, "Position Error: Invalid Input",
                                                   "There are 3 arguments required for this function. The chromosome, start position, and end position.",
                                                   QtWidgets.QMessageBox.Ok)
                    self.progressBar.setValue(0)
                    return

                # make sure user inputs digits
                if not searchIndicies[0].isdigit() or not searchIndicies[1].isdigit() or not searchIndicies[2].isdigit():
                    QtWidgets.QMessageBox.question(self, "Position Error: Invalid Input",
                                                   "The positions given must be integers. Please try again.",
                                                   QtWidgets.QMessageBox.Ok)
                    self.progressBar.setValue(0)
                    return
                # make sure start is less than end
                elif int(searchIndicies[1]) >= int(searchIndicies[2]):
                    QtWidgets.QMessageBox.question(self, "Position Error: Start must be less than End",
                                                   "The start index must be less than the end index.",
                                                   QtWidgets.QMessageBox.Ok)
                    self.progressBar.setValue(0)
                    return
                # append the data into the checked_info
                tempString = 'chrom: ' + str(searchIndicies[0]) + ' start: ' + str(searchIndicies[1]) + ' end: ' + str(searchIndicies[2])
                self.checked_info[tempString] = (int(searchIndicies[0]), int(searchIndicies[1]), int(searchIndicies[2]))

            self.progressBar.setValue(50)


            self.Results.transfer_data(self.shortHand[full_org], [str(self.endoChoice.currentText())], os.getcwd(),
                                   self.checked_info, self.check_ntseq_info, "")
            self.progressBar.setValue(100)
            self.pushButton_ViewTargets.setEnabled(True)
        # sequence code below
        if inputtype == "sequence":
            checkString = 'AGTCN'
            self.checked_info.clear()
            self.progressBar.setValue(10)
            inputstring = inputstring.upper()

            # check to make sure that the use gave a long enough sequence
            if len(inputstring) < 100:
                QtWidgets.QMessageBox.question(self, "Error", "The sequence given is too small. At least 100 characters is required.", QtWidgets.QMessageBox.Ok)
                self.progressBar.setValue(0)
                return

            # give a warning if the length of the sequence is long
            if len(inputstring) > 30000:
                error = QtWidgets.QMessageBox.question(self, "Large Sequence Detected",
                                                       "The sequence given is a large one, and could slow down the process.\n\n"
                                                       "Do you wish to continue?",
                                                       QtWidgets.QMessageBox.Yes |
                                                       QtWidgets.QMessageBox.No,
                                                       QtWidgets.QMessageBox.No)
                if (error == QtWidgets.QMessageBox.No):
                    self.progressBar.setValue(0)
                    return

            # make sure all the chars are one of A, G, T, C, or N
            for letter in inputstring:
                # skip the end line character
                if letter == '\n':
                    continue
                if letter not in checkString:
                    QtWidgets.QMessageBox.question(self, "Sequence Error",
                                                   "The sequence must consist of A, G, T, C, or N. No other characters are allowed.",
                                                   QtWidgets.QMessageBox.Ok)
                    self.progressBar.setValue(0)
                    return
            self.progressBar.setValue(30)

            # build the CSPR file, and go into results
            fna_file_path = GlobalSettings.CSPR_DB + '/temp.fna'
            self.checked_info['Sequence Finder'] = (1, 0, len(inputstring))
            self.check_ntseq_info['Sequence Finder'] = inputstring.replace('\n', '')
            outFile = open(fna_file_path, 'w')
            outFile.write('>temp org here\n')
            outFile.write(inputstring)
            outFile.write('\n\n')
            outFile.close()
            self.progressBar.setValue(55)
            self.build_cspr_for_seq(fna_file_path)

    # this function actually builds the cspr file
    # I am not pulling out the data for the progress bar here, but it could be added easily. Figured it would be best without it since the progress should run quickly
    def build_cspr_for_seq(self, fna_file_path):
        seq_search_process = QtCore.QProcess()
        # this function decides what to do when the process is finished.
        # what it does:
        #       Deletes the temp FNA file, and calls transfer_data in results.
        def finished():
            # get the file name and process all of the targets
            cspr_file_name = GlobalSettings.CSPR_DB + '/' + code + '_' + myEndoChoice + '.cspr'
            os.remove(fna_file_path)
            seq_search_process.kill()
            self.Results.transfer_data('tempCode', [str(self.endoChoice.currentText())], os.getcwd(),
                                   self.checked_info, self.check_ntseq_info, "")

            self.Results.displayGeneViewer.setEnabled(True)
            self.Results.gene_viewer_settings_button.setEnabled(False)
            self.progressBar.setValue(100)
            self.pushButton_ViewTargets.setEnabled(True)
            self.Results.seq_finder_cspr_file = cspr_file_name

        #--------------------getting all of the arguments-------------------
        myEndoChoice = str(self.endoChoice.currentText())
        my_seq = SeqTranslate()
        pam = my_seq.endo_info[myEndoChoice][0].split(',')[0]
        code = 'tempCode'
        if int(my_seq.endo_info[myEndoChoice][3]) == 3:
            pamdir = False
        else:
            pamdir = True

        output_location = GlobalSettings.CSPR_DB
        path_to_info = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'CASPERinfo')
        orgName = 'temp org'
        gRNA_length = my_seq.endo_info[myEndoChoice][2]
        seed_length = my_seq.endo_info[myEndoChoice][1]
        secondCode = 'test second code'
        #----------------------------------------------------------------------

        #------------compile the executable line------------------------------------
        program = '"' + GlobalSettings.appdir + '/Casper_Seq_Finder_Windows' + '" '
        args =  '"' + myEndoChoice + '" '
        args = args + '"' + pam + '" '
        args = args + '"' + code + '" '
        args = args + str(pamdir) + ' '
        args = args + '"' + output_location + '" '
        args = args + '"' + path_to_info + '" '
        args = args + '"' + fna_file_path + '" '
        args = args + '"' + orgName + '" '
        args = args + gRNA_length + ' '
        args = args + seed_length + ' '
        args = args + '"' + secondCode + '"'
        program = program + args
        #----------------------------------------------------------------------------

        # start the process
        seq_search_process.start(program)
        seq_search_process.finished.connect(finished)

    def launch_newGenome(self):
        self.hide()
        self.newGenome.show()

    def launch_newEndonuclease(self):
       self.newEndonuclease.show()

    def launch_newGenomeBrowser(self):
       print("creating graph")
       createGraph(self)

    # this function does the same stuff that the other collect_table_data does, but works with the other types of files
    def collect_table_data_nonkegg(self):
        # start out the same as the other collect_table_data
        self.checked_info.clear()
        self.check_ntseq_info.clear()
        full_org = str(self.orgChoice.currentText())
        holder = ()


        for item in self.checkBoxes:
            if item[2].isChecked():
                # if they searched base on Locus Tag
                if item[0] in self.annotation_parser.reg_dict:
                    # go through the dictionary, and if they match, store the item in holder
                    for match in self.annotation_parser.reg_dict[item[0]]:
                        if item[1] == match:
                            holder = (match[1], match[3], match[4])
                            self.checked_info[item[0]] = holder
                else:
                    # now we need to go through the para_dict
                    for i in range(len(self.annotation_parser.para_dict[item[0]])):
                        # now go through the matches in the normal dict's data
                        for match in self.annotation_parser.reg_dict[self.annotation_parser.para_dict[item[0]][i]]:
                            # if they match, store it in holder
                            if item[1] == match:
                                holder = (match[1], match[3], match[4])
                                self.checked_info[item[0]] = holder

        # now call transfer data
        self.progressBar.setValue(95)
        self.Results.transfer_data(self.shortHand[full_org], [str(self.endoChoice.currentText())], os.getcwd(),
                                   self.checked_info, self.check_ntseq_info, "")
        self.progressBar.setValue(100)
        self.pushButton_ViewTargets.setEnabled(True)


    # this function is very similar for the Kegg version of collect_table_data
    # difference is that it doesn't use the checkboxes from Annotation Window, just the search values
    def getKeggDataForGenLib(self):
        # variables
        k = Kegg()
        full_org = self.orgChoice.currentText()
        organism = self.Annotations_Organism.currentText().split(" ")[1]
        nameFull = ""
        holder = ()
        temp_dict = dict()

        # for each search value
        for item in self.searches:
            # for each gene in the search results
            for gene in self.searches[item]:
                # name is the ID tag from KEGG
                name = self.searches[item][gene][0]

                gene_info = k.gene_locator(organism+":"+name)

                if gene_info != -1:
                    if gene_info[1] == False:
                        holder = (gene_info[0],'-',gene_info[2],gene_info[3])
                    elif gene_info[1] == True:
                        holder = (gene_info[0],'+',gene_info[2],gene_info[3])
                    # append the data
                    if item not in temp_dict:
                        temp_dict[item] = dict()
                        if gene not in temp_dict[item]:
                            temp_dict[item][gene] = [holder]
                        else:
                            temp_dict[item][gene].append(holder)
                    else:
                        temp_dict[item][gene] = [holder]

        # return this new dict
        return temp_dict


    def collect_table_data(self):
        # need to change this code so that it works with other types of files
        self.checked_info.clear()
        self.check_ntseq_info.clear()

        k = Kegg()
        full_org = str(self.orgChoice.currentText())
        organism= self.Annotations_Organism.currentText().split(" ")[1]
        nameFull = ""
        holder = ()
        for item in self.checkBoxes:
            if item[1].isChecked() ==True:
                nameFull = item[0].split(" ")
                name  = nameFull[len(nameFull)-1]
                gene_info = k.gene_locator(organism+":"+name)
                #print(nameFull)

                # get kegg's ntsequence and store it
                if gene_info != -1:
                    nt_sequence = k.get_nt_sequence(organism+":"+name)

                    #print(item[0])
                    holder = (gene_info[0],gene_info[2],gene_info[3],gene_info[4])
                    self.checked_info[item[0]]=holder
                    self.check_ntseq_info[item[0]] = nt_sequence

        if len(self.checked_info) > 0:
            self.progressBar.setValue(80)
            self.Results.transfer_data(self.shortHand[full_org], [str(self.endoChoice.currentText())] ,os.getcwd(),self.checked_info, self.check_ntseq_info, "")
            self.progressBar.setValue(100)
            self.pushButton_ViewTargets.setEnabled(True)
        else:
            print("No items were found. Please search again.")


# ------------------------------------------------------------------------------------------------------ #

# ----- Following Code is helper functions for processing input data ----- #
    def separate_line(self, input_string):
        export_array = []
        while True:
            index = input_string.find('\n')
            if index==-1:
                if len(input_string)==0:
                    return export_array
                else:
                    export_array.append(input_string)
                    return export_array
            export_array.append(input_string[:index])
            input_string= input_string[index+1:]

    def removeWhiteSpace(self, strng):
        while True:
            if len(strng)==0 or (strng[0] != " " and strng[0]!="\n"):
                break
            strng = strng[1:]
        while True:
            if len(strng)==0 or (strng[len(strng)-1]!=" " and strng[0]!="\n"):
                return strng
            strng = strng[:len(strng)-1]

    # Function to enable and disable the Annotation function if searching by position or sequence
    def toggle_annotation(self):
        if self.radioButton_Gene.isChecked():
            s = True
        else:
            s = False

        # check to see if the sequence button is pressed, and act accordingly
        if self.radioButton_Sequence.isChecked():
            mySeq = SeqTranslate()
            seq_checker = False
            # time to reset the endo's
            self.endoChoice.clear()
            for item in mySeq.endo_info:
                self.endoChoice.addItem(item)
        else:
            seq_checker = True
            self.changeEndos()
        current = self.selected_annotation()
        #print(current)
        if current == "Own":
            self.Search_Button.setText("Browse")
            self.Search_Button.setEnabled(s)
            self.Search_Label.setText("Select an annotation file...")
        elif current == "Kegg":
            self.Search_Button.setText("Search")
            self.Search_Button.setEnabled(s)
            self.Search_Label.setText("Search KEGG Database for genes")
        else:
            self.Search_Button.setText("Search")
            self.Search_Button.setEnabled(s)
            self.Search_Label.setText("Search NCBI Database for genes")
        self.Annotations_Organism.setEnabled(s)
        self.Annotation_Ownfile.setEnabled(s)
        self.Annotation_Kegg.setEnabled(s)
        self.orgChoice.setEnabled(seq_checker)
        self.NCBI_Select.setEnabled(s)

    def selected_annotation(self):
        if self.Annotation_Ownfile.isChecked():
            return "Own"
        elif self.Annotation_Kegg.isChecked():
            return "Kegg"
        else:
            return "NCBI"



    def change_annotation(self):
        if self.Annotation_Ownfile.isChecked():
            self.refseq_button.hide()
            self.genbank_button.hide()
            self.feature_table_button.hide()
            self.gff_button.hide()
            self.gbff_button.hide()
            self.ncbi_ret_max_line_edit.hide()
            self.ncbi_ret_max_label.hide()
            self.Search_Input.setEnabled(False)
            self.Search_Button.setText("Browse")
            self.Search_Label.setText("Select an annotation file...")
        elif self.Annotation_Kegg.isChecked():
            self.refseq_button.hide()
            self.genbank_button.hide()
            self.feature_table_button.hide()
            self.gff_button.hide()
            self.gbff_button.hide()
            self.ncbi_ret_max_line_edit.hide()
            self.ncbi_ret_max_label.hide()
            self.Search_Input.setEnabled(False)
            self.Search_Button.setText("Search")
            self.Search_Label.setText("Search KEGG Database for genome annotation")
            self.Search_Input.setText(self.orgChoice.currentText())
        elif self.NCBI_Select.isChecked():
            self.refseq_button.show()
            self.genbank_button.show()
            self.feature_table_button.show()
            self.gff_button.show()
            self.gbff_button.show()
            self.ncbi_ret_max_line_edit.show()
            self.ncbi_ret_max_label.show()
            self.Search_Input.setEnabled(True)
            self.Search_Button.setText("Search")
            self.Search_Label.setText("Search NCBI Database for genome annotation")
            self.Search_Input.setText(self.orgChoice.currentText())



    # This function works as a way to look up a search term in the Kegg database to potentially get the code
    # for the gene


    def search_kegg_ncbi_browse_own(self):
        # code that lets the user input their own Annotation File
        if self.Annotation_Ownfile.isChecked():
            # have the user choose an annotaiton file, set Search_Input's text to that file
            filed = QtWidgets.QFileDialog()
            myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose an Annotation File")
            if(myFile[0] != ""):
                self.Annotations_Organism.clear()
                self.Search_Input.setText(myFile[0])
                self.Annotations_Organism.addItem(myFile[0])

        # code that uses Kegg
        elif self.Annotation_Kegg.isChecked():
            #make sure user actually inputs something
            if (self.Search_Input.text() == ""):
                QtWidgets.QMessageBox.question(self, "Error",
                                               "Please eneter a search parameter.",
                                               QtWidgets.QMessageBox.Ok)
            else:
                self.Annotations_Organism.clear()
                k = KEGG()

                print("Searching Kegg for: ", self.Search_Input.text())

                #search the kegg database. If the orignal search input returns nothing, split it up on spaces and search each individual term
                All_org = k.lookfor_organism(self.Search_Input.text())
                """ self.lineEdit_search.text()"""
                for item in All_org:
                    hold = self.organism_finder(item)
                    self.Annotations_Organism.addItem(hold)
                    self.TNumbers[hold] = item[:6]
                #if the main string returned nothing, split it up and search each individual term
                if(len(self.Annotations_Organism) <= 0):
                    stringList = self.Search_Input.text().split(" ")
                    for i in range(len(stringList)):
                        print("Searching Kegg for: ", stringList[i])
                        if stringList[i] != "":
                            All_org = k.lookfor_organism(stringList[i])

                            for item in All_org:
                                hold = self.organism_finder(item)
                                #make sure that there are no repeats
                                if(self.Annotations_Organism.findText(hold) == -1):
                                    self.Annotations_Organism.addItem(hold)
                                    self.TNumbers[hold] = item[:6]

                #print("Done searching.\n")
                if(len(self.Annotations_Organism) <= 0):
                    QtWidgets.QMessageBox.question(self, "Error",
                                                  "No matches found with that search parameter",
                                                   QtWidgets.QMessageBox.Ok)

        # code that uses NCBI
        elif self.NCBI_Select.isChecked():
            # error check
            if not self.refseq_button.isChecked() and not self.genbank_button.isChecked():
                QtWidgets.QMessageBox.question(self, "Error", "Please select either RefSeq or GenBank databases.",
                                               QtWidgets.QMessageBox.Ok)
                return

            # if there is nothing set in the line edit, set the ret max to 20
            if self.ncbi_ret_max_line_edit.displayText() != "":

                # if the string is only digits (and non negative), set ret max to that, otherwise just set it to 20
                if self.ncbi_ret_max_line_edit.displayText().isdigit() and int(self.ncbi_ret_max_line_edit.displayText()) > 0:
                    ret_max = int(self.ncbi_ret_max_line_edit.displayText())
                else:
                    ret_max = 20
            else:
                ret_max = 20

            # make sure that the retmax value is not too large
            if ret_max > 100:
                QtWidgets.QMessageBox.question(self, "Error",
                                           "Return Max number is too high, please set it to something below 100",
                                           QtWidgets.QMessageBox.Ok)
                return

            # clear all the things
            self.link_list = list()
            self.organismDict = dict()
            self.Annotations_Organism.clear()

            # now we can finally search NCBI
            database_type = ""
            if self.refseq_button.isChecked():
                database_type ="RefSeq"
            else:
                database_type = "GenBank"

            # actually search, if nothing is returned, break out
            self.link_list, self.organismDict = self.ncbi_searcher.getDataBaseURL(self.Search_Input.displayText(), database_type, ret_max)
            if len(self.link_list) == 0 and len(self.organismDict) == 0:
                QtWidgets.QMessageBox.question(self, "Error", "Search yielded 0 results. Please try again.",
                                               QtWidgets.QMessageBox.Ok)
                return
            # add each item found into the dropdown menu
            for item in self.organismDict:
                print(item)
                self.Annotations_Organism.addItem(item)
           # print("Done searching NCBI")



    def make_dictonary(self):
        url = "https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+gn:"+self.TNumbers[self.Annotations_Organism.currentText()]
        source_code = requests.get(url)
        plain_text = source_code.text
        buf = io.StringIO(plain_text)

        while True:
            line = buf.readline()
            if line[0] == "-":
                break
        while True:
            line = buf.readline()
            if line[1] != "a":
                return
            line = line[line.find(">")+1:]
            seq = line[line.find(":")+1:line.find("<")]
            line = line[line.find(">")+1:]

            i =0
            while True:
                if line[i] == " ":
                    i = i+1
                else:
                    break
            key = line[i:line.find("\n") - 1]
            if key in self.gene_list:
                if seq not in self.gene_list[key]:
                    self.gene_list[key].append(seq)
            else:
                self.gene_list[key] = [seq]
            z=5

    def organism_finder(self, long_str):
        semi = long_str.find(";")
        index =1
        while True:
            if long_str[semi-index] == " ":
                break
            index = index+1
        return long_str[:semi-index]

    def search_own_file(self):
        print("searching for own file")

    # This method is for testing the execution of a button call to make sure the button is linked properly
    def testexe(self):
        choice = QtWidgets.QMessageBox.question(self, "Extract!", "Are you sure you want to Quit?",
                                            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if choice == QtWidgets.QMessageBox.Yes:
            #print(self.orgChoice.currentText())
            sys.exit()
        else:
            pass

    # This method checks if a check button or radio button works appropriately by printing the current organism
    def testcheckandradio(self):
         print(str(self.orgChoice.currentText()))

    def addOrgoCombo(self):
        self.Add_Orgo_Combo.addItem("Select Organism")
        for item in self.data:
            if (self.endoChoice.currentText() in self.data[item]) and (item != str(self.orgChoice.currentText())):
                self.Add_Orgo_Combo.addItem(item)


    # ----- CALLED IN STARTUP WINDOW ------ #
    def getData(self):
        mypath = os.getcwd()
        found = False;
        self.dbpath = mypath
        onlyfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
        orgsandendos = {}
        shortName = {}
        for file in onlyfiles:
            if file.find('.cspr')!=-1:
                found=True;
                newname = file[0:-4]
                s = newname.split('_')
                hold = open(file)
                buf = (hold.readline())
                species = buf[8:buf.find('\n')]
                endo = str(s[1][:len(s[1])-1])
                if species not in shortName:
                    shortName[species] = s[0]
                if species in orgsandendos:
                    orgsandendos[species].append(endo)
                else:
                    orgsandendos[species] = [endo]
                    if self.orgChoice.findText(species) == -1:
                        self.orgChoice.addItem(species)

        #auto fill the kegg search bar with the first choice in orgChoice
        self.Search_Input.setText(self.orgChoice.currentText())
        if found==False:
            return False
        self.data = orgsandendos
        self.shortHand= shortName
        self.endoChoice.clear()
        self.endoChoice.addItems(self.data[str(self.orgChoice.currentText())])
        self.orgChoice.currentIndexChanged.connect(self.changeEndos)


    def changeEndos(self):

        self.endoChoice.clear()
        self.endoChoice.addItems(self.data[str(self.orgChoice.currentText())])
        self.Search_Input.setText(self.orgChoice.currentText())

    def change_directory(self):
        filed = QtWidgets.QFileDialog()
        mydir = QtWidgets.QFileDialog.getExistingDirectory(filed, "Open a Folder",
                                                           self.dbpath, QtWidgets.QFileDialog.ShowDirsOnly)

        if(os.path.isdir(mydir)):
            os.chdir(mydir)
        self.getData()

    #Tanner - added this function to allow the Tools->Multitargeting button to work
    #Function launches the multitargeting window and closing the current one
    def changeto_multitargeting(self):
        os.chdir(os.getcwd())
        GlobalSettings.MTWin.show()
        GlobalSettings.mainWindow.hide()

    def changeto_population_Analysis(self):
        GlobalSettings.pop_Analysis.launch(GlobalSettings.CSPR_DB)
        GlobalSettings.pop_Analysis.show()
        GlobalSettings.mainWindow.hide()

    def add_org_popup(self):
        info = "This functionality will allow users to use different organisms for off-target analysis in a future " \
               "version of the software. If you need to run analysis on multiple organisms, please use the Population " \
               "Analysis feature."
        QtWidgets.QMessageBox.information(self, "Add Organism Information", info, QtWidgets.QMessageBox.Ok)

    def annotation_information(self):
        info = "Annotation files for searching for targets on a gene/locus basis can be selected here using either KEGG " \
               "or NCBI databases, or uploading your own file. Note that KEGG searches are best done with exact matches " \
               "(e.g include strain designation), whereas NCBI will often return multiple assemblies of the same species. " \
               "If you have trouble deciding on the search returned annotation, go to the website, download the annotation " \
               "file to your local computer and choose ‘Choose Annotation File’"
        QtWidgets.QMessageBox.information(self, "Annotation Information", info, QtWidgets.QMessageBox.Ok)

    def open_ncbi_blast_web_page(self):
        webbrowser.open('https://blast.ncbi.nlm.nih.gov/Blast.cgi', new=2)

    def open_ncbi_web_page(self):
        webbrowser.open('https://www.ncbi.nlm.nih.gov/', new=2)

    def open_casper2_web_page(self):
        webbrowser.open('http://casper2.org/', new=2)

    def visit_repo_func(self):
        webbrowser.open('https://github.com/TrinhLab/CASPERapp')

    @QtCore.pyqtSlot()
    def view_results(self):
        self.hide()

        # set Results endo combo box
        self.Results.endonucleaseBox.clear()

        # set the results window endoChoice box menu
        # set the mainWindow's endoChoice first, and then loop through and set the rest of them
        self.Results.endonucleaseBox.addItem(self.endoChoice.currentText())
        for item in self.data[str(self.orgChoice.currentText())]:
            if item != self.Results.endonucleaseBox.currentText():
                self.Results.endonucleaseBox.addItem(item)

        self.Results.show()

    # this function calls the closingWindow class.
    def closeEvent(self, event):
        self.closeFunction()
        event.accept()

    # this if the function that is called when the user closes the program entirely.
    # so far I only know of 4 spots that can do this
    #       1. mainWindow
    #       2. annotationsWindow
    #       3. Results
    #       4. Multitargetting
    def closeFunction(self):
        self.myClosingWindow.get_files()
        self.myClosingWindow.show()

    def close_app(self):
        self.closeFunction()
        self.close()



# ----------------------------------------------------------------------------------------------------- #
# =========================================================================================
# CLASS NAME: StartupWindow
# Inputs: Takes information from the main application window and displays the gRNA results
# Outputs: The display of the gRNA results search
# =========================================================================================


class StartupWindow(QtWidgets.QDialog):

    def __init__(self):

        super(StartupWindow, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'startupCASPER.ui'), self)
        self.setWindowModality(2)  # sets the modality of the window to Application Modal
        #self.make_window = annotations_Window()
        #---Button Modifications---#
        self.setWindowIcon(QtGui.QIcon(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "cas9image.png")))
        pixmap = QtGui.QPixmap(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'mainart.jpg'))
        self.labelforart.setPixmap(pixmap)
        self.pushButton_2.setDefault(True)
        # Check to see the operating system you are on and change this in Global Settings:
        GlobalSettings.OPERATING_SYSTEM_ID = platform.system()
        self.info_path = os.getcwd()
        self.gdirectory = self.check_dir()
        GlobalSettings.CSPR_DB = self.gdirectory
        self.lineEdit.setText(self.gdirectory)

        self.pushButton_3.clicked.connect(self.changeDir)
        self.pushButton_2.clicked.connect(self.show_window)
        self.pushButton.clicked.connect(self.errormsgmulti)
        self.show()

    ####---FUNCTIONS TO RUN EACH BUTTON---####
    def changeDir(self):
        filed = QtWidgets.QFileDialog()
        mydir = QtWidgets.QFileDialog.getExistingDirectory(filed, "Open a Folder",
                                                       self.gdirectory, QtWidgets.QFileDialog.ShowDirsOnly)
        if(os.path.isdir(mydir) == False):
            return

        self.lineEdit.setText(mydir)
        cdir = self.lineEdit.text()
        self.gdirectory = mydir
        GlobalSettings.CSPR_DB = cdir
        #print(mydir)
        #print(cdir)

    def errormsgmulti(self):
            self.gdirectory = str(self.lineEdit.text())
            #print(self.gdirectory)
            if "Please select a directory that contains .cspr files" in self.gdirectory:
                QtWidgets.QMessageBox.question(self, "Must select directory", "You must select your directory",
                                               QtWidgets.QMessageBox.Ok)

            elif (os.path.isdir(self.gdirectory)):
                os.chdir(self.gdirectory)
                #change dir, still load main window, still load MT data, and then open main window and newGenome window
                GlobalSettings.filedir = self.gdirectory
                GlobalSettings.CASPER_FOLDER_LOCATION = self.info_path
                self.re_write_dir()
                GlobalSettings.mainWindow.launch_newGenome()
                self.close()
            else:
                QtWidgets.QMessageBox.question(self, "Not a directory", "The directory you selected does not exist.",
                                                                                    QtWidgets.QMessageBox.Ok)

    def check_dir(self):
        if GlobalSettings.OPERATING_SYSTEM_ID == "Windows":
            cspr_info = open(self.info_path+"\\CASPERinfo",'r+')
        else:
            cspr_info = open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "CASPERinfo"), 'r+')
        cspr_info = cspr_info.read()
        lines = cspr_info.split('\n')
        line = ""
        for item in lines:
            if 'DIRECTORY:' in item:
                line = item
                break
        if len(line)<11:
            return os.path.expanduser("~\Documents").replace('\\','/')
        else:
            return line[10:]

    def re_write_dir(self):
        if GlobalSettings.OPERATING_SYSTEM_ID == "Windows":
            cspr_info = open(self.info_path+"\\CASPERinfo",'r+')
        else:
            cspr_info = open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "CASPERinfo"), 'r+')
        cspr_info_text = cspr_info.read()
        cspr_info_text = cspr_info_text.split('\n')
        full_doc = ""
        for item in cspr_info_text:
            if 'DIRECTORY:' in item:
                line = item
                break
        line_final  = "DIRECTORY:"+self.gdirectory
        for item in cspr_info_text:
            if item == line:
                full_doc= full_doc+"\n"+line_final
            else:
                full_doc = full_doc+"\n" + item
        full_doc = full_doc[1:]
        cspr_info.close()
        if GlobalSettings.OPERATING_SYSTEM_ID == "Windows":
            cspr_info = open(self.info_path+"\\CASPERinfo",'w+')
        else:
            cspr_info = open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "CASPERinfo"), 'r+')
        cspr_info.write(full_doc)

        cspr_info.close()

    def show_window(self):

        self.gdirectory = str(self.lineEdit.text())
        #print(self.gdirectory)
        if "Please select a directory that contains .capr files" in self.gdirectory:
            QtWidgets.QMessageBox.question(self, "Must select directory", "You must select your directory",
                                                                                      QtWidgets.QMessageBox.Ok)
        elif(os.path.isdir(self.gdirectory)):

            os.chdir(self.gdirectory)
            found = GlobalSettings.mainWindow.getData()
            if found==False:
                QtWidgets.QMessageBox.question(self, "No Cspr files", "Please select a directory that contains cspr files.",
                                               QtWidgets.QMessageBox.Ok)
                return
            GlobalSettings.filedir = self.gdirectory
            self.re_write_dir()
            GlobalSettings.CASPER_FOLDER_LOCATION = self.info_path
            GlobalSettings.mainWindow.show()
            #Tanner - still setup data for MT
            GlobalSettings.MTWin.launch(self.gdirectory)
            self.close()
        else:
            QtWidgets.QMessageBox.question(self, "Not a directory", "The directory you selected does not exist.",
                                                                                      QtWidgets.QMessageBox.Ok)



if __name__ == '__main__':
    #enable DPI scaling
    GlobalSettings.appdir = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0]))) #used as global constant
    print(GlobalSettings.appdir)

    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

    app = Qt.QApplication(sys.argv)
    app.setOrganizationName("TrinhLab-UTK")
    app.setApplicationName("CASPER")

    startup = StartupWindow()
    GlobalSettings.mainWindow = CMainWindow(os.getcwd())
    GlobalSettings.MTWin = multitargeting.Multitargeting()
    GlobalSettings.pop_Analysis = populationAnalysis.Pop_Analysis()

    sys.exit(app.exec_())
