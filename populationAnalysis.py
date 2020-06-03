from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
import GlobalSettings
import os
from CSPRparser import CSPRparser
import Algorithms
from functools import partial
from mpl_toolkits.mplot3d import axes3d
import numpy as np
from matplotlib_venn import venn3_unweighted
import matplotlib.pyplot as plt
import itertools
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib_venn import venn3_unweighted
import matplotlib.pyplot as plt
import show_nams_ui
import show_names2_ui
import sys
from statistics import mode


class Pop_Analysis(QtWidgets.QMainWindow):
    def __init__(self):
        super(Pop_Analysis, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'populationanalysis.ui'), self)
        self.setWindowIcon(QtGui.QIcon(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "cas9image.png")))
        self.goBackButton.clicked.connect(self.go_back)
        self.analyze_button.clicked.connect(self.pre_analyze)
        self.clear_Button.clicked.connect(self.clear)
        self.ncbi_search_button.clicked.connect(self.launch_ncbi_seacher)
        self.meta_genomic_cspr_checkbox.stateChanged.connect(self.get_data)
        self.parser = CSPRparser("")
        self.Endos = dict()
        self.fna_files = dict()
        self.cspr_files = {}
        self.sq = Algorithms.SeqTranslate()
        self.ref_para_list = list()
        self.mode = 0
        self.find_locs_button.clicked.connect(self.find_locations)
        self.clear_loc_button.clicked.connect(self.clear_loc_table)
        self.directory = ""
        self.names = []
        self.names_venn = []
        self.show_names.clicked.connect(self.show_names_func)
        self.show_names2.clicked.connect(self.show_names_func2)
        self.name_form = show_nams_ui.show_names_table()
        self.name_form2 = show_names2_ui.show_names_table2()


        #organism table
        self.org_Table.setColumnCount(1)
        self.org_Table.setShowGrid(False)
        self.org_Table.setHorizontalHeaderLabels(["Organism"])
        self.org_Table.horizontalHeader().setSectionsClickable(True)
        self.org_Table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.org_Table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.org_Table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.org_Table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.org_Table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)

        #top-right table
        self.table2.setColumnCount(9)
        self.table2.setShowGrid(False)
        self.table2.setHorizontalHeaderLabels(["Seed","% Coverage","Total Repeats","Avg. Repeats/Scaffold", "Consensus Sequence", "% Consensus", "Score","PAM", "Strand"])
        self.table2.horizontalHeader().setSectionsClickable(True)
        self.table2.horizontalHeader().sectionClicked.connect(self.table2_sorting)
        self.table2.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table2.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.table2.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
        self.table2.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.table2.resizeColumnsToContents()

        # Finder table
        self.loc_finder_table.setColumnCount(5)
        self.loc_finder_table.setShowGrid(False)
        self.loc_finder_table.setHorizontalHeaderLabels(["Seed ID", "Sequence", "Organism", "Scaffold", "Location"])
        self.loc_finder_table.horizontalHeader().setSectionsClickable(True)
        self.loc_finder_table.horizontalHeader().sectionClicked.connect(self.loc_table_sorter)
        self.loc_finder_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.loc_finder_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.loc_finder_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch) #This keeps the organism column from being too small.
        self.loc_finder_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)
        self.loc_finder_table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeToContents)
        self.loc_finder_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.loc_finder_table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.loc_finder_table.setSelectionBehavior(QtWidgets.QTableView.SelectRows)
        self.loc_finder_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.loc_finder_table.resizeColumnsToContents()

        # action buttons
        self.actionMetaGenome_Parser.triggered.connect(self.launch_chrom_selector)

        self.combinerWindow = fna_and_cspr_combiner()

        self.total_org_number = 0

        self.switcher_table2 = [1, 1, 1, 1, 1, 1, 1, 1,
                                1]  # for keeping track of where we are in the sorting clicking for each column
        self.switcher_loc_table = [1, 1, 1, 1, 1]

    def launch_ncbi_seacher(self):
        GlobalSettings.mainWindow.ncbi_search_dialog.searchProgressBar.setValue(0)
        GlobalSettings.mainWindow.ncbi_search_dialog.show()

    # launches the chromesome selector window
    def launch_chrom_selector(self):
        GlobalSettings.mainWindow.cspr_selector.launch(self.cspr_files)

    def launch(self, path):
        os.chdir(path)
        self.directory = path
        self.get_data()

    # this function clears the loc_finder_table
    def clear_loc_table(self):
        self.loc_finder_table.clearContents()
        self.loc_finder_table.setRowCount(0)

    def find_locations(self):
        selectedList = self.table2.selectedItems()
        tableIndex = 0

        org_chrom_num = list()
        tempD = dict()

        # error checking
        if len(selectedList) == 0:
            QtWidgets.QMessageBox.question(self, "Error", "Please select at least 1 seed to find locations of.",
                                           QtWidgets.QMessageBox.Ok)
            self.loc_finder_table.setRowCount(0)
            return

        # get the chromosome numbers for each organism.
        counter = 0
        org_index = 0
        for item in self.ref_para_list:
            if item[0] not in tempD:
                if counter != 0:
                    org_chrom_num.append(counter)

                counter = 1
                tempD[item[0]] = org_index
                org_index += 1
            else:
                counter += 1
        if counter != 0:
            org_chrom_num.append(counter)

        # loop through and get the data
        for i in range(len(selectedList)):
            # we only want the first column's data for the popParser key
            if i % 9 == 0:
                currentSeed = selectedList[i].text()

                for item in self.parser.popData[currentSeed]:
                    self.loc_finder_table.setRowCount(tableIndex + 1)
                    tempSeq = item[3]
                    tempOrg = item[0]
                    tempChrom = int(item[1])
                    tempLoc = item[2]

                    # calculate the true chromosome number
                    tempIndex = tempD[tempOrg]
                    while tempIndex > 0:
                        tempChrom = tempChrom - org_chrom_num[int(tempIndex - 1)]
                        tempIndex -= 1

                    tabSeed = QtWidgets.QTableWidgetItem()
                    tabSeed.setData(QtCore.Qt.EditRole, currentSeed)
                    tabSeq = QtWidgets.QTableWidgetItem()
                    tabSeq.setData(QtCore.Qt.EditRole, tempSeq)
                    tabOrg = QtWidgets.QTableWidgetItem()
                    tabOrg.setData(QtCore.Qt.EditRole, tempOrg)
                    tabChrom = QtWidgets.QTableWidgetItem()
                    tabChrom.setData(QtCore.Qt.EditRole, int(tempChrom))
                    tabLoc = QtWidgets.QTableWidgetItem()
                    tabLoc.setData(QtCore.Qt.EditRole, int(tempLoc))

                    self.loc_finder_table.setItem(tableIndex, 0, tabSeed)
                    self.loc_finder_table.setItem(tableIndex, 1, tabSeq)
                    self.loc_finder_table.setItem(tableIndex, 2, tabOrg)
                    self.loc_finder_table.setItem(tableIndex, 3, tabChrom)
                    self.loc_finder_table.setItem(tableIndex, 4, tabLoc)
                    tableIndex += 1

        self.loc_finder_table.resizeColumnsToContents()

        # this function builds the Select Organisms table

    def get_data(self):
        if self.directory == '':
            return

        onlyfiles = [f for f in os.listdir(self.directory) if os.path.isfile(os.path.join(self.directory, f))]
        self.fna_files.clear()

        # show/hide the stuff that isn't needed
        if self.meta_genomic_cspr_checkbox.isChecked():
            self.endoBox.hide()
            self.ncbi_search_button.hide()
            self.label_3.hide()
            self.label_2.setText('Select a Metagenomic CSPR File')
        elif not self.meta_genomic_cspr_checkbox.isChecked():
            self.endoBox.show()
            self.ncbi_search_button.show()
            self.label_3.show()
            self.label_2.setText('Select organism(s) and endonuclease:')

        # if the user wants the FNA/Fast files to be sown
        if not self.meta_genomic_cspr_checkbox.isChecked():
            index = 0
            for file in onlyfiles:
                if file.find('.fna') != -1 or file.find('.fasta') != -1:
                    # find the organism name
                    f = open(file, 'r')
                    hold = f.readline()
                    f.close()
                    # get the organism name
                    spaceIndex = hold.find(' ') + 1
                    commaIndex = hold.find(',')
                    buf = hold[spaceIndex:commaIndex]

                    # store the name in the dict of fna_files, that keys the name with the file path
                    self.fna_files[buf] = file

                    # store the data in the table
                    tabWidget = QtWidgets.QTableWidgetItem(buf)
                    self.org_Table.setRowCount(index + 1)
                    self.org_Table.setItem(index, 0, tabWidget)
                    index += 1
        # if the user wants the metagenomic cspr files to be shown
        else:
            index = 0
            for file in onlyfiles:
                if file.find('.cspr') != -1:
                    f = open(file, 'r')
                    hold = f.readline()
                    f.close()
                    # only show the files that are metagenomic
                    if '(meta)' in hold:
                        colonIndex = hold.find(':') + 1
                        commaIndex = hold.find(',')
                        orgName = hold[colonIndex:commaIndex]

                        self.fna_files[orgName] = file

                        tabWidget = QtWidgets.QTableWidgetItem(orgName)
                        self.org_Table.setRowCount(index + 1)
                        self.org_Table.setItem(index, 0, tabWidget)
                        index += 1
            if index == 0:
                self.org_Table.clearContents()
                self.org_Table.setRowCount(0)

#        self.org_Table.resizeColumnsToContents() ##Commenting this out allows the header to remain full sized

        self.fillEndo()
        # self.changeEndos()

    # this function opens CASPERinfo and builds the dropdown menu of selectable endonucleases
    def fillEndo(self):
        if GlobalSettings.OPERATING_SYSTEM_ID == "Windows":
            f = open(GlobalSettings.appdir + "\\CASPERinfo")
        else:
            f = open(GlobalSettings.appdir + '/CASPERinfo')
        while True:
            line = f.readline()
            if line.startswith('ENDONUCLEASES'):
                while True:
                    line = f.readline()
                    if (line[0] == "-"):
                        break
                    line_tokened = line.split(";")
                    endo = line_tokened[0]
                    # Checking to see if there is more than one pam sequence in the list
                    if line_tokened[1].find(",") != -1:
                        p_pam = line_tokened[1].split(",")[0]
                    else:
                        p_pam = line_tokened[1]
                    default_seed_length = line_tokened[2]
                    default_tot_length = line_tokened[3]
                    self.Endos[endo + " PAM: " + p_pam] = (endo, p_pam, default_seed_length, default_tot_length)

                break
        f.close()
        self.endoBox.addItem("None Selected")
        self.endoBox.addItems(self.Endos.keys())
        # self.endoBox.currentIndexChanged.connect(self.changeEndos)

    # this function displays all of the organisms of which the user has that endo in their DB
    def changeEndos(self):
        endo_box = str(self.endoBox.currentText())
        endo_box = endo_box[:endo_box.find(" ")]
        self.org_Table.setRowCount(0)
        index = 0
        for keys in self.cspr_files.keys():
            filename = str(self.cspr_files[keys])
            endo = filename[filename.find("_") + 1:]
            endo = endo.replace(".cspr", "")
            if (endo == endo_box or endo_box == "None"):
                self.org_Table.setRowCount(index + 1)
                name = QtWidgets.QTableWidgetItem(str(keys))
                self.org_Table.setItem(index, 0, name)
                index += 1
        self.org_Table.resizeColumnsToContents()

    # this function calls the popParser function and fills all the tables
    def pre_analyze(self):
        # if the user is wanting to go with 1 meta genomic cspr file
        if self.meta_genomic_cspr_checkbox.isChecked():
            selectedList = self.org_Table.selectedItems()

            # error check
            if len(selectedList) == 0 or len(selectedList) > 1:
                QtWidgets.QMessageBox.question(self, "Error", "Please select no more than 1 CSPR file for analysis.",
                                               QtWidgets.QMessageBox.Ok)
                return

            # get the cspr_file name, the endochoice, and call the popParser
            orgName = selectedList[0].text()
            cspr_file_name = self.fna_files[orgName]
            # split the file name by '_', then take that second index, split by '.', and then take the first index. Thus giving the Endo Choice
            endoChoice = cspr_file_name.split('_')[1].split('.')[0]
            # call the parser and the call fill_data
            cspr_file_name = GlobalSettings.CSPR_DB + '/' + cspr_file_name
            self.total_org_number, self.ref_para_list = self.parser.popParser(cspr_file_name, endoChoice)
            self.fill_data(endoChoice)
        # if the user is wanting to go with creating a new meta genomic cspr file
        else:
            selectedList = self.org_Table.selectedItems()

            # if the table is showing only fna/fasta files
            if not self.meta_genomic_cspr_checkbox.isChecked():
                # rules for selecting FNA/Fasta files
                # check to make sure that the user selected at least 2 organisms, and 1 endonuclease
                if len(selectedList) < 1 or self.endoBox.currentText() == 'None Selected':
                    QtWidgets.QMessageBox.question(self, "Nothing Seleted", "No items selected. Please select at least 1 organism and only 1 endonuclease.",
                                                QtWidgets.QMessageBox.Ok)
                    return
                if len(selectedList) == 1:
                    error = QtWidgets.QMessageBox.question(self, "Only 1 Organism Selected",
                                                            "Population Analysis works with multiple organisms, or a metagenome. If the file selected is not a metagenome, the program may not function correctly. Do you wish to continue? \n\n"
                                                            "Do you wish to continue?",
                                                            QtWidgets.QMessageBox.Yes |
                                                            QtWidgets.QMessageBox.No,
                                                            QtWidgets.QMessageBox.No)
                    if (error == QtWidgets.QMessageBox.No):
                        return

                submitList = list()
                for item in selectedList:
                    submitList.append(self.fna_files[item.text()])
                self.combinerWindow.launch(submitList)
            # rules for selecting cspr files
            elif self.meta_genomic_cspr_checkbox.isChecked():
                if len(selectedList) == 0:
                    QtWidgets.QMessageBox.question(self, "Nothing Seleted", "No items selected. Please select one metagenome for Population Analysis.",
                                                QtWidgets.QMessageBox.Ok)
                    return
                elif len(selectedList) > 1:
                    QtWidgets.QMessageBox.question(self, "Too Many Selected", "Only 1 metagenomic CSPR file is allowed to be selected.",
                                                QtWidgets.QMessageBox.Ok)
                    return

    # this function calculates the percentConserved for the table
    # it runs through and finds out how many different organisms each seed is repeated in
    # if it's equal to the total_org_number, it then returns 1, otherwise it returns a double
    def findPercentConserved(self, seed):
        tempSet = set()
        for item in self.parser.popData[seed]:
            tempSet.add(item[0])

        if self.total_org_number == len(tempSet):
            return 1
        else:
            return len(tempSet) / self.total_org_number

    # this function calculates the average repeats per chromosome for a seed
    # runs through the sequences in a seed and calculates it
    # returns the average
    def findAvgRepeats(self, seed):
        firstChrom = 0
        secondChrom = 0
        divideBy = 1
        tempSum = 0
        for item in self.parser.popData[seed]:
            firstChrom = item[1]

            if firstChrom != secondChrom and secondChrom != 0:
                divideBy += 1

            tempSum += 1
            secondChrom = item[1]

        return tempSum / divideBy

    # this function fills the top-right table
    def fill_data(self, endoChoice):
        self.table2.setRowCount(0)
        index = 0
        for seeds in self.parser.popData:
            self.table2.setRowCount(index + 1)

            seed = QtWidgets.QTableWidgetItem()
            total_repeats = QtWidgets.QTableWidgetItem()
            total_repeats.setData(QtCore.Qt.EditRole, len(self.parser.popData[seeds]))
            seed.setData(QtCore.Qt.EditRole, str(seeds))

            self.table2.setItem(index, 0, seed)
            self.table2.setItem(index, 2, total_repeats)
            tempPercentConserved = self.findPercentConserved(seeds) * 100
            tempPercentConserved = float("%.2f" % tempPercentConserved)
            percentTab = QtWidgets.QTableWidgetItem(str(tempPercentConserved) + '%')
            self.table2.setItem(index, 1, percentTab)

            # get the avg repeats per chromosome
            tempAvgRepeatsPerChrom = self.findAvgRepeats(seeds)
            rounded = float("%.2f" % tempAvgRepeatsPerChrom)
            avgTab = QtWidgets.QTableWidgetItem()
            avgTab.setData(QtCore.Qt.EditRole, rounded)
            self.table2.setItem(index, 3, avgTab)

            # loop through the tuples for each seed
            sequences = []
            for tuples in self.parser.popData[seeds]:
                sequences.append(tuples[3])

            # set consensus seq
            con_seq_temp = str(max(set(sequences), key=sequences.count))
            conIndex = sequences.index(con_seq_temp)
            consensus_seq = QtWidgets.QTableWidgetItem()
            consensus_seq.setData(QtCore.Qt.EditRole, con_seq_temp)
            self.table2.setItem(index, 4, consensus_seq)

            # get the data for the rest
            tabScore = QtWidgets.QTableWidgetItem()
            tabScore.setData(QtCore.Qt.EditRole, int(self.parser.popData[seeds][conIndex][5]))
            tabPAM = QtWidgets.QTableWidgetItem(self.parser.popData[seeds][conIndex][4])
            tabStrand = QtWidgets.QTableWidgetItem(self.parser.popData[seeds][conIndex][6])

            # set all that data
            self.table2.setItem(index, 6, tabScore)
            self.table2.setItem(index, 7, tabPAM)
            self.table2.setItem(index, 8, tabStrand)

            # set consensus percentage
            consensus_percentage = sequences.count(con_seq_temp) / len(sequences) * 100
            consensus_percentage = round(consensus_percentage, 2)
            consensus_percentage = format(consensus_percentage, '.2f')
            consensus_perc = QtWidgets.QTableWidgetItem()
            consensus_perc.setData(QtCore.Qt.EditRole, consensus_percentage)
            self.table2.setItem(index, 5, consensus_perc)

            index += 1
        self.table2.resizeColumnsToContents()
        self.plot_repeats_vs_seeds(endoChoice)
        self.plot_3D_graph(endoChoice)
        self.plot_venn()

    def clear(self):
        self.table2.setRowCount(0)

    # this function graphs the repeats_vs_seeds graph
    def plot_repeats_vs_seeds(self, endoChoice):
        data = {}
        for seed in self.parser.popData:
            number = 0
            for repeat in self.parser.popData[seed]:
                number += 1
            if number in data:
                data[number] += 1
            else:
                data[number] = 1

        max = 0
        y1 = []
        x1 = []
        plots = []
        time = 0
        index = 0

        for number in self.order(data):
            time += 1

            if int(data[number]) > max:
                max = int(data[number])
                self.mode = number

            hold = 0
            while hold < data[number]:
                if index == int(round(len(self.parser.popData) / 2)):
                    self.median = number
                x1.append(index)
                y1.append(number)

                index = index + 1
                hold += 1

        # now plot the stuff
        self.pop_analysis_repeats_graph.canvas.axes.clear()
        # set everything up
        self.pop_analysis_repeats_graph.canvas.axes.plot(x1, y1)
        self.pop_analysis_repeats_graph.canvas.axes.set_xlabel('Seed ID Number')
        self.pop_analysis_repeats_graph.canvas.axes.set_ylabel('Number of Repeats')
        self.pop_analysis_repeats_graph.canvas.axes.set_title('Number of Repeats per Seed ID Number')
        # now draw
        self.pop_analysis_repeats_graph.canvas.draw()

    # this function is for the 3D bar graph
    def plot_3D_graph(self, endoChoice):
        if (self.total_org_number > 2):
            rows, cols = (self.total_org_number, self.total_org_number)

            arr = [[0 for i in range(cols)] for j in range(rows)]

            x3 = []
            y3 = []
            z3 = np.zeros(int((self.total_org_number * (self.total_org_number - 1)) / 2))
            dz = []
            self.names = []
            dx = np.ones(int((self.total_org_number * (self.total_org_number - 1)) / 2))
            dy = np.ones(int((self.total_org_number * (self.total_org_number - 1)) / 2))

            for keys in self.parser.popData:
                for items in self.parser.popData[keys]:
                    if items[0] not in self.names:
                        self.names.append(items[0])

            for keys in self.parser.popData:
                temp_names = []
                for items in self.parser.popData[keys]:
                    if items[0] not in temp_names:
                        temp_names.append(items[0])

                if len(temp_names) >= 2:
                    for i in range(len(temp_names) - 1):
                        j = i + 1
                        while j != len(temp_names):
                            arr[self.names.index(temp_names[i])][self.names.index(temp_names[j])] += 1
                            arr[self.names.index(temp_names[j])][self.names.index(temp_names[i])] += 1
                            j += 1

            for j in range(cols):
                i = len(self.names) - 1
                while i != j:
                    x3.append(i)
                    y3.append(j)
                    dz.append(arr[i][j])
                    i -= 1

            self.pop_analysis_3dgraph.canvas.axes.clear()
            self.pop_analysis_3dgraph.canvas.axes.bar3d(x3, y3, z3, dx, dy, dz)

            new_names = []

            for n in range(len(self.names)):
                new_names.append(n)

            self.pop_analysis_3dgraph.canvas.axes.set_xlabel('x')
            self.pop_analysis_3dgraph.canvas.axes.set_ylabel('y')
            self.pop_analysis_3dgraph.canvas.axes.set_zlabel('z')

            self.pop_analysis_3dgraph.canvas.axes.set_xticks(np.arange(1, self.total_org_number + 1, 1))
            self.pop_analysis_3dgraph.canvas.axes.set_yticks(np.arange(0, self.total_org_number, 1))

            self.pop_analysis_3dgraph.canvas.axes.tick_params(labelsize=8)
            self.pop_analysis_3dgraph.canvas.axes.set_xticklabels(new_names, rotation=45)
            self.pop_analysis_3dgraph.canvas.axes.set_yticklabels(new_names, rotation=-45)

            self.pop_analysis_3dgraph.canvas.draw()

    def plot_venn(self):
        self.pop_analysis_venn_diagram.canvas.figure.clf()
        # self.pop_analysis_venn_diagram.canvas = FigureCanvas(plt.figure(figsize=(7.5,7.5)))

        rows, cols = (self.total_org_number, self.total_org_number)
        arr = [[0 for i in range(cols)] for j in range(rows)]
        all_3 = 0
        singles = [0 for i in range(cols)]
        counter = 0
        self.names_venn = []

        for keys in self.parser.popData:
            for items in self.parser.popData[keys]:
                if items[0] not in self.names_venn:
                    self.names_venn.append(items[0])

        if len(self.names_venn) >= 3:
            for keys in self.parser.popData:
                temp_names = []
                for items in self.parser.popData[keys]:
                    if items[0] not in temp_names:
                        temp_names.append(items[0])

                if len(temp_names) >= 2:
                    for i in range(len(temp_names) - 1):
                        j = i + 1
                        while j != len(temp_names):
                            arr[self.names_venn.index(temp_names[i])][self.names_venn.index(temp_names[j])] += 1
                            arr[self.names_venn.index(temp_names[j])][self.names_venn.index(temp_names[i])] += 1
                            j += 1
                else:
                    if temp_names[0] == 'Corynebacterium tuberculostearicum SK141 contig00004':
                        counter += 1
                    singles[self.names_venn.index(temp_names[0])] += 1

            # all 3 orgs
            for keys in self.parser.popData:
                temp_names = []
                for items in self.parser.popData[keys]:
                    if items[0] not in temp_names:
                        temp_names.append(items[0])
                if all(x in temp_names for x in [self.names_venn[0], self.names_venn[1], self.names_venn[2]]):
                    all_3 += 1

            venn3_unweighted(subsets=(singles[0], singles[1], arr[0][1], singles[2], arr[0][2],
                                      arr[1][1], all_3), set_labels=('0', '1', '2'))

            self.pop_analysis_venn_diagram.canvas.draw()
        else:
            self.pop_analysis_venn_diagram.canvas.figure.clf()
            self.pop_analysis_venn_diagram.canvas.draw()

            # del self.pop_analysis_venn_diagram.canvas.axes
            # self.pop_analysis_venn_diagram.canvas.draw()

            # self.pop_analysis_venn_diagram.canvas.axes = venn3_unweighted(
            # subsets=('null', 'null', 'null', 'null', 'null', 'null', 'null'), set_labels=('Org. 1', 'Org. 2', 'Org. 3'))

    def show_names_func(self):
        # print(self.names)
        self.name_form.fill_table(self.names)
        self.name_form.show()

    def show_names_func2(self):
        # print(self.names)
        if len(self.names_venn) >= 3:
            self.name_form2.fill_table(self.names_venn[0:3])
            self.name_form2.show()
        else:
            self.name_form2.name_table2.setRowCount(0)

    def order(self, data_par):
        data = dict(data_par)
        data2 = []
        while len(data) > 0:
            max = 0
            for item in data:
                if item > max:
                    max = item
            data2.append(max)
            if len(data2) == 1:
                self.max_repeats = max
            del data[max]
        return data2

    def go_back(self):
        GlobalSettings.mainWindow.getData()
        GlobalSettings.mainWindow.show()
        self.hide()

    # this function calls the close window class. Allows the user to choose what files they want to keep/delete
    def closeEvent(self, event):
        GlobalSettings.mainWindow.closeFunction()
        event.accept()

    def table_sorting(self, logicalIndex):
        self.switcher[logicalIndex] *= -1
        if self.switcher[logicalIndex] == -1:
            self.table2.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
        else:
            self.table2.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)

    # sorting to table2: IE the table in top-right
    def table2_sorting(self, logicalIndex):
        self.switcher_table2[logicalIndex] *= -1
        if self.switcher_table2[logicalIndex] == -1:
            self.table2.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
        else:
            self.table2.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)

    # sorting for location table: IE table in bottom right
    def loc_table_sorter(self, logicalIndex):
        self.switcher_loc_table[logicalIndex] *= -1
        if (self.switcher_loc_table[logicalIndex] == -1):
            self.loc_finder_table.sortItems(logicalIndex, QtCore.Qt.DescendingOrder)
        else:
            self.loc_finder_table.sortItems(logicalIndex, QtCore.Qt.AscendingOrder)


###############################################################################
# class name: fna_and_cspr_combiner
# this opens a window for the user to select the Organism Name, FNA File Name, and the organism code.
# It combines the FNA files, and then runs the sequencer. This way it creates 1 cspr file.
# after the new cspr file is created, it will run the population analysis
###############################################################################
class fna_and_cspr_combiner(QtWidgets.QDialog):
    def __init__(self):
        # Qt init stuff
        super(fna_and_cspr_combiner, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "pop_analysis_fna_combiner.ui"), self)
        self.sequencer_prog_bar.setValue(0)

        # button connections
        self.cancel_button.clicked.connect(self.cancel_function)
        self.start_button.clicked.connect(self.run_analysis)

        # variables below
        self.ref_para_list = list()
        self.fna_file_names = list()
        self.sq = Algorithms.SeqTranslate()
        self.proc_running = False
        self.combined_fna_file = ''
        self.cancelled = False

    # for when the user clicks the 'x' button
    def closeEvent(self, event):
        closeWindow = self.cancel_function()

        # if the user is doing OT and does not decide to cancel it ignore the event
        if closeWindow == -2:
            event.ignore()
        else:
            event.accept()

    # this is the function that will eventually run the analysis. For now, it is just running the combine/build new cspr function
    def run_analysis(self):
        # make sure the process isn't already running
        if self.proc_running:
            return

        # make sure the user inputs whats needed
        if self.orgName_line_edit.text() == '' or self.org_code_line_edit.text() == '' or self.orgNum_lineEdit.text() == '':
            QtWidgets.QMessageBox.question(self, "Missing Information",
                                           "Please input an organism name, an organism code, and the number of organisms you are analyzing.",
                                           QtWidgets.QMessageBox.Ok)
            return

        # make sure the user inputs an integer for the number of organisms
        if not self.orgNum_lineEdit.text().isdigit():
            QtWidgets.QMessageBox.question(self, "Error",
                                           "Organism number must be integers only!",
                                           QtWidgets.QMessageBox.Ok)
            return

        self.combine_fna_files()
        self.build_new_cspr_file()

    # cancel function, just clears everything and closes the window
    def cancel_function(self):
        # check to see if the sequencer is running. If so ask the user if they wish to close out
        if self.proc_running:
            error = QtWidgets.QMessageBox.question(self, "Sequencer Is Running",
                                                   "Sequencer is running. Closing this window will cancel that process and return to the Population Analysis window. \n\n"
                                                   "Do you wish to continue?",
                                                   QtWidgets.QMessageBox.Yes |
                                                   QtWidgets.QMessageBox.No,
                                                   QtWidgets.QMessageBox.No)
            if (error == QtWidgets.QMessageBox.No):
                return -2
            else:
                self.cancelled = True
                self.process.kill()

        # close out and leave
        self.orgName_line_edit.setText('')
        self.org_code_line_edit.setText('')
        self.sequencer_prog_bar.setValue(0)
        self.hide()

    # open the window and get the fna files that the user wishes to use
    def launch(self, filenames):
        self.fna_file_names = filenames
        self.process = QtCore.QProcess()
        self.show()

    # this function takes a list of file paths, which are FNA or Fasta files, and combines them into one file
    # this function also builds a parallel/reference list of the chromosomes
    def combine_fna_files(self):
        self.ref_para_list.clear()

        # open the output file (currently just a test file)
        self.combined_fna_file = GlobalSettings.CSPR_DB + '/' + 'temp.fna'
        out_stream = open(self.combined_fna_file, 'w')

        # for each file in the list
        for file in self.fna_file_names:
            file_stream = open(file, 'r')

            buf = file_stream.readline()

            spaceIndex = buf.find(' ') + 1
            commaIndex = buf.find(',')
            orgName = buf[spaceIndex:commaIndex]

            # read the whole file
            while buf != "":

                # if buf is empty, break out
                if buf == '' or buf == '\n':
                    break

                # if buf starts with a '>', it's a chromosome so make sure to store it in the ref_para list
                if buf.startswith('>'):
                    self.ref_para_list.append((orgName, buf))
                # right the buf to the new file, and re-read it
                out_stream.write(buf)
                buf = file_stream.readline()
            file_stream.close()
        out_stream.close()

    # this function builds a new cspr file from the combined FNA file
    # very similar to New Genome
    def build_new_cspr_file(self):
        self.num_chromo_next = False
        self.num_chromo = 0

        # this function reads output. Just used to populate the progress bar
        def output_stdout(p):
            line = str(p.readAll())
            line = line[2:]
            line = line[:len(line) - 1]
            line = line.strip('\\n')
            line = line.strip('\n')
            for lines in filter(None, line.split(r'\n')):
                lines = lines.strip('\n')
                lines = lines.strip('\\n')
                print(lines)
                if (lines == 'Finished reading in the genome file.'):
                    self.num_chromo_next = True
                elif (self.num_chromo_next == True):
                    lines = lines.strip('\\n')
                    self.num_chromo_next = False
                    self.num_chromo = int(lines)
                elif (lines.find('Chromosome') != -1 and lines.find('complete.') != -1):
                    temp = lines
                    temp = temp.replace('Chromosome ', '')
                    temp = temp.replace(' complete.', '')
                    temp = temp.strip('\n')
                    if (int(temp) == self.num_chromo):
                        self.sequencer_prog_bar.setValue(99)
                    else:
                        self.sequencer_prog_bar.setValue(int(temp) / self.num_chromo * 100)
                elif (lines == 'Finished Creating File.'):
                    self.sequencer_prog_bar.setValue(100)

        # this function will end up doing stuff when the process is finished.
        def finished():
            self.proc_running = False

            if self.cancelled == False:
                # get the file name
                cspr_file_name = GlobalSettings.CSPR_DB + '/' + self.org_code_line_edit.text() + '_' + \
                                 GlobalSettings.pop_Analysis.endoBox.currentText().split(' ')[0] + '.cspr'
                self.process.kill()
                endoChoice = GlobalSettings.pop_Analysis.endoBox.currentText().split(' ')[0]
                GlobalSettings.pop_Analysis.total_org_number, GlobalSettings.pop_Analysis.ref_para_list = GlobalSettings.pop_Analysis.parser.popParser(
                    cspr_file_name, endoChoice)
                GlobalSettings.pop_Analysis.fill_data(endoChoice)
            os.remove(self.combined_fna_file)
            self.cancelled = False
            self.cancel_function()

        # --------------getting the arugments---------------------------------
        # get the file path to the combined fna file
        fna_file_path = self.combined_fna_file
        path_to_fna = fna_file_path
        # get the endo_choice, hard coded for now, will eventually from with the user
        endo_choice = GlobalSettings.pop_Analysis.endoBox.currentText().split(' ')[0]
        # get the pam itself, taken from the endo_choice. Also, if there's multiple, always take the first one
        pam = self.sq.endo_info[endo_choice][0].split(',')[0]
        # get the code. currently hard coded to test_code until I know what to put it to later one
        code = self.org_code_line_edit.text()
        # check to see if the seq_finder should do 5' prime or no
        if int(self.sq.endo_info[endo_choice][3]) == 3:
            pamdir = False
        else:
            pamdir = True
        # get the output folder location
        output_location = GlobalSettings.CSPR_DB
        # get the path to CASPERinfo
        path_to_info = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'CASPERinfo')
        # make org name something that will make more sense, this is just for testing right now
        orgName = self.orgName_line_edit.text() + '  , (meta), ' + self.orgNum_lineEdit.text()
        # get the seed and RNA length, based on the endo choice
        gRNA_length = self.sq.endo_info[endo_choice][2]
        seed_length = self.sq.endo_info[endo_choice][1]

        secondCode = ''
        # get the notes here
        for i in range(len(self.ref_para_list)):
            secondCode = secondCode + self.ref_para_list[i][0] + ',' + self.ref_para_list[i][1] + '|'
            secondCode = secondCode.replace('\n', '')
        #------------------done getting the arguments------------------------------

        # get the program path
        program = '"' + os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'Casper_Seq_Finder') + '" '

        # compile all of the arguments into one line
        args = '"' + endo_choice + '" '
        args = args + '"' + pam + '" '
        args = args + '"' + code + '" '
        args = args + str(pamdir) + ' '
        args = args + '"' + output_location + '/" '
        args = args + '"' + path_to_info + '" '
        args = args + '"' + path_to_fna + '" '
        args = args + '"' + orgName + '" '
        args = args + gRNA_length + ' '
        args = args + seed_length + ' '
        args = args + '"' + secondCode + '"'
        print(args)
        # combine the program and arguments into 1
        program = program + args

        self.process.readyReadStandardOutput.connect(partial(output_stdout, self.process))
        self.proc_running = True
        self.process.start(program)
        self.process.finished.connect(finished)
