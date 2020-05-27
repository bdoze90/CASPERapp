import GlobalSettings
from PyQt5 import QtWidgets, Qt, uic
from CSPRparser import CSPRparser

####################################################################################################
# Class: cspr_chromosome_selection
# This class takes a filename as the constructor
# It opens a window that shows the chromosomes that are in the CSPR file provided
# the user can then choose which chromosome they want to put into a new CSPR file
# it creates a new CSPR file that is compatible with the parser from the chromosomes selected
####################################################################################################

class cspr_chromosome_selection(QtWidgets.QDialog):
    # the init function takes the cspr file name
    # it sets up the window
    # it searches the CSPR file for the chromosomes
    # it loads those chromosomes into the table
    def __init__(self):
        # qt stuff
        super(cspr_chromosome_selection, self).__init__()
        uic.loadUi("cspr_chromosome_selection.ui", self)
        self.setWindowTitle("Choose which chromosomes to pull data from!")
        self.setWindowIcon(Qt.QIcon("cas9image.png"))

        # button connections
        self.cancel_button.clicked.connect(self.cancel_function)
        self.submit_button.clicked.connect(self.submit_function)
        self.select_button.clicked.connect(self.load_chrom_names)

        # chrom_table stuff
        self.chromosome_table.setColumnCount(1)
        self.chromosome_table.setShowGrid(True)
        self.chromosome_table.setHorizontalHeaderLabels("chromosome;".split(";"))
        self.chromosome_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.chromosome_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.chromosome_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

        # cspr_files_available_table
        self.cspr_files_available_table.setColumnCount(1)
        self.cspr_files_available_table.setShowGrid(True)
        self.cspr_files_available_table.setHorizontalHeaderLabels("Files;".split(";"))
        self.cspr_files_available_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.cspr_files_available_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.cspr_files_available_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

        # variables
        self.cspr_file  = ''
        self.gene = ''
        self.misc = ''
        self.avail_cspr = dict()
        self.myParser = CSPRparser(self.cspr_file)
        self.orgName = ''


    # launch function
    # parameter is a dictionary. Key is a organism name, value is the file name
    # parameter is passed from the popAnalysis class, it is the cspr_files dict
    def launch(self, cspr_file_data):
        self.avail_cspr = cspr_file_data
        self.cspr_files_available_table.setRowCount(len(self.avail_cspr))

        # populate the cspr_files_available_table
        loopCount = 0
        for item in self.avail_cspr:
            tabWidget = QtWidgets.QTableWidgetItem(item)
            self.cspr_files_available_table.setItem(loopCount, 0, tabWidget)
            loopCount += 1

        self.cspr_files_available_table.resizeColumnsToContents()
        self.show()

    # load_chrom_names makes sure that the user only selects 1 file
    # then it loads all of the chroms in that file, and populates the chromosome_table
    def load_chrom_names(self):
        selected = self.cspr_files_available_table.selectedItems()

        # error checking
        if len(selected) == 0:
            QtWidgets.QMessageBox.question(self, "Nothing Selected",
                                           "No organism selected. Please selected 1 File",
                                           QtWidgets.QMessageBox.Ok)
            return
        elif len(selected) > 1:
            QtWidgets.QMessageBox.question(self, "Too Many Selected",
                                           "Please only choose 1 File",
                                           QtWidgets.QMessageBox.Ok)
            return

        # get the data from the cspr file
        self.chromosome_table.clearContents()
        self.orgName = selected[0].text()
        self.file_name_line_edit.setText(self.avail_cspr[self.orgName])
        self.cspr_file = self.avail_cspr[self.orgName]
        self.myParser.fileName = self.avail_cspr[self.orgName]
        self.gene, self.misc = self.myParser.get_chromosome_names()

        # loop through and set the table
        loopCount = 0
        self.chromosome_table.setRowCount(len(self.myParser.chromosomesSelectedList))
        for item in self.myParser.chromosomesSelectedList:
            tabWidget = QtWidgets.QTableWidgetItem(item)
            self.chromosome_table.setItem(loopCount, 0, tabWidget)
            loopCount += 1
        self.chromosome_table.resizeColumnsToContents()

    # this function builds a new CSPR file from the selected
    # gets the Genome/Misc line from the main file, also gets the correct Karystats numbers
    # gets the chromosome data from the chromosomes selected
    # gets the repeats data for the chromosomes selected
    # it also updates the org_table and corresponding dict in the pop_analysis class
    def submit_function(self):
        selectedList = self.chromosome_table.selectedItems()
        # make sure at least one chrom is selected
        if len(selectedList) == 0:
            QtWidgets.QMessageBox.question(self, "Nothing Selected",
                                           "No items selected, please select at least 1 chromosome to pull out.",
                                           QtWidgets.QMessageBox.Ok)
            return

        # now make sure they have something in the subtext field. if not, throw a warning
        if self.output_file_name.text() == '':
            QtWidgets.QMessageBox.question(self, "Error! No file subtext given.",
                                           "Please give a subtext for the new CSPR file to differentiate the files.",
                                           QtWidgets.QMessageBox.Ok)
            return
        # get the output file name info
        dotIndex = self.cspr_file.find('.')
        outputFile = self.cspr_file[:dotIndex]
        outputFile = outputFile + self.output_file_name.text() + '.cspr'

        output_list = list()
        keep_list = list()
        for item in selectedList:
            output_list.append(item.text())
            keep_list.append(item.row())

        # remove the ones the user does not want from the karystats list
        for i in range(len(self.myParser.karystatsList), 0, -1):
            if i - 1 not in keep_list:
                self.myParser.karystatsList.pop(i - 1)


        # write the first three lines
        outputStream = open(outputFile, 'w')
        self.gene = self.gene.replace('\n', '')
        outputStream.write(self.gene + self.output_file_name.text())
        outputStream.write('\n')
        outputStream.write('KARYSTATS: ')
        # write all of the karystats data
        for item in self.myParser.karystatsList:
            outputStream.write(str(item))
            outputStream.write(',')
        outputStream.write('\n')
        outputStream.write(self.misc)

        # now go through and write the chromosomes
        for item in output_list:
            outputStream.write(item)
            csprData = open(self.cspr_file, 'r')
            buffer = csprData.readline()
            while True:
                if buffer == item:
                    buffer = csprData.readline()
                    while '>' not in buffer:
                        if buffer == 'REPEATS\n':
                            break
                        outputStream.write(buffer)
                        buffer = csprData.readline()
                    break
                buffer = csprData.readline()
            csprData.close()
        outputStream.write("REPEATS\n")
# now to go through and write all of the repeats
        csprFile = open(self.cspr_file, 'r')
        csprData = csprFile.read().split('\n')

        # skip to the repeats part of the file
        index = csprData.index('REPEATS')
        index += 1

        storeRepeats = []
       # loop through the rest of the file
        while(index + 1 < len(csprData)):
            storeRepeats.clear()
            seed = csprData[index]
            repeats = csprData[index + 1].split('\t')

            # loop through and store the repeats i want
            for repeat in repeats:
                checkData = repeat.split(',')
                if len(checkData) > 1:
                    if int(checkData[0]) - 1 in keep_list:
                        storeRepeats.append(repeat)

            # if i found repeats I want, print it all to the file
            if len(storeRepeats) > 0:
                outputStream.write(seed)
                outputStream.write('\n')
                for item in storeRepeats:
                    outputStream.write(item)
                    outputStream.write('\t')
                outputStream.write('\n')
            index += 2

        # write the end of file, and close the file
        outputStream.write('END_OF_FILE')
        outputStream.close()

        # set the table and dict in the pop Analysis class
        GlobalSettings.pop_Analysis.cspr_files[self.orgName + self.output_file_name.text()] = outputFile
        GlobalSettings.pop_Analysis.org_Table.setRowCount(GlobalSettings.pop_Analysis.org_Table.rowCount() + 1)
        tabWidget = QtWidgets.QTableWidgetItem(self.orgName + self.output_file_name.text())
        GlobalSettings.pop_Analysis.org_Table.setItem(GlobalSettings.pop_Analysis.org_Table.rowCount() - 1, 0, tabWidget)
        # call the cancel function to close out the window and clear everything
        self.cancel_function()




    # cancel function
    # hides the window
    # clears table contents
    # clears the text variables as well
    def cancel_function(self):
        self.chromosome_table.clearContents()
        self.chromosome_table.setRowCount(0)
        self.cspr_files_available_table.clearContents()
        self.cspr_files_available_table.setRowCount(0)
        self.avail_cspr = dict()
        self.file_name_line_edit.setText(" ")
        self.cspr_file = " "
        self.output_file_name.setText("_Abridged")
        self.myParser.fileName = ' '
        self.gene = ' '
        self.misc = ' '
        self.orgName = ' '
        self.hide()
