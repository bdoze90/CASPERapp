import sys, os
from PyQt5 import QtWidgets, uic, QtGui, QtCore, Qt
from bioservices import KEGG
import GlobalSettings
from functools import partial
from Algorithms import SeqTranslate
from NCBI_Search_Window import NCBI_Search_File

def iter_except(function, exception):
    """Works like builtin 2-argument `iter()`, but stops on `exception`."""
    try:
        while True:
            yield function()
    except exception:
        return

class NewGenome(QtWidgets.QMainWindow):
    def __init__(self, info_path):
        super(NewGenome, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'NewGenome.ui'), self)
        self.setWindowTitle('New Genome')
        self.k = KEGG()
        self.info_path = info_path

        #---Style Modifications---#

        groupbox_style = """
        QGroupBox:title{subcontrol-origin: margin;
                        left: 10px;
                        padding: 0 5px 0 5px;}
        QGroupBox#Step1{border: 2px solid rgb(111,181,110);
                        border-radius: 9px;
                        font: 11pt "Sans Serif";
                        font: bold;
                        margin-top: 10px;}"""

        self.Step1.setStyleSheet(groupbox_style)
        self.Step2.setStyleSheet(groupbox_style.replace("Step1","Step2").replace("rgb(111,181,110)","rgb(77,158,89)"))
        self.Step3.setStyleSheet(groupbox_style.replace("Step1","Step3").replace("rgb(111,181,110)","rgb(53,121,93)"))

        #---Button Modifications---#

        self.setWindowIcon(Qt.QIcon(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "cas9image.png")))
        self.whatsthisButton.clicked.connect(self.whatsthisclicked)
        self.KeggSearchButton.clicked.connect(self.updatekegglist)
        self.resetButton.clicked.connect(self.reset)
        self.submitButton.clicked.connect(self.submit)
        self.browseForFile.clicked.connect(self.selectFasta)
        self.NCBI_File_Search.clicked.connect(self.prep_ncbi_search)
        self.JobsQueueBox.setReadOnly(True)
        self.output_browser.setText("Waiting for program initiation...")
        self.CompletedJobs.setText(" ")
        self.contButton.clicked.connect(self.continue_to_main)

        self.comboBoxEndo.currentIndexChanged.connect(self.endo_settings)

        self.runButton.clicked.connect(self.run_jobs)
        self.clearButton.clicked.connect(self.clear_job_queue)

        self.viewStatButton.setEnabled(False)

        self.JobsQueue = []  # holds Job classes.
        self.Endos = dict()
        self.file = ""

        self.process = QtCore.QProcess()
        self.process.setProcessChannelMode(QtCore.QProcess.MergedChannels)
        self.process.finished.connect(self.upon_process_finishing)
        self.seqTrans = SeqTranslate()
        self.exit = False



        self.first = False
        #show functionalities on window
        self.fillEndo()
        #self.show()

        self.num_chromo_next = False

    ####---FUNCTIONS TO RUN EACH BUTTON---####
    def selectFasta(self):

        filed = QtWidgets.QFileDialog()
        myFile = QtWidgets.QFileDialog.getOpenFileName(filed, "Choose a File")
        if (myFile[0] != ""):

            if not myFile[0].endswith(".fa") and not myFile[0].endswith(".fna") and not myFile[0].endswith(".gbff") and not myFile[0].endswith(".fasta"):
                QtWidgets.QMessageBox.question(self, "File Selection Error",
                                               "You have selected an incorrect type of file. "
                                               "Please choose a GenBank, .fasta, GBFF, or .fna file.",
                                               QtWidgets.QMessageBox.Ok)
                return
            else:
                self.file = myFile[0]
                self.s_file.append(myFile[0])
        """cdir = self.lineEdit.text()
        os.chdir(mydir)
        self.gdirectory = mydir
        print(mydir)
        print(cdir)"""

    # this function figures out which type of file the user is searching for, and then shows the
    # ncbi_search_dialog window
    # connected to the button: self.NCBI_File_Search
    def prep_ncbi_search(self):
        GlobalSettings.mainWindow.ncbi_search_dialog.searchProgressBar.setValue(0)
        GlobalSettings.mainWindow.ncbi_search_dialog.organismLineEdit.setText(self.lineEdit_1.displayText())
        GlobalSettings.mainWindow.ncbi_search_dialog.show()

    def submit(self):
        warning = ""
        if len(self.lineEdit_1.text())==0:
            warning = warning + "\nYou need to include the organism's name."
        if len(self.file) == 0:
            warning = warning + "\nYou need to select a file."
        if len(warning) !=0:
            QtWidgets.QMessageBox.information(self, "Required Information", warning, QtWidgets.QMessageBox.Ok)
            return

        if len(self.lineEdit_2.text()) == 0:
            warning = warning + "\nIt is recommended to include the organism's subspecies/strain."
        if len(self.lineEdit_3.text()) == 0:
            warning = warning + "\nYou must include an organism code (KEGG code recommended)."
        if len(warning)!=0:
            hold = QtWidgets.QMessageBox.question(self, "Missing Information", warning+
                                        "\n\nDo you wish to continue without including this information?"
                                       , QtWidgets.QMessageBox.Yes |
                                       QtWidgets.QMessageBox.No,
                                       QtWidgets.QMessageBox.No)
            if hold == QtWidgets.QMessageBox.No:
                return

        myjob = CasperJob(self.lineEdit_1.text() + " " + self.lineEdit_2.text(), self.lineEdit_2.text(),
                          self.Endos[self.comboBoxEndo.currentText()], self.lineEdit_3.text(), self.file,
                          self.tot_len_box.text(), self.seed_len_box.text(), self.pamBox.isChecked())
        self.JobsQueue.append(myjob)
        nxtLine=""
        if len(self.JobsQueueBox.toPlainText())!=0:
            nxtLine = "\n"
        self.JobsQueueBox.setPlainText(self.JobsQueueBox.toPlainText()+nxtLine+myjob.name)

    def fillEndo(self):
        f = open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "CASPERinfo"))
        while True:
            line = f.readline()
            if line.startswith('ENDONUCLEASES'):
                while True:
                    line = f.readline()
                    if(line[0]=="-"):
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
                    self.Endos[endo + "PAM: " + p_pam] = (endo, p_pam, default_seed_length, default_tot_length)

                break
        f.close()
        self.comboBoxEndo.addItems(self.Endos.keys())

    def endo_settings(self):
        # check the if it's 3' or 5', and check the box accordingly
        if int(self.seqTrans.endo_info[self.Endos[self.comboBoxEndo.currentText()][0]][3]) == 3:
            self.pamBox.setChecked(0)
        elif int(self.seqTrans.endo_info[self.Endos[self.comboBoxEndo.currentText()][0]][3]) == 5:
            self.pamBox.setChecked(1)

        self.tot_len_box.setText(self.Endos[self.comboBoxEndo.currentText()][3])
        self.seed_len_box.setText(self.Endos[self.comboBoxEndo.currentText()][2])

    def findFasta(self):
        choice = QtWidgets.QMessageBox.question(self, "Extract!", "Are you sure you want to quit?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtWidgets.QMessageBox.Yes:
            sys.exit()
        else:
            pass

    def testcheckandradio(self, state):
        if state == QtCore.Qt.Checked:
            pass

    def whatsthisclicked(self):
        QtWidgets.QMessageBox.information(self, "Organism Code", "The organism code is the manner in which CASPER will"
                                                                 " label its data files and references for the organism"
                                                                 " you are importing here. It is HIGHLY RECOMMENDED that"
                                                                 " you use the 3-4 letter code used by KEGG as this will"
                                                                 " aid in automatic accession of annotations from the"
                                                                 " database.", QtWidgets.QMessageBox.Ok)

    def updatekegglist(self):

        if len(self.lineEdit_1.text()) == 0:
            QtWidgets.QMessageBox.question(self, "Error!", "Please enter an organism into the Organism Name!", QtWidgets.QMessageBox.Ok)
            return

        self.keggSuggested.clear()
        kegg_orglist = self.k.lookfor_organism(self.lineEdit_1.text())
        holder = 0
        self.keggSuggested.setColumnCount(2)

        for item in kegg_orglist:

            second_space = item[item.find(" ") + 1:].find(" ")+item.find(" ")
            code = item[item.find(" ")+1:second_space+1]
            item = item[second_space+2:]
            semi = item.find(";")
            index = 1
            while True:
                if item[semi - index] == " ":
                    break
                index = index + 1
            organism = item[:semi - index]
            self.keggSuggested.setRowCount(holder+1)
            table_code = QtWidgets.QTableWidgetItem(code)
            table_organism = QtWidgets.QTableWidgetItem(organism)
            self.keggSuggested.setItem(holder, 0, table_organism)
            self.keggSuggested.setItem(holder, 1, table_code)
            # self.keggsearchresults.insertPlainText(item)
            holder+=1
        self.keggSuggested.resizeColumnsToContents()

    def run_jobs(self):
        if len(self.JobsQueue) > 0:
            self.progressBar.setValue(0)
            def output_stdout(p):
                line = str(p.readAll())
                line = line[2:]
                line = line[:len(line) - 1]
                for lines in filter(None, line.split(r'\n')):
                    lines.strip(r'\n')
                    if (lines == 'Finished reading in the genome file.'):
                        self.num_chromo_next = True
                    elif (self.num_chromo_next == True):
                        self.num_chromo_next = False
                        self.num_chromo = int(lines)
                    elif (lines.find('Chromosome') != -1 and lines.find('complete.') != -1):
                        temp = lines
                        temp = temp.replace('Chromosome ', '')
                        temp = temp.replace(' complete.', '')
                        temp = temp.strip(r'\\n')
                        if temp.isnumeric():
                            if (int(temp) == self.num_chromo):
                                self.progressBar.setValue(99)
                            else:
                                self.progressBar.setValue(int(temp) / self.num_chromo * 100)
                    elif (lines == 'Finished Creating File.'):
                        self.progressBar.setValue(100)

                    self.output_browser.append(lines)


            # Top layer for loop to go through all of the jobs in the queue:
            job = self.JobsQueue[0]
            program = '"' + os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'Casper_Seq_Finder" ')
            self.JobInProgress.setText(job.name)
            self.process.readyReadStandardOutput.connect(partial(output_stdout, self.process))
            program += job.get_arguments()
            print(program)
            self.process.start(program)
            self.JobsQueueBox.clear()
            for jobs in self.JobsQueue:
                if(job.name != jobs.name):
                    self.JobsQueueBox.append(jobs.name)
        else:
            error = QtWidgets.QMessageBox.critical(self, "No Jobs To Run",
                                                   "No jobs are in the Queue to run. Please submit a job before Running.",
                                                   QtWidgets.QMessageBox.Ok)

    def upon_process_finishing(self):
        self.CompletedJobs.append(self.JobsQueue[0].name)
        self.JobsQueue.pop(0)
        self.process.close()
        self.num_chromo = 0
        if len(self.JobsQueue) != 0:
            self.progressBar.setValue(0)
            self.run_jobs()
        else:
            self.JobInProgress.clear()
            self.JobsQueue = []
            self.JobsQueueBox.clear()
            self.output_browser.clear()

    def clear_job_queue(self):
        self.process.kill()
        self.JobsQueue = []
        self.JobsQueueBox.clear()
        self.lineEdit_1.clear()
        self.lineEdit_2.clear()
        self.lineEdit_3.clear()
        self.keggSuggested.setRowCount(0)
        self.output_browser.clear()
        self.JobInProgress.clear()
        self.CompletedJobs.clear()
        self.s_file.clear()
        self.progressBar.setValue(0)
        self.first = False
        self.s_file.clear()

    def reset(self):
        self.lineEdit_1.clear()
        self.lineEdit_2.clear()
        self.lineEdit_3.clear()
        self.keggSuggested.clear()
        self.first = False
        self.s_file.clear()

    def closeEvent(self, event):
        # make sure that there are cspr files in the DB
        file_names = os.listdir(GlobalSettings.CSPR_DB)
        noCSPRFiles = True
        for file in file_names:
            if 'cspr' in file:
                noCSPRFiles = False
                break
        if noCSPRFiles == True:
            error = QtWidgets.QMessageBox.question(self, "No CSPR File generated",
                                                    "No CSPR file has been generated, thus the main program cannot run. Please create a CSPR file."
                                                    "Alternatively, you could quit the program. Would you like to quit?",
                                                    QtWidgets.QMessageBox.Yes |
                                                    QtWidgets.QMessageBox.No,
                                                    QtWidgets.QMessageBox.No)
            if (error == QtWidgets.QMessageBox.No):
                event.ignore()
                return
            else:
                self.exit = False
                event.accept()

        else:
            self.process.kill()
            self.JobsQueue = []
            self.JobsQueueBox.clear()
            self.lineEdit_1.clear()
            self.lineEdit_2.clear()
            self.lineEdit_3.clear()
            self.keggSuggested.setRowCount(0)
            self.output_browser.clear()
            self.JobInProgress.clear()
            self.CompletedJobs.clear()
            self.s_file.clear()
            self.progressBar.setValue(0)
            self.first = False
            GlobalSettings.CASPER_FOLDER_LOCATION = self.info_path
            GlobalSettings.mainWindow.show()
            if GlobalSettings.mainWindow.orgChoice.currentText() != '':
                GlobalSettings.mainWindow.orgChoice.currentIndexChanged.disconnect()
            GlobalSettings.mainWindow.orgChoice.clear()
            GlobalSettings.mainWindow.endoChoice.clear()
            GlobalSettings.mainWindow.getData()
            GlobalSettings.MTWin.launch(GlobalSettings.CSPR_DB)
            GlobalSettings.pop_Analysis.launch(GlobalSettings.CSPR_DB)
            event.accept()

    def continue_to_main(self):
        # make sure that there are cspr files in the DB
        file_names = os.listdir(GlobalSettings.CSPR_DB)
        noCSPRFiles = True
        for file in file_names:
            if 'cspr' in file:
                noCSPRFiles = False
                break
        if noCSPRFiles == True:

            error = QtWidgets.QMessageBox.question(self, "No CSPR File generated",
                                                   "No CSPR file has been generated, thus the main program cannot run. Please create a CSPR file."
                                                   "Alternatively, you could quit the program. Would you like to quit?",
                                                   QtWidgets.QMessageBox.Yes |
                                                   QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.Yes)

            if (error == QtWidgets.QMessageBox.Yes):
                self.exit = True
                self.close()

        else:
            self.process.kill()
            self.JobsQueue = []
            self.JobsQueueBox.clear()
            self.lineEdit_1.clear()
            self.lineEdit_2.clear()
            self.lineEdit_3.clear()
            self.keggSuggested.setRowCount(0)
            self.output_browser.clear()
            self.JobInProgress.clear()
            self.CompletedJobs.clear()
            self.s_file.clear()
            self.progressBar.setValue(0)
            self.first = False
            GlobalSettings.CASPER_FOLDER_LOCATION = self.info_path
            GlobalSettings.mainWindow.show()
            if GlobalSettings.mainWindow.orgChoice.currentText() != '':
                GlobalSettings.mainWindow.orgChoice.currentIndexChanged.disconnect()
            GlobalSettings.mainWindow.orgChoice.clear()
            GlobalSettings.mainWindow.endoChoice.clear()
            GlobalSettings.mainWindow.getData()
            GlobalSettings.MTWin.launch(GlobalSettings.CSPR_DB)
            GlobalSettings.pop_Analysis.launch(GlobalSettings.CSPR_DB)
            self.hide()


class CasperJob:
    def __init__(self, org, suborg, endo, org_code, ref_file, tot_len, seed_len, pamdir):
        self.name = endo[0] + " targets in " + org
        self.organism_name = org
        self.substrain = suborg
        self.organism_code = org_code
        self.endo_name = endo[0]
        self.endo_pam = endo[1]
        self.reference_file = ref_file

        # These are endonuclease specific settings that should be pulled from CASPERinfo
        self.anti = pamdir
        self.sequence_length = tot_len
        self.seed_length = seed_len

    def get_arguments(self):
        db_location = GlobalSettings.CSPR_DB
        if(GlobalSettings.OPERATING_SYSTEM_ID == "Windows"):
            db_location = db_location.replace('/','\\')
            ref = str(self.reference_file).replace('/','\\')
        else:
            db_location = db_location
            ref = str(self.reference_file)
        cmd = str()
        cmd += '"' + str(self.endo_name) + '" '
        cmd += '"' + str(self.endo_pam) + '" '
        if(self.organism_code == ""):
            cmd += '""'
        else:
            cmd += '"' + str(self.organism_code) + '" '

        if self.anti:
            cmd += '"' + "TRUE" + '" '
        else:
            cmd += '"' + "FALSE" + '" '


            cmd += '"' + db_location + "/" + '" '
            cmd += '"' + os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "CASPERinfo") + '" '

        cmd += '"' + ref + '" '
        cmd += '"' + self.organism_name + '" '
        cmd += '"' + self.sequence_length + '" '
        cmd += '"' + self.seed_length + '" '
        if(self.substrain == ""):
            cmd += '" "'
        else:
            cmd += '"' + self.substrain + '"'


        print(cmd)

        if (GlobalSettings.OPERATING_SYSTEM_ID == "Windows"):
            db_location = db_location.replace('/', '\\')
            ref = str(self.reference_file).replace('/', '\\')
        else:
            db_location = db_location
            ref = str(self.reference_file)
        ret_array = [self.endo_name, self.endo_pam, self.organism_code]
        # attach the 5' or 3' direction
        if self.anti:
            ret_array.append("TRUE")
        else:
            ret_array.append("FALSE")

        ret_array.append(GlobalSettings.CSPR_DB + "/")
        ret_array.append(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "CASPERinfo"))
        ret_array.append(self.reference_file)

        ret_array.append(self.organism_name)
        ret_array.append(self.sequence_length)
        ret_array.append(self.seed_length)
        ret_array.append(self.substrain)

        return cmd
