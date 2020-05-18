
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
import GlobalSettings

######################################################
# CoTargeting class: This class is a window that lets the user select which endonucleases to co-target with
# inputs are from the user and from results
# from results: the organism name and the list of endonucleases for that organism
# from user: which endonucleases to co-target
######################################################
class CoTargeting(QtWidgets.QDialog):
    def __init__(self, path):
        # pyqt stuff
        super(CoTargeting, self).__init__()
        uic.loadUi('co_targeting.ui', self)
        self.setWindowIcon(QtGui.QIcon("cas9image.png"))

        # endo_table stuff
        self.endo_table.setColumnCount(1)  # hardcoded because there will always be 1 columns
        self.endo_table.setShowGrid(True)
        self.endo_table.setHorizontalHeaderLabels("Endonuclease;".split(";"))
        self.endo_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.endo_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.endo_table.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

        # variables
        self.info_path = path

        # button connections
        self.cancel_button.clicked.connect(self.cancel_function)
        self.submit_button.clicked.connect(self.submission_function)

    # launches the window
    # it is expecting endo_choices in the form of a list, and the orgName in the form of a string
    # it sets the organism name, and sets the table as well
    def launch(self, endo_choices, orgName):
        self.orgName.setText(orgName)
        setTableList = list()

        # only get the endo choices that were original
        for item in endo_choices:
            checkList = item.split(",")
            if len(checkList) == 1:
                setTableList.append(item)

        # go through and set each table item, but also set the row count
        self.endo_table.setRowCount(len(setTableList))
        loopCount = 0
        for item in setTableList:
            tabWidget = QtWidgets.QTableWidgetItem(item)
            self.endo_table.setItem(loopCount, 0, tabWidget)
            loopCount += 1
        self.endo_table.resizeColumnsToContents()

        # now show
        self.show()

    # this is the cancel function
    # it clears the table, sets the organism to nothing, and hides the window
    def cancel_function(self):
        self.endo_table.clearContents()
        self.endo_table.setRowCount(0)
        self.orgName.setText("")
        self.hide()

    # this is the submission function.
    # it makes sure the user selects at least 2 endonucleases
    # then it goes through and returns the endonucleases selected
    # once it gets those, it then calls a function in Results that repopulates the table correctly
    def submission_function(self):
        # set the selected_list, and make sure they select at least 2 endonucleases
        selected_list = self.endo_table.selectedItems()
        if len(selected_list) <= 1:
            QtWidgets.QMessageBox.question(self, "Nothing Selected", "No endonucleases selected. Please select at least 2 endonucleases.",
                                           QtWidgets.QMessageBox.Ok)
            return

        # go through and get which endonuclease's have been selected
        ret_endo_list = list()
        for i in range(self.endo_table.rowCount()):
            if self.endo_table.item(i, 0).isSelected():
                ret_endo_list.append(self.endo_table.item(i, 0).text())

        GlobalSettings.mainWindow.Results.co_target_endo_list = ret_endo_list
        GlobalSettings.mainWindow.Results.populate_cotarget_table()
        self.cancel_function()

