import sys, os
import GlobalSettings
from PyQt5 import QtWidgets, uic, QtGui, QtCore, Qt
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QDir, QUrl
from Bio import Entrez, SeqIO
from bioservices.kegg import KEGG

def main():
	createGraph()


def splitStringNCBI(longString):
	try:
		return (longString).split(':')[2]
	except:
		print(longString.split(':'))
		print("Please search for a genome")

def splitStringLocal(longString):
	try:
		print(longString)
		print(longString.split('/'))
		return (longString.split('/').pop()).split('.')[0]
	except:
		print(longString.split('/').pop()).split('.')[0]
		print("Please search for a genome")

def ncbiAPI(filename):

	gb_file = filename
	genomeList = []
	for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
		genomeList.append(gb_record.id)

	return genomeList


def createHtml(genomeList):

	htmlString1 = """
	<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
	<ihtml>
	<head>
		<title>NCBI Sequence Viewer - programmatic initialization</title>
		<script type="text/javascript" src="https://www.ncbi.nlm.nih.gov/projects/sviewer/js/sviewer.js"></script>
		<script>
			function loadSV(id) {
				var svapp = SeqView.App.findAppByDivId("mySeqViewer1");
				if (!svapp)
					svapp = new SeqView.App("mySeqViewer1");
				params = 'appname=UniversityofTennesseeCasper&amp;id=' + id;
				svapp.reload(params);
			}
		</script>
	</head>
	<body>
		Select a chromosome from the list:<br/>

	<select onchange="loadSV(event.target.value);">
		<option value="">-</option>
	"""

	htmlString2 = """
		<option value="{}">Chromosome {}</option>
	"""

	htmlString3 = """
		</select>
		<br/>
	
		<div id="mySeqViewer1" class="SeqViewerApp" data-autoload> 
			<a href="?embedded=true&appname=testapp1&id={}"></a> 
		</div>	
		
	</body>
	</html>
	""".format(genomeList[0])

	# Find file path for template
	# seek to beginning and truncate
	genomeBrowserTemplateFilePath = os.path.abspath(os.path.join(os.path.dirname(__file__), "genomeBrowserTemplate.html"))

	raw = open(genomeBrowserTemplateFilePath, "r+")
	raw.seek(0)
	raw.truncate()

	#write the 3 part string format
	raw.write(htmlString1)

	for index,genome in enumerate(genomeList):
		raw.write(htmlString2.format(genome,index+1))

	raw.write(htmlString3)


def createGraph(self):
	annotationWindow = self.findChild(QtWidgets.QComboBox,'Annotations_Organism')
	NCBIFileBoolean = self.findChild(QtWidgets.QRadioButton, 'NCBI_Select')
	localFileBoolean = self.findChild(QtWidgets.QRadioButton, 'Annotation_Ownfile')

	selectedGenome = annotationWindow.currentText()

	#returns gci String value that can be used to retrieve the file
	print("selected Genome is ", selectedGenome)
	gciVariable = ""

	if(localFileBoolean.isChecked() == True):
		print("local")
		gciVariable = splitStringLocal(selectedGenome)
	elif(NCBIFileBoolean.isChecked() == True):
		print("NCBI")
		gciVariable = splitStringNCBI(selectedGenome)

	print("gci variable is ", gciVariable)

	if(gciVariable == None):
		print("Please select NCBI Genome")
		return

	directory = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])))
	onlyfiles = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
	for file in onlyfiles:
		print("One file is ",file)
		if gciVariable in file and ".fna" not in file:
			fileToSearch = file
			print("searching file:", fileToSearch)
			break

	try:
		print("file to search is ", fileToSearch)
		genomeList = ncbiAPI(fileToSearch)

	except:
		print("No gbff file found")
		QtWidgets.QMessageBox.question(self, "GBFF_FileNotFound", "GBFF file is not in selected directory", QtWidgets.QMessageBox.Ok)
		return

	createHtml(genomeList)
	self.browser = QWebEngineView()
	file_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "genomeBrowserTemplate.html"))
	local_url = QUrl.fromLocalFile(file_path)
	self.browser.load(local_url)
	self.browser.show()


if __name__ == "__main__":
	main()
