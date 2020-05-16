import sys, os
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QDir, QUrl

def main():
	createGraph()

def createGraph(self):
	
	self.browser = QWebEngineView()
	file_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "kishanTestHtml.html"))
	local_url = QUrl.fromLocalFile(file_path)
	self.browser.load(local_url)
	self.browser.show()


if __name__ == "__main__":
	main()
