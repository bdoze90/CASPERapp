from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvas
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.figure import Figure

class pop_analysis_3dgraph(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.canvas = FigureCanvas(Figure(figsize=(0.3,0.8)))
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        self.canvas.axes = self.canvas.figure.add_subplot(111, projection='3d')
        self.canvas.figure.subplots_adjust(left=0.0, right=0.9, bottom=0.3, top=0.9)
        self.setLayout(vertical_layout)
        # self.canvas.figure.tight_layout()