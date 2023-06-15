   # Copyright (C) 2023  Sergio Frasca
   #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

'''
      Module GUISNAG_QT
GUI QT features for SnagPy


'''

import sys
from PyQt6.QtWidgets import QApplication,QPushButton,QWidget, QTableWidget, QTableWidgetItem, QVBoxLayout
import BASIC



class TableWindow(QWidget):
    '''
    Class for GUI table for snag_table or simp_dict
    '''
    def __init__(self,stordic):
        super().__init__()

        # global qt_data
        if isinstance(stordic,BASIC.snag_table):
            self.setWindowTitle(stordic.name)
            stab=stordic
            print('snag_table with titles ',stab.titles)
            nr=stab.nr
            nc=stab.nc
            titles=stab.titles
            qt_data=stab.data[1:]
        elif isinstance(stordic,BASIC.simp_dict) or isinstance(stordic,dict): 
            if isinstance(stordic,dict):
                print('dictionary converted to simple_dict')
                stordic=BASIC.simp_dict(stordic)
            self.setWindowTitle(stordic.name)
            sdic=stordic
            print('simp_dict with length ',sdic.len)
            nr=sdic.len
            nc=2
            keys=list(sdic.keys)
            values=list(sdic.values)
            titles=['keys','values']
            qt_data=[]
            for i in range(nr):
                qt_data.append([keys[i],values[i]])
        else:
            print('*** ERROR : ',type(stordic))
        
        # print(qt_data)

        # Creazione del layout verticale per la finestra
        layout = QVBoxLayout()

        self.push = QPushButton('Close')

        # Creazione del widget QTableWidget
        self.table = QTableWidget()

        # Impostazione del numero di righe e colonne
        self.table.setRowCount(nr)
        self.table.setColumnCount(nc)

        # Impostazione dei nomi delle colonne
        self.table.setHorizontalHeaderLabels(titles)

        # Impostazione dei dati delle celle
        for row in range(nr):
            for col in range(nc):
                dat=str(qt_data[row][col])
                item = QTableWidgetItem(dat)
                self.table.setItem(row, col, item)

        # Aggiunta del widget QTableWidget al layout 
        
        layout.addWidget(self.table)
        layout.addWidget(self.push)

        self.push.clicked.connect(self.close_table)

        self.setLayout(layout)
        self.setWindowTitle("Table Example")
        self.show()


# tableWidget.cellChanged.connect(on_cell_changed)

    def close_table(self):
        self.close()
        


def crea_table(stordic):
    app = QApplication(sys.argv)
    window = TableWindow(stordic)
    window.show()
    app.exec()


# def on_cell_changed(row, column):
#     item = tableWidget.item(row, column)
#     value = item.text()
#     # Aggiorna i dati nella struttura dati
#     qt_data[row][column] = value


import sys
from PyQt6.QtWidgets import QApplication, QMainWindow, QTableWidget, QTableWidgetItem
from PyQt6.QtCore import Qt

import BASIC

def gui_table(stordic):
    '''
    Offers GUI for snag_table and simp_dict 
    '''
    if isinstance(stordic, BASIC.snag_table):
        typ = 1
        stab = stordic
        print('snag_table with titles', stab.titles)
        name = stab.name
        nr = stab.nr
        nc = stab.nc
        titles = stab.titles
        data = stab.data[1:]
    elif isinstance(stordic, BASIC.simp_dict) or isinstance(stordic, dict):
        if isinstance(stordic, dict):
            print('dictionary converted to simple_dict')
            stordic = BASIC.simp_dict(stordic)
        sdic = stordic
        print('simp_dict with length', sdic.len)
        name = sdic.name
        nr = sdic.len
        nc = 2
        keys = list(sdic.keys)
        values = list(sdic.values)
        titles = ['keys', 'values']
    else:
        print('*** ERROR:', type(stordic))
    
    print('table with nr, nc', nr, nc)
    print(len(data), len(data[0]))
    
    app = QApplication(sys.argv)
    window = QMainWindow()
    window.setWindowTitle(name)
    
    table = QTableWidget(window)
    table.setRowCount(nr)
    table.setColumnCount(nc)
    table.setHorizontalHeaderLabels(titles)
    
    for row in range(nr):
        for col in range(nc):
            item = QTableWidgetItem(str(data[row][col]))
            table.setItem(row, col, item)
    
    table.setEditTriggers(QTableWidget.EditTriggers.NoEditTriggers)
    table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
    table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
    table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
    
    window.setCentralWidget(table)
    window.show()
    
    sys.exit(app.exec())
    
