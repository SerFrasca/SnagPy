   # Copyright (C) 2023  Sergio Frasca
   #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

'''
      Module GUISNAG_TK
GUI features for SnagPy version Tk

Sections:
'''

import tkinter as tk
import tksheet
import copy

import BASIC



def on_cell_change(event):
    get_edited_data()

def get_edited_data():
    global sheet_9399,dat_9399
    sheet=sheet_9399
    edited_data = sheet.get_sheet_data()
    dat_9399=edited_data
    print("Dati modificati:")
    print(dat_9399)


def gui_table(stordic):
    '''
    Offers GUI for snag_table and simp_dict 
    '''
    global sheet_9399,dat_9399
    if isinstance(stordic,BASIC.snag_table):
        typ=1
        stab=stordic
        print('snag_table with titles ',stab.titles)
        name=stab.name
        nr=stab.nr
        nc=stab.nc
        titles=stab.titles
        data=stab.data[1:]
    elif isinstance(stordic,BASIC.simp_dict) or isinstance(stordic,dict): 
        if isinstance(stordic,dict):
            print('dictionary converted to simple_dict')
            stordic=BASIC.simp_dict(stordic)
        sdic=stordic
        print('simp_dict with length ',sdic.len)
        name=sdic.name
        nr=sdic.len
        nc=2
        keys=list(sdic.keys)
        values=list(sdic.values)
        titles=['keys','values']
            # qt_data=[]
            # for i in range(nr):
            #     qt_data.append([keys[i],values[i]])
    else:
        print('*** ERROR : ',type(stordic))
    print('table with nr,nc',nr,nc)
    print(len(data),len(data[0]))
    root = tk.Tk()
    root.title=name
    sheet = tksheet.Sheet(root)
    sheet.pack(fill="both", expand=True)
    
    # sheet.grid()

    sheet.set_sheet_data([[data[ri][cj] for cj in range(nc)] for ri in range(nr)])# table enable choices listed below:

    sheet.enable_bindings(("single_select",

                       "row_select",

                       "column_width_resize",

                       "arrowkeys",

                       "right_click_popup_menu",

                       "rc_select",

                       "rc_insert_row",

                       "rc_delete_row",

                       "copy",

                       "cut",

                       "paste",

                       "delete",

                       "undo",

                       "edit_cell"))

    # sheet.set_sheet_data(data=data)
    sheet_9399=sheet
    sheet.extra_bindings([("cell_select", on_cell_change)])

    root.mainloop()

    print('ciao')

    return dat_9399

