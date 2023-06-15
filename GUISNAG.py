   # Copyright (C) 2023  Sergio Frasca
   #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

'''
      Module GUISNAG
GUI features for SnagPy

Sections:
> Cursor                -> dummy_curs
> MPL toolbox           -> dummy_tbox
'''

"""
=================
Cross-hair cursor
=================

This example adds a cross-hair as a data cursor.  The cross-hair is
implemented as regular line objects that are updated on mouse move.

We show three implementations:

1) A simple cursor implementation that redraws the figure on every mouse move.
   This is a bit slow, and you may notice some lag of the cross-hair movement.
2) A cursor that uses blitting for speedup of the rendering.
3) A cursor that snaps to data points.

Faster cursoring is possible using native GUI drawing, as in
:doc:`/gallery/user_interfaces/wxcursor_demo_sgskip`.

The mpldatacursor__ and mplcursors__ third-party packages can be used to
achieve a similar effect.

__ https://github.com/joferkington/mpldatacursor
__ https://github.com/anntzer/mplcursors

.. redirect-from:: /gallery/misc/cursor_demo
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pylab as pyl
from matplotlib.backend_bases import MouseEvent
import copy

import wx
import wx.grid

def get_posclick(n=1):
    '''mouse click position
    '''
    pos=plt.ginput(n=n)

    return pos


class SnappingCursor:
    """
    A cross-hair cursor that snaps to the data point of a line, which is
    closest to the *x* position of the cursor.

    For simplicity, this assumes that *x* values of the data are sorted.
    """
    def __init__(self, ax, line):
        self.ax = ax
        self.horizontal_line = ax.axhline(color='k', lw=0.8, ls='--')
        self.vertical_line = ax.axvline(color='k', lw=0.8, ls='--')
        self.x, self.y = line.get_data()
        self._last_index = None
        # text location in axes coords
        self.text = ax.text(0.72, 0.9, '', transform=ax.transAxes)

    def set_cross_hair_visible(self, visible):
        need_redraw = self.horizontal_line.get_visible() != visible
        self.horizontal_line.set_visible(visible)
        self.vertical_line.set_visible(visible)
        self.text.set_visible(visible)
        return need_redraw

    def on_mouse_move(self, event):
        if not event.inaxes:
            self._last_index = None
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.draw()
        else:
            self.set_cross_hair_visible(True)
            x, y = event.xdata, event.ydata
            index = min(np.searchsorted(self.x, x), len(self.x) - 1)
            if index == self._last_index:
                return  # still on the same data point. Nothing to do.
            self._last_index = index
            x = self.x[index]
            y = self.y[index]
            # update the line positions
            self.horizontal_line.set_ydata([y])
            self.vertical_line.set_xdata([x])
            self.text.set_text('x=%1.2f, y=%1.2f' % (x, y))
            self.ax.figure.canvas.draw()


def snap_curs(x,y):
    fig, ax = plt.subplots()
    ax.set_title('Snapping cursor')
    line, = ax.plot(x, y, 'o')
    snap_cursor = SnappingCursor(ax, line)
    fig.canvas.mpl_connect('motion_notify_event', snap_cursor.on_mouse_move)
    t = ax.transData
    MouseEvent("motion_notify_event", ax.figure.canvas, *t.transform((0.5, 0.5)))._process()
    
    plt.show()


data=[]

class STabFrame(wx.Frame):
    def __init__(self, parent, stab):
        print(stab.titles)
        nr=stab.nr
        nc=stab.nc
        titles=stab.titles
        data=stab.data[1:]

        super(STabFrame, self).__init__(parent)

        # Creazione del pannello principale
        panel = wx.Panel(self)

        # Creazione di un oggetto wx.grid.Grid per visualizzare la tabella
        table_grid = wx.grid.Grid(panel)

        # Imposta il numero di righe e colonne della tabella
        table_grid.CreateGrid(numRows=nr, numCols=nc)

        # Imposta i nomi delle colonne
        ii=0
        for cl in titles:
            table_grid.SetColLabelValue(ii, cl)
            ii+=1

        # Imposta i dati delle celle
        for row in range(nr):
            for col in range(nc):
                # print(row,col)
                # table_grid.SetCellValue(row, col, f"Riga {row+1}, Colonna {col+1}")
                table_grid.SetCellValue(row, col,str( data[row][col]))

        # Imposta la dimensione automatica delle colonne
        table_grid.AutoSizeColumns()

        # Imposta il sizer per adattare la tabella al pannello
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(table_grid, 1, wx.EXPAND)
        panel.SetSizer(sizer)

        # Mostra la finestra
        self.Show()


def creatable(stab):
    app = wx.App()
    
    app.MainLoop()


def on_cell_changed(event):
    row = event.GetRow()
    col = event.GetCol()
    value = wx.grid.GetCellValue(row, col)
    # Aggiorna i dati nella struttura dati
    data[row+1][col] = value


# frame = STabFrame(None,stab)
# grid = wx.grid.Grid(frame)
# grid.Bind(wx.grid.EVT_GRID_CELL_CHANGED, on_cell_changed)


def go_stab(stab):
    data=copy.copy(stab.data)

    creatable(stab)

    stab.data=data
    return stab


def wx_stable(stab):
    '''
    GUI for snag_table
    '''