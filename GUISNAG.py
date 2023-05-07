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


def get_posclick(n=1):
    '''mouse click position
    '''
    pos=plt.ginput(n=n)

    return pos

