# -*- coding: utf-8 -*-

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from mueller_matrices import stokes2ellipse


def plotEllipse(ellipt, orient, size=5, axes=None):
    def add_arrow(line, size=20, color='k'):
        """
        add an arrow to a line.

        line:       Line2D object
        position:   x-position of the arrow. If None, mean of xdata is taken
        direction:  'left' or 'right'
        size:       size of the arrow in fontsize points
        color:      if None, line color is taken.

        adapted from: https://stackoverflow.com/questions/34017866/arrow-on-a-line-plot-with-matplotlib
        """
        if color is None:
            color = line.get_color()

        xdata = line.get_xdata()
        ydata = line.get_ydata()

        position = xdata.mean()
        # find closest index
        start_ind = np.argmin(np.absolute(xdata - position))
        end_ind = start_ind - 1

        line.axes.annotate('',
                           xytext=(xdata[start_ind], ydata[start_ind]),
                           xy=(xdata[end_ind], ydata[end_ind]),
                           arrowprops=dict(arrowstyle="->", color=color),
                           size=size)

    theta = np.linspace(0, 2 * np.pi, 360)
    rotmat = [[np.cos(orient), -np.sin(orient)], [np.sin(orient), np.cos(orient)]]

    x = np.cos(theta)
    y = ellipt * np.sin(theta)
    (xr, yr) = np.matmul(rotmat, [x.flatten("C"), y.flatten("C")])

    if axes == None:
        plt.figure(figsize=[size, size])
        axes = plt.axes()
    else:
        plt.sca(axes)
    axes.axis('square')
    line = plt.plot(xr, yr)[0]
    add_arrow(line)
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    
    
def plotEllipseStokes(S0, S1, S2, S3, size=5, axes=None):
    def add_arrow(line, size=20, color='k'):
        """
        add an arrow to a line.

        line:       Line2D object
        position:   x-position of the arrow. If None, mean of xdata is taken
        direction:  'left' or 'right'
        size:       size of the arrow in fontsize points
        color:      if None, line color is taken.

        adapted from: https://stackoverflow.com/questions/34017866/arrow-on-a-line-plot-with-matplotlib
        """
        if color is None:
            color = line.get_color()

        xdata = line.get_xdata()
        ydata = line.get_ydata()

        position = xdata.mean()
        # find closest index
        start_ind = np.argmin(np.absolute(xdata - position))
        end_ind = start_ind - 1

        line.axes.annotate('',
                           xytext=(xdata[start_ind], ydata[start_ind]),
                           xy=(xdata[end_ind], ydata[end_ind]),
                           arrowprops=dict(arrowstyle="->", color=color),
                           size=size)

    Stokes = np.array([S0,S1,S2,S3])
    
    ellipt, orient = stokes2ellipse(Stokes)
    
    theta = np.linspace(0, 2 * np.pi, 360)
    rotmat = [[np.cos(orient), -np.sin(orient)], [np.sin(orient), np.cos(orient)]]

    x = np.cos(theta)
    y = ellipt * np.sin(theta)
    (xr, yr) = np.matmul(rotmat, [x.flatten("C"), y.flatten("C")])

    if axes == None:
        plt.figure(figsize=[size, size])
        axes = plt.axes()
    else:
        plt.sca(axes)
    axes.axis('square')
    line = plt.plot(xr, yr)[0]
    add_arrow(line)
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)