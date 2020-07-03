# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 17:32:00 2020

@author: ctorti
"""



""" Plot result of interpolation """


def PlotInterpolation(Points1, Points2, Interp1to2, Interp2to1, ExportPlot):
    import matplotlib.pyplot as plt
    import time
    
    
    MarkerSize = 5

    # Create a figure with two subplots and the specified size:
    plt.subplots(1, 1, figsize=(14, 14))
    
    
    # Group all points:
    Points = [Points1, Points2, Interp1to2, Interp2to1]
    
    # Plot colours:
    colours = ['blue', 'green', 'red', 'magenta']
    
    for i in range(len(Points)):
        # Unpack tuple and store each x,y tuple in arrays 
        # X and Y:
        X = []
        Y = []
        
        for x, y, z in Points[i]:
            X.append(x)
            Y.append(y)
        
        plt.plot(X, Y, linewidth=1, c=colours[i]);
        plt.plot(X, Y, '.', markersize=MarkerSize, c=colours[i]);
    
    
    plt.title(f'Interpolation (red and magenta) of contours in blue and green')
    
    # Create filename for exported figure:
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + f'_Interpolation_of_contours_both_directions.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
        
        print('\nPlot exported as', FigFname)
        
        
    return