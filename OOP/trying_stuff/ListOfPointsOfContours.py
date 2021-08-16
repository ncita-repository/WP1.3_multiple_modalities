# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 14:30:29 2020

@author: ctorti
"""



""" 
Create list of contour points containing N points from up-sampled Contour1
(USContour1) followed by N points from the modified up-sampled Contour2 
(NewUPContour2), ready to use as input to Transformix, where:

- USContour1 is up-sampled Contour1
- USContour2 is up-sampled Contour2
- NewUSContour2 is re-ordered points from USContour2 that minimises surface 
area formed by joining points in USContour1 and USContour2

(Combining all points from both contours into a single list ready to transform.)
"""

def ListOfPointsOfContours(Contour1, Contour2, N):
    # Import package:
    import time
    
    # Import functions:
    from UpsampleContour import UpsampleContour
    from FindIndForMinCumLength import FindIndForMinCumLength
    from ShiftContourPoints import ShiftContourPoints


    print(f'\nThere are {len(Contour1)} points in Contour1 and ',
          f'{len(Contour2)} points in Contour2')
    
    
    
    # Start timing:
    times = []
    times.append(time.time())
    
    if N > 1:
        print(f'\nUp-sampling contours by {N}x...')
        
        # Up-sample Contour1 and Contour2:
        Contour1 = UpsampleContour(Contour=Contour1, N=N)
        Contour2 = UpsampleContour(Contour=Contour2, N=N)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        print(f'Took {Dtime} s to up-sample contours.')
    
    
    print(f'\nMinimising area...')
    
    # Re-order Contour2:
    ROContour2 = ShiftContourPoints(Contour=Contour2, 
                                    StartInd=FindIndForMinCumLength(Contour1=Contour1,
                                                                    Contour2=Contour2))
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'Took {Dtime} s to minimise area and re-order Contour2',
          f'({round(Dtime/len(Contour2), 3)} s per point).')
    
        
    # Initialise list of points, then extend USContour1 and ROUSContour2:
    Points = []
    Points.extend(Contour1)
    Points.extend(ROContour2)
    
    print(f'\nThere are {len(Contour1)} points in Contour1, {len(Contour2)} in Contour2',
          f'and {len(Points)} in the combined list of points.')
    
    return Points