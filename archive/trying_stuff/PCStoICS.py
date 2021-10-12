# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 13:58:43 2020

@author: ctorti
"""





"""
Function:
    PCStoICS()
    
Purpose:
    Convert points Pts_PCS from Patient Coordinate System (PCS) to Image 
    Coordinate System (ICS).
    
    Pts_PCS are in the form:
        
    [[point0_x, point0_y, point0_z],
     [point1_x, point1_y, point1_z],
     ...
     ]
    
    The conversion requires the image origin (ImagePositionPatient), the 
    direction cosines along x, y and z, and the pixel spacings d = [di, dj, dk].


Input:
    Pts_PCS    - Points in Patient Coordinate System 
    
    Origin     - The image origin (or ImagePositionPatient) (e.g. [x0, y0, z0])
    
    Directions - The direction cosine along x (rows), y (columns) and z (slices)
                 (e.g. [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz])
    
    Spacings   - The pixel spacings along x, y and z (e.g. [di, dj, dk])
    
    
    
Returns:
    Pts_PCS - Points in Image Coordinate System
    



Equations:

P = S + (di*i).X + (dj*j).Y + (dk*k).Z

where P = (Px, Py, Pz) is the point in the Patient Coordinate System

      S = (Sx, Sy, Sz) is the origin in the Patient Coordinate System (i.e.
      the ImagePositionPatient)
      
      X = (Xx, Xy, Xz) is the direction cosine vector along rows (x) 
      
      Y = (Yx, Yy, Yz) is the direction cosine vector along columns (y)
      
      Z = (Zx, Zy, Zz) is the direction cosine vector along slices (z)
      
      di, dj and dk are the pixel spacings along x, y and z
      
      i, j, and k are the row (x), column (y) and slice (z) indices
      
      
Solve for i, j and k, i.e.:
    
    Vx = Xx*di*i + Yx*dj*j + Zx*dk*k
    
    Vy = Xy*di*i + Yy*dj*j + Zy*dk*k
    
    Vz = Xz*di*i + Yz*dk*k + Zz*dk*k
    
where Vx = Px - Sx

      Vy = Py - Sy
      
      Vz = Pz - Sz
      
Solving for i, j and k results in the expressions:
    
    i = (1/di)*d/e
        
    j = (1/dj)*( (a/b)*di*i + c/b )
    
    k = (1/dk)*(1/Zz)*(Vz - Xz*di*i - Yz*dj*j)
    or
    k = (1/dk)*(1/Zz)*( Vz - (c/b)*Yz - (Xz + (a/b)*Yz)*di*i )
    
    (either expression for k should be fine)
    
In the case where X = [1, 0, 0]
                  Y = [0, 1, 0]
                  Z = [0, 0, 1]
                  
The expressions simplify to:
    
    a = 0
    
    b = 1
    
    c = Vy
    
    d = Vx
    
    e = 1
    
    i = Vx/di
    
    j = Vy/dj
    
    k = Vz/dk
    
"""






""" Note: 
    
    The input Pts_PCS for PCStoICS() needs to be an array of points, 
    so if only converting a single point, the Pts_PCS must have 
    square brackets around it to turn it into an array of length 1,
    
    e.g. Pts_PCS = [[x,y,z]] for point = [x,y,z]
    
    but the output must be indexed with [0],
    
    e.g. Pts_ICS[0]
    
    to get the point in the form [x,y,z] (otherwise it will be [[x,y,z]]).
    
    
"""
        

def PCStoICS(Pts_PCS, Origin, Directions, Spacings):
    # Import packages:
    import numpy as np
    #import SimpleITK as sitk
    
    """ 
    Depending on whether Pts_PCS is a list of points or a single point things
    will need to be dealt with differently. So first determine the number of
    points in Pts_PCS.
    """
    
    # If Pts_PCS is a list of points, the first item will be a list of 
    # [x, y, z] coordinates:
    #if isinstance(Pts_PCS[0], list):
    #    # Pts_PCS is a list of P points:
    #    P = len(Pts_PCS)
    #    
    #else:
    #    # Pts_PCS is a single point:
    #    P = 1
    #    
    #    # Since the for loop below is intended to act on a list of points, 
    #    # convert the point into a list of length 1:
    #    Pts_PCS = [Pts_PCS]
        
        
    # Get the shape of the list of points:
    ListShape = np.array(Pts_PCS).shape
    
    #print('\nListShape =', ListShape)
    
    #  Get the length of ListShape:
    LenListShape = len(ListShape)
    
    #print('\nLenListShape =', LenListShape)
    
    if LenListShape == 1:
        # Pts_PCS is a single point. Turn into a list with one point:
        Pts_PCS = [Pts_PCS]
        
        # The number of points:
        P = 1
        
    else:
        # The number of points:
        P = ListShape[0]
    
    #print(f'P = {P}')
    
    if False:
        if P == 1:
            # Since the for loop below is intended to act on a list of points, 
            # convert the point into a list of length 1:
            Pts_PCS = [Pts_PCS]
    
    
    
    # Define S, X, Y and Z:
    S = Origin # the origin
    
    X = Directions[0:3] # the row (x) direction cosine vector
    Y = Directions[3:6] # the column (y) direction cosine vector
    Z = Directions[6:] # the slice (z) direction cosine vector
    
    # The indeces of the largest direction cosines along rows, columns and 
    # slices:
    ind_i = X.index(max(X)) # typically = 0
    ind_j = Y.index(max(Y)) # typically = 1
    ind_k = Z.index(max(Z)) # typically = 2

    # The pixel spacings:
    di = Spacings[0]
    dj = Spacings[1]
    dk = Spacings[2] 
    
    # Simplifying expressions:
    Xx = X[ind_i]
    Xy = X[ind_j]
    Xz = X[ind_k]
    
    Yx = Y[ind_i]
    Yy = Y[ind_j]
    Yz = Y[ind_k]
    
    Zx = Z[ind_i]
    Zy = Z[ind_j]
    Zz = Z[ind_k]
    
    
    # Initialise Pts_ICS:
    Pts_ICS = []
    
    # Loop through each point in Pts_PCS:
    for p in range(P):
        #print(f'Pts_PCS[{p}] = {Pts_PCS[p]}')

        # Get the point:
        point = Pts_PCS[p]
        
        #print(f'\npoint = {point}')
        
        Vx = point[0] - S[0]
        Vy = point[1] - S[1]
        Vz = point[2] - S[2]
        
            
        """ Equations from 17/06 Attempt #2: """
        # Define simplifying expressions:
        a = Xz*Zy/Zz - Xy
        
        b = Yy - Yz*Zy/Zz
        
        c = Vy - Vz*Zy/Zz
        
        d = Vx - (c/b)*Yx - (Zx/Zz)*(Vz - (c/b)*Yz)
        
        e = Xx + (a/b)*Yx - (Zx/Zz)*(Xz + (a/b)*Yz)
        
        # Solve for i, j and k:
        i = (1/di)*d/e
        
        j = (1/dj)*( (a/b)*di*i + c/b )
        
        k = (1/dk)*(1/Zz)*(Vz - Xz*di*i - Yz*dj*j)
        #k = (1/dk)*(1/Zz)*( Vz - (c/b)*Yz - (Xz + (a/b)*Yz)*di*i )
        """ Above two expressions should be the same """

        
        #Pts_ICS.append([di*i, dj*j, dk*k])
        Pts_ICS.append([i, j, k])

     
    """ Note 14/09:
        For the points to be interpretted properly by further functions, even
        a single point should be a list of length 1, so commenting out the 
        lines below. 
    """
    #if P == 1:
    #    # Convert Pts_ICS from a list of length 1 to a single point (as was
    #    # Pts_PCS):
    #    Pts_ICS = Pts_ICS[0]
    
    #if P == 1:
    #    Pts_ICS = [Pts_ICS]
    
    return Pts_ICS





def PCStoICSnested(Pts_PCS, Origin, Directions, Spacings):
    """
    Since PCStoICS cannot currently deal with a nested list of points (e.g. of
    the form in the dictionary ContourData), this will be handled by calling
    PCStoICS in this function.
    """
    
    # Initialise the converted points:
    Pts_ICS = []
    
    for s in range(len(Pts_PCS)):
        Contours_PCS = Pts_PCS[s]
        
        if Contours_PCS:
            #print(f'len(ContoursPCS) = {len(ContoursPCS)}')
            
            Contours_ICS = []
            
            for c in range(len(Contours_PCS)):
                #print(f'ContoursPCS[{c}] =', ContoursPCS[c])
                
                if Contours_PCS[c]:
                    Contours_ICS.append(PCStoICS(Pts_PCS=Contours_PCS[c], 
                                                Origin=Origin, 
                                                Directions=Directions, 
                                                Spacings=Spacings))
                else:
                    Contours_ICS.append([])
                    
            Pts_ICS.append(Contours_ICS)
            
        else:
            Pts_ICS.append([])
            
    return Pts_ICS
    
    
    
    
    
     


    """ PREVIOUS DERIVATIONS: """

    """ 
    Equations from 15/06 Attempt #3 (same results as 17/06 Attempt #2): 
    # Define simplifying expressions:
    a = Yz*Zy/(Yy*Zz) - 1
    
    b = Vz*Zy/(Vy*Zz) - 1
    
    c = 1 - Xz*Zy/(Xy*Zz)
    
    d = 1 - b*Vy*Yz/(a*Vz*Yy)
    
    e = 1 + c*Xy*Yz/(a*Xz*Yy)
    
    f = 1 + c*Xy*Yx/(a*Xx*Yy) - e*Xz*Zx/(Xx*Zz)
    
    g = 1 - b*Vy*Yx/(a*Vx*Yy) - d*Vz*Zx/(Vx*Zz)
    
    # Solve for i, j and k:
    i = g*Vx / (f*Xx*di)

    j = 1/(a*dj*Yy) * (b*Vy + c*Xy*di*i)
    
    k = ( 1/(dk*Zz) ) * (Vz - Xz*di*i - Yz*dj*j)
    """
        
    
        
    """ 
    (INCORRECT) Equations from 17/06 Attempt #1:
    # Define simplifying expressions:
    a = (Xy - Xz*Zy/Zz) / (Yz*Zy/Zz - Yy)
                
    b = (Vz*Zy/Zz - Vy) / (Yz*Zy/Zz - Yy)
    
    c = Vx - b*Yx - (Zx/Zz)*(Vz - b*Yz)
    
    d = Xx + a*Yx - (Zx/Zz)*(Xz + a*Yz)

    # Solve for i, j and k:
    i = (1/di)*c/d
                
    j = (1/dj)*(a*di*i + b)
    
    k = (1/dk)*(1/Zz)*(Vz - b*Yz - (Xz + a*Yz)*di*i)
    #k = (1/dk)*(1/Zz)*(Vz - di*Xz*i - dj*Yz*j)
    # Above two expressions should be the same
    """
    
    
    
    """ Equations derived for 2D but with approximation for k: 
    a = 1 - (Xy*Yx)/(Xx*Yy)
            
    i = (1/di)*( (Vx - Vy*Yx/Yy) / (a*Xx) )
    
    j = (1/dj)*( (Vy - Xy*di*i) / Yy )
    
    k = (1/dk)*Vz # approximation    
    """
    
    
    
    """ Equations adapted from James Petts' equations: 
    # The direction cosines along rows and columns:
    Theta_r = Directions[0:3]
    Theta_c = Directions[3:6]
    
    # The indeces of the largest direction cosines along rows and columns:
    IndMaxTheta_r = Theta_r.index(max(Theta_r)) # typically = 0
    IndMaxTheta_c = Theta_c.index(max(Theta_c)) # typically = 1
    
    # The pixel spacings:
    dx = Spacings[0]
    dy = Spacings[1]
        
    v = [point[0] - S[0], point[1] - S[1], point[2] - S[2]]
         
    i = (v[IndMaxTheta_r] - v[IndMaxTheta_c]*Theta_c[IndMaxTheta_r]/Theta_c[IndMaxTheta_c]) / \
    (
    Theta_r[IndMaxTheta_r]*dx *
    (1 - (Theta_c[IndMaxTheta_r]*Theta_r[IndMaxTheta_c]) / (Theta_r[IndMaxTheta_r]*Theta_c[IndMaxTheta_c]))
    )
    
    # Equation above copied below but without line breaks:
    # i = (v[IndMaxTheta_r] - Theta_c[IndMaxTheta_r]*v[IndMaxTheta_c]/Theta_c[IndMaxTheta_c]) / ( Theta_r[IndMaxTheta_r]*dx * (1 - (Theta_c[IndMaxTheta_r]*Theta_r[IndMaxTheta_c]) / (Theta_r[IndMaxTheta_r]*Theta_c[IndMaxTheta_c])) )
    
    j = (v[IndMaxTheta_c] - Theta_r[IndMaxTheta_c]*i*dx) / (Theta_c[IndMaxTheta_c]*dy)
    
    # Use approximation for j:
    k = v[2]/dz
    
    # Equations above convered to variables used in the function definition above:
    # i = (Vx - Yx*Vy/Yy) / ( Xx*di * (1 - (Yx*Xy) / (Xx*Yy)) )
    
    # j = (Vy - Xy*i*di) / (Yy*dj)
    
    # k = Vz/dk
    
    # The above equations are identical to the derivations I came up with - see
    # "Equations derived for 2D but with approximation for k" above.
    """