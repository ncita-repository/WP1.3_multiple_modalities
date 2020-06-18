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
    
    The conversion uses the Image Plane Attributes from the input image (image
    that the points belong to, as an SimpleITK image)


Input:
    Pts_PCS - Points in Patient Coordinate System 
    
    SitkIm  - SimpleITK image which the points relate to
    
    
    
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
      
Solving for di*i, dj*j and dk*k results in the expressions:
    
    dixi = d/e
        
    djxj = (a/b)*di*i + c/b
    
    dkxk = (1/Zz)*(Vz - Xz*di*i - Yz*dj*j)
    or
    dkxk = (1/Zz)*( Vz - (c/b)*Yz - (Xz + (a/b)*Yz)*di*i )
    
    (either expression for dk*k should be fine)
    
In the case where X = [1, 0, 0]
                  Y = [0, 1, 0]
                  Z = [0, 0, 1]
                  
The expressions simplify to:
    
    a = 0
    
    b = 1
    
    c = Vy
    
    d = Vx
    
    e = 1
    
    dixi = Vx
    
    djxj = Vy
    
    dkxk = Vz
    
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
        

def PCStoICS(Pts_PCS, SitkIm):
    # Import packages:
    #import SimpleITK as sitk
    
    """ Convert points in PCS (Pts_PCS) to points in ICS using Image Plane
    Attributes from the Sitk image they belong to. """
    
    Origin = SitkIm.GetOrigin() 
    # e.g. (-104.92378234863281, -152.4906463623047, -9.22148609161377)
    Direction = SitkIm.GetDirection() 
    # e.g. (0.9972648738498597, 1.8634260454120088e-10, 0.07391056342109259,
    #       -0.022262997957353512, 0.9535561337886633, 0.3003915089278791,
    #       -0.07047787104598309, -0.301215370979001, 0.9509480374756598)
    Spacing = SitkIm.GetSpacing() # e.g. (0.8984375, 0.8984375, 5.0)


    # Define S, X, Y and Z:
    S = Origin # the origin
    
    X = Direction[0:3] # the direction cosine along rows (X)
    Y = Direction[3:6] # the direction cosine along columns (Y)
    Z = Direction[6:9] # the direction cosine along slices (Z)
    
    # The indeces of the largest direction cosines along rows, columns and 
    # slices:
    ind_i = X.index(max(X)) # typically = 0
    ind_j = Y.index(max(Y)) # typically = 1
    ind_k = Z.index(max(Z)) # typically = 2

    # The pixel spacings:
    di = Spacing[0]
    dj = Spacing[1]
    dk = Spacing[2] 
    
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
    
    for p in range(len(Pts_PCS)):
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
        
        # Solve for di*i, dj*j and dk*k:
        dixi = d/e
        
        djxj = (a/b)*di*i + c/b
        
        dkxk = (1/Zz)*(Vz - Xz*di*i - Yz*dj*j)
        #dkxk = (1/Zz)*( Vz - (c/b)*Yz - (Xz + (a/b)*Yz)*di*i )
        """ Above two expressions should be the same """
        
        Pts_ICS.append([dixi, djxj, dkxk])
        
    return Pts_ICS



        """ 
        Equations from 15/06 Attempt #3: 
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
        Equations from 17/06 Attempt #1:
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