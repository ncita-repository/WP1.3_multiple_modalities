# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:32:13 2020

@author: ctorti
"""



def ResampleImage(Image, RefImage, Method, Interpolation):#, RefImageSpacings):
    
    #NewSrcImage = sitk.Resample(image=SrcImage, 
    #                            size=TrgSitkSize,
    #                            transform=sitk.Transform(),
    #                            interpolator=sitk.sitkLinear,
    #                            outputOrigin=TrgSitkOrigin,
    #                            outputSpacing=TrgSitkSpacing,
    #                            outputDirection=TrgSitkDirection,
    #                            defaultPixelValue=0,
    #                            outputPixelType=TrgImage.GetPixelID())
    
    
    # Define which interpolator to use:
    """ 
    Use the linear (or BSpline) interpolator for intensity images, and
    NearestNeighbour for binary images (e.g. segmentation) so that no new
    labels are introduced.
    """
    
    if 'inear' in Interpolation:
        Interpolator = sitk.sitkLinear
    if 'pline' in Interpolation:
        Interpolator = sitk.sitkBSpline
    if ('earest' and 'eighbo' or 'nn' or 'NN') in Interpolation:    
        Interpolator = sitk.sitkNearestNeighbor
            
        
    #print('\nUsing', Interpolation, 'interpolation\n')
    
    Dimension = Image.GetDimension()
    

    
    if 'ResampleImageFilter' in Method:
        Resampler = sitk.ResampleImageFilter()
        
        Resampler.SetReferenceImage(RefImage)
        
        Resampler.SetInterpolator(Interpolator)
        
        if 'Identity' in Method:
            Resampler.SetTransform(sitk.Transform())
        if 'Affine' in Method:
            Resampler.SetTransform(sitk.AffineTransform(Dimension))
            
            
        #print(f'\nImage.GetPixelIDValue() = {Image.GetPixelIDValue()}')
        
        Resampler.SetOutputSpacing(RefImage.GetSpacing())
        Resampler.SetSize(RefImage.GetSize())
        Resampler.SetOutputDirection(RefImage.GetDirection())
        Resampler.SetOutputOrigin(RefImage.GetOrigin())
        #Resampler.SetDefaultPixelValue(Image.GetPixelIDValue())
        Resampler.SetDefaultPixelValue(0)
        
        #print('\n')
        #Resampler.AddCommand(sitk.sitkProgressEvent, lambda: print("\rProgress: {0:03.1f}%...".format(100*Resampler.GetProgress()),end=''))
        #Resampler.AddCommand(sitk.sitkProgressEvent, lambda: sys.stdout.flush())
        
        ResImage = Resampler.Execute(Image)
    
    
    
    if 'Resample2' in Method:
        """
        Use the 2nd method listed in the SimpleITK Transforms and Resampling 
        tutorial (the new Size, Spacing, Origin and Direction will be 
        determined intrinsically from RefImage):
            
        https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/21_Transforms_and_Resampling.ipynb
        
        Also follow the example here:
            
        https://gist.github.com/zivy/79d7ee0490faee1156c1277a78e4a4c4
        """
        
        # Define the transform which maps RefImage to Image with the 
        # translation mapping the image origins to each other:
        Transform = sitk.AffineTransform(Dimension)
        #Transform.SetMatrix(Image.GetDirection())                                      # <-- 9:55 on 30/10
        Transform.SetMatrix(RefImage.GetDirection())                                    # <-- 9:55 on 30/10
        #Transform.SetTranslation(np.array(Image.GetOrigin()) - RefImage.GetOrigin())   # <-- 9:55 on 30/10
        Transform.SetTranslation(np.array(RefImage.GetOrigin()) - Image.GetOrigin())    # <-- 9:55 on 30/10
        
        if 'Centres' in Method:
            # Use the TransformContinuousIndexToPhysicalPoint to compute an 
            # indexed point's physical coordinates as this takes into account  
            # size, spacing and direction cosines. 
            """
            For the vast majority of images the direction cosines are the 
            identity matrix, but when this isn't the case simply multiplying 
            the central index by the spacing will not yield the correct 
            # coordinates resulting in a long debugging session. 
            """
            RefCentre = np.array(RefImage.TransformContinuousIndexToPhysicalPoint(np.array(RefImage.GetSize())/2.0))
            
            # Modify the transformation to align the centres of the original 
            # and reference image instead of their origins:
            CenteringTransform = sitk.TranslationTransform(Dimension)
            
            ImageCentre = np.array(Image.TransformContinuousIndexToPhysicalPoint(np.array(Image.GetSize())/2.0))
            
            CenteringTransform.SetOffset(np.array(Transform.GetInverse().TransformPoint(ImageCentre) - RefCentre))
            
            CenteredTransform = sitk.Transform(Transform)
            
            CenteredTransform.AddTransform(CenteringTransform)

            ResImage = sitk.Resample(Image, RefImage, CenteredTransform, 
                                     Interpolator, Image.GetPixelIDValue())
            
            
        if 'Origins' in Method:
            ResImage = sitk.Resample(Image, RefImage, Transform, Interpolator, 
                                     Image.GetPixelIDValue())
        
        
    
    
    if 'Resample3' in Method:
        """
        Use the 3rd method listed in the SimpleITK Transforms and Resampling 
        tutorial (explicitly define the new Size, Spacing, Origin and 
        Direction):
            
        https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/21_Transforms_and_Resampling.ipynb
        """
        
        NewSize = RefImage.GetSize()
        NewOrigin = RefImage.GetOrigin()
        NewSpacing = RefImage.GetSpacing()
        NewDirection = RefImage.GetDirection()
        
        
        # Define the transform which maps RefImage to Image with the 
        # translation mapping the image origins to each other:
        Transform = sitk.AffineTransform(Dimension)
        #Transform.SetMatrix(Image.GetDirection())                                     # <-- 10:03 on 30/10
        Transform.SetMatrix(RefImage.GetDirection())                                   # <-- 10:03 on 30/10
        #Transform.SetTranslation(np.array(Image.GetOrigin()) - RefImage.GetOrigin())  # <-- 10:03 on 30/10
        Transform.SetTranslation(np.array(RefImage.GetOrigin()) - Image.GetOrigin())   # <-- 10:03 on 30/10
        
        if 'Centres' in Method:
            # Use the TransformContinuousIndexToPhysicalPoint to compute an 
            # indexed point's physical coordinates as this takes into account  
            # size, spacing and direction cosines. 
            """
            For the vast majority of images the direction cosines are the 
            identity matrix, but when this isn't the case simply multiplying 
            the central index by the spacing will not yield the correct 
            # coordinates resulting in a long debugging session. 
            """
            RefCentre = np.array(RefImage.TransformContinuousIndexToPhysicalPoint(np.array(RefImage.GetSize())/2.0))
            
            # Modify the transformation to align the centres of the original 
            # and reference image instead of their origins:
            CenteringTransform = sitk.TranslationTransform(Dimension)
            
            ImageCentre = np.array(Image.TransformContinuousIndexToPhysicalPoint(np.array(Image.GetSize())/2.0))
            
            CenteringTransform.SetOffset(np.array(Transform.GetInverse().TransformPoint(ImageCentre) - RefCentre))
            
            CenteredTransform = sitk.Transform(Transform)
            
            CenteredTransform.AddTransform(CenteringTransform)
        
            ResImage = sitk.Resample(Image, NewSize, CenteredTransform, 
                                     Interpolator, NewOrigin, NewSpacing, 
                                     NewDirection, Image.GetPixelIDValue())  
            
            
        if 'Origins' in Method:
            ResImage = sitk.Resample(Image, NewSize, Transform, 
                                     Interpolator, NewOrigin, NewSpacing, 
                                     NewDirection, Image.GetPixelIDValue())  
            
            
            
            
    if Method == 'NotCorrect':
        """ 
        Use the 3rd method listed in the SimpleITK Transforms and Resampling 
        tutorial (explicitly define the new Size, Spacing, Origin and 
        Direction):
            
        https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/21_Transforms_and_Resampling.ipynb
        """
        
        ExtremePts = GetExtremePoints(RefImage)
        
        Affine = sitk.Euler3DTransform(RefImage.TransformContinuousIndexToPhysicalPoint([size/2 for size in RefImage.GetSize()]),
                                       RefImage.GetAngleX(), 
                                       RefImage.GetAngleY(), 
                                       RefImage.GetAngleZ(), 
                                       RefImage.GetTranslation())
        
        InvAffine = Affine.GetInverse()
        
        # Transformed ExtremePts:
        TxExtremePts = [InvAffine.TransformPoint(pt) for pt in ExtremePts]
        
        min_x = min(TxExtremePts)[0]
        min_y = min(TxExtremePts, key=lambda p: p[1])[1]
        min_z = min(TxExtremePts, key=lambda p: p[2])[2]
        max_x = max(TxExtremePts)[0]
        max_y = max(TxExtremePts, key=lambda p: p[1])[1]
        max_z = max(TxExtremePts, key=lambda p: p[2])[2]
        
        NewSpacing = RefImage.GetSpacing()
        #NewDirection = np.eye(3).flatten()
        NewDirection = RefImage.GetDirection()
        #NewOrigin = [min_x, min_y, min_z]
        NewOrigin = RefImage.GetOrigin()
        
        #NewSize = [math.ceil((max_x - min_x)/NewSpacing[0]), 
        #           math.ceil((max_y - min_y)/NewSpacing[1]),
        #           math.ceil((max_z - min_z)/NewSpacing[2])]
        NewSize = RefImage.GetSize()
        
        ResImage = sitk.Resample(Image, NewSize, Affine, Interpolator, 
                                 NewOrigin, NewSpacing, NewDirection)
    
    
    #sitk.Show(ResImage)

    return ResImage







