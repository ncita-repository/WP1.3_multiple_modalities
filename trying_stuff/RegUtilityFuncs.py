# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:18:10 2020

@author: ctorti
"""



""" Define utility functions for plotting of registration results """


""" Source:
https://github.com/InsightSoftwareConsortium/SimpleITKWorkshopImageJ2015/blob/master/6_Similarity_Metrics_and_Optimization.ipynb
"""


# Import packages:
import matplotlib.pyplot as plt


# callback invoked by the interact ipython method for scrolling through the image stacks of
# the two images (moving and fixed)
""" Original: """
def display_images(fixed_image_z, moving_image_z, fixed_npa, moving_npa):
    # create a figure with two subplots and the specified size
    plt.subplots(1,2,figsize=(10,8))

    # draw the fixed image in the first subplot
    plt.subplot(1,2,1)
    plt.imshow(fixed_npa[fixed_image_z,:,:],cmap=plt.cm.Greys_r);
    plt.title('fixed image')
    plt.axis('off')

    # draw the moving image in the second subplot
    plt.subplot(1,2,2)
    plt.imshow(moving_npa[moving_image_z,:,:],cmap=plt.cm.Greys_r);
    plt.title('moving image')
    plt.axis('off')


# callback invoked by the interact ipython method for scrolling through the image stacks of
# the two images (moving and fixed)
""" Slight modification to display_images: Add figure titles as arguments: """
def display_images_v2(fixed_image_z, moving_image_z, fixed_npa, moving_npa, fixed_title, moving_title):
    # Get min, max:
    #fixed_min = np.min(fixed_npa)
    #fixed_max = np.max(fixed_npa)
    #moving_min = np.min(moving_npa)
    #moving_max = np.max(moving_npa)

    # create a figure with two subplots and the specified size
    plt.subplots(1,2,figsize=(10,8))

    # draw the fixed image in the first subplot
    plt.subplot(1,2,1)
    plt.imshow(fixed_npa[fixed_image_z,:,:],cmap=plt.cm.Greys_r);
    #plt.imshow(fixed_npa[fixed_image_z,:,:], cmap=plt.cm.Greys_r, vmin=fixed_min, vmax=fixed_max);
    #plt.imshow(fixed_npa[fixed_image_z,:,:], cmap=plt.cm.Greys_r, norm=True);
    #plt.imshow(fixed_npa[fixed_image_z,:,:], cmap=plt.cm.Greys_r, vmin=0, vmax=255);
    plt.title(fixed_title)
    plt.axis('off')

    # draw the moving image in the second subplot
    plt.subplot(1,2,2)
    plt.imshow(moving_npa[moving_image_z,:,:],cmap=plt.cm.Greys_r);
    #plt.imshow(moving_npa[moving_image_z,:,:], cmap=plt.cm.Greys_r, vmin=moving_min, vmax=moving_max);
    #plt.imshow(moving_npa[moving_image_z,:,:], cmap=plt.cm.Greys_r, norm=True);
    #plt.imshow(moving_npa[moving_image_z,:,:], cmap=plt.cm.Greys_r, vmin=0, vmax=255);
    plt.title(moving_title)
    plt.axis('off')

# callback invoked by the interact ipython method for scrolling through the image stacks of
# the two images (moving and fixed)
""" Slight modification to display_images_v2: Plot registered image as well: """
def display_images_v3(fixed_image_z, moving_image_z, registered_image_z, 
                      fixed_npa, moving_npa, registered_npa, 
                      fixed_title, moving_title, registered_title):

    # create a figure with two subplots and the specified size
    plt.subplots(1,3,figsize=(9,6))

    # draw the fixed image in the first subplot
    plt.subplot(1,3,1)
    plt.imshow(fixed_npa[fixed_image_z,:,:],cmap=plt.cm.Greys_r);
    plt.title(fixed_title)
    plt.axis('off')

    # draw the moving image in the second subplot
    plt.subplot(1,3,2)
    plt.imshow(moving_npa[moving_image_z,:,:],cmap=plt.cm.Greys_r);
    plt.title(moving_title)
    plt.axis('off')
    
    # draw the registered image in the third subplot
    plt.subplot(1,3,3)
    plt.imshow(registered_npa[registered_image_z,:,:],cmap=plt.cm.Greys_r);
    plt.title(registered_title)
    plt.axis('off')
    
# callback invoked by the interact ipython method for scrolling through the image stacks of
# the two images (moving and fixed)
""" Add contours to display_images_v3 """
def display_images_and_contours(fixed_ind, moving_ind, reg_ind, 
                                fixed_npa, moving_npa, reg_npa, 
                                fixed_im, moving_im, reg_im,
                                fixed_contour_pts,
                                fixed_title, moving_title, reg_title):
    
    # Get the number of slices in the images:
    fixed_Nslices = fixed_im.GetSize()[2]
    moving_Nslices = moving_im.GetSize()[2]
    reg_Nslices = reg_im.GetSize()[2]

    # Get the z positions of the origins:
    fixed_origin_z = fixed_im.GetOrigin()[2]
    moving_origin_z = moving_im.GetOrigin()[2]
    reg_origin_z = reg_im.GetOrigin()[2]

    # Get the pixel spacings along the z axis:
    fixed_dz = fixed_im.GetSpacing()[2]
    moving_dz = moving_im.GetSpacing()[2]
    reg_dz = reg_im.GetSpacing()[2]

    # Use the z positions of the origins, the pixel spacings along the z axis,
    # and the indeces to get the z positions of the slices in the stacks:
    fixed_zpos = fixed_origin_z + fixed_ind*fixed_dz
    moving_zpos = moving_origin_z + moving_ind*moving_dz
    reg_zpos = reg_origin_z + reg_ind*reg_dz

    # Get the number of contours for this index:
    fixed_Ncontours = len(fixed_contour_pts[fixed_ind])

    print(f'\nFixed slice  {fixed_ind}/{fixed_Nslices}',
          f'at z = {round(fixed_zpos, 2)} mm')
    print(f'Moving slice {moving_ind}/{moving_Nslices}',
          f'at z = {round(moving_zpos, 2)} mm')
    print(f'Reg slice    {reg_ind}/{reg_Nslices}',
          f'at z = {round(reg_zpos, 2)} mm')
    
    #print(f'\nFixed depth position  {round(fixed_zpos, 2)} mm')
    #print(f'Moving depth position {round(moving_zpos, 2)} mm')
    #print(f'Reg depth position    {round(reg_zpos, 2)} mm')
    
    print(f'\nFixed no. of contours  = {fixed_Ncontours}')
    #print(f'Moving no. of contours = {fixed_Ncontours}')
    #print(f'Reg no. of contours    = {fixed_Ncontours}')


    # create a figure with two subplots and the specified size
    plt.subplots(1, 3, figsize=(14,5))

    # draw the fixed image in the first subplot
    plt.subplot(1, 3, 1, aspect='equal')
    plt.imshow(fixed_npa[fixed_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(fixed_npa[fixed_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(fixed_title)
    plt.axis('off')
    
    # Plot contours for fixed image if Npts is not empty:
    if fixed_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        #for x, y, z in fixed_contour_pts[fixed_ind][0]:
        for x, y, z in fixed_contour_pts[fixed_ind]:
            X.append(x)
            Y.append(y)

        # Flip Y pixel coordinates?:
        if False:
            N_r, N_c = np.shape(fixed_npa[fixed_ind,:,:])

            Y = N_c - np.array(Y)

        plt.plot(X, Y, linewidth=0.5, c='red');

    # draw the moving image in the second subplot
    plt.subplot(1, 3, 2, aspect='equal')
    plt.imshow(moving_npa[moving_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(moving_npa[moving_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(moving_title)
    plt.axis('off')
    
    # draw the registered image in the third subplot
    plt.subplot(1, 3, 3, aspect='equal')
    plt.imshow(reg_npa[reg_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(reg_npa[reg_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(reg_title)
    plt.axis('off')
    
    # Using reg_inds and fixed_contour_pts, plot contours belonging to the 
    # fixed image for the registered image if Npts is not empty:
    if fixed_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        #for x, y, z in fixed_contour_pts[reg_ind][0]:
        for x, y, z in fixed_contour_pts[reg_ind]:
            X.append(x)
            Y.append(y)

        # Flip Y pixel coordinates?:
        if False:
            N_r, N_c = np.shape(reg_npa[reg_ind,:,:])

            Y = N_c - np.array(Y)

        plt.plot(X, Y, linewidth=0.5, c='red');
        
        
""" Similar to display_images_and_contours() but include MovingContourPts: """
def display_images_and_contours_v2(fixed_Zind, moving_Zind, reg_Zind, 
                                   fixed_npa, moving_npa, reg_npa, 
                                   fixed_im, moving_im, reg_im,
                                   fixed_contour_pts, moving_contour_pts,
                                   fixed_title, moving_title, reg_title):
    
    # Get the number of slices in the images:
    fixed_Nslices = fixed_im.GetSize()[2]
    moving_Nslices = moving_im.GetSize()[2]
    reg_Nslices = reg_im.GetSize()[2]

    # Get the z positions of the origins:
    fixed_origin_z = fixed_im.GetOrigin()[2]
    moving_origin_z = moving_im.GetOrigin()[2]
    reg_origin_z = reg_im.GetOrigin()[2]

    # Get the pixel spacings along the z axis:
    fixed_dz = fixed_im.GetSpacing()[2]
    moving_dz = moving_im.GetSpacing()[2]
    reg_dz = reg_im.GetSpacing()[2]

    # Use the z positions of the origins, the pixel spacings along the z axis,
    # and the indeces to get the z positions of the slices in the stacks:
    fixed_zpos = fixed_origin_z + fixed_Zind*fixed_dz
    moving_zpos = moving_origin_z + moving_Zind*moving_dz
    reg_zpos = reg_origin_z + reg_Zind*reg_dz

    # Get the number of contours for this index:
    fixed_Ncontours = len(fixed_contour_pts[fixed_Zind])
    moving_Ncontours = len(moving_contour_pts[moving_Zind])

    print(f'\nFixed slice  {fixed_Zind}/{fixed_Nslices}',
          f'at z = {round(fixed_zpos, 2)} mm')
    print(f'Moving slice {moving_Zind}/{moving_Nslices}',
          f'at z = {round(moving_zpos, 2)} mm')
    print(f'Reg slice    {reg_Zind}/{reg_Nslices}',
          f'at z = {round(reg_zpos, 2)} mm')
    
    #print(f'\nFixed depth position  {round(fixed_zpos, 2)} mm')
    #print(f'Moving depth position {round(moving_zpos, 2)} mm')
    #print(f'Reg depth position    {round(reg_zpos, 2)} mm')
    
    print(f'\nFixed no. of contours  = {fixed_Ncontours}')
    print(f'Moving no. of contours = {moving_Ncontours}')
    #print(f'Reg no. of contours    = {fixed_Ncontours}')


    # create a figure with two subplots and the specified size
    plt.subplots(1, 3, figsize=(14,5))

    # draw the fixed image in the first subplot
    plt.subplot(1, 3, 1, aspect='equal')
    plt.imshow(fixed_npa[fixed_Zind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(fixed_npa[fixed_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(fixed_title)
    plt.axis('off')
    
    # Plot contours for fixed image if fixed_Ncontours is not 0:
    if fixed_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        #for x, y, z in fixed_contour_pts[fixed_ind][0]:
        for x, y, z in fixed_contour_pts[fixed_Zind]:
            X.append(x)
            Y.append(y)

        # Flip Y pixel coordinates?:
        if False:
            N_r, N_c = np.shape(fixed_npa[fixed_Zind,:,:])

            Y = N_c - np.array(Y)

        plt.plot(X, Y, linewidth=0.5, c='red');

    # draw the moving image in the second subplot
    plt.subplot(1, 3, 2, aspect='equal')
    plt.imshow(moving_npa[moving_Zind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(moving_npa[moving_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(moving_title)
    plt.axis('off')
    
    # Plot contours for moving image if moving_Ncontours is not 0:
    if moving_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        for x, y, z in moving_contour_pts[moving_Zind]:
            X.append(x)
            Y.append(y)

        # Flip Y pixel coordinates?:
        if False:
            N_r, N_c = np.shape(moving_npa[moving_Zind,:,:])

            Y = N_c - np.array(Y)

        plt.plot(X, Y, linewidth=0.5, c='red');
    
    # draw the registered image in the third subplot
    plt.subplot(1, 3, 3, aspect='equal')
    plt.imshow(reg_npa[reg_Zind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(reg_npa[reg_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(reg_title)
    plt.axis('off')
    
    # Using reg_inds and fixed_contour_pts, plot contours belonging to the 
    # fixed image for the registered image if fixed_Ncontours is not 0:
    if fixed_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        #for x, y, z in fixed_contour_pts[reg_ind][0]:
        for x, y, z in fixed_contour_pts[reg_Zind]:
            X.append(x)
            Y.append(y)

        # Flip Y pixel coordinates?:
        if False:
            N_r, N_c = np.shape(reg_npa[reg_Zind,:,:])

            Y = N_c - np.array(Y)

        plt.plot(X, Y, linewidth=0.5, c='red');

# callback invoked by the ipython interact method for scrolling and modifying the alpha blending
# of an image stack of two images that occupy the same physical space. 
def display_images_with_alpha(image_z, alpha, fixed, moving):
    img = (1.0 - alpha)*fixed[:,:,image_z] + alpha*moving[:,:,image_z] 
    plt.imshow(sitk.GetArrayFromImage(img),cmap=plt.cm.Greys_r);
    plt.axis('off')
    
    
# callback invoked when the StartEvent happens, sets up our new data
def start_plot():
    global metric_values, multires_iterations
    
    metric_values = []
    multires_iterations = []

# callback invoked when the EndEvent happens, do cleanup of data and figure
def end_plot():
    global metric_values, multires_iterations
    
    del metric_values
    del multires_iterations
    #close figure, we don't want to get a duplicate of the plot latter on
    plt.close()

# callback invoked when the IterationEvent happens, update our data and display new figure    
def plot_values(registration_method):
    global metric_values, multires_iterations
    
    metric_values.append(registration_method.GetMetricValue())                                       
    #clear the output area (wait=True, to reduce flickering), and plot current data
    clear_output(wait=True)
    #plot the similarity metric values
    plt.plot(metric_values, 'r')
    plt.plot(multires_iterations, [metric_values[index] for index in multires_iterations], 'b*')
    plt.xlabel('Iteration Number',fontsize=12)
    plt.ylabel('Metric Value',fontsize=12)
    plt.show()
    
# callback invoked when the sitkMultiResolutionIterationEvent happens, update the index into the 
# metric_values list. 
def update_multires_iterations():
    global metric_values, multires_iterations
    multires_iterations.append(len(metric_values))