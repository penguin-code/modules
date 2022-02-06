#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 13:14:38 2021

@author: prowe

Notes:
    
From the documentation
    matplotlib.pyplot.get_cmap(name): Gets a colormap instance
    If name is a matplotlib.colors.Colormap instance, it will be returned.
"""
    

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import imageio

'''

color2index = {
    (255, 255, 255) : 0,
    (0,     0, 255) : 1,
    (0,   255, 255) : 2,
    (0,   255,   0) : 3,
    (255, 255,   0) : 4,
    (255,   0,   0) : 5
}

def rgb2mask(img):

    assert len(img.shape) == 3
    height, width, ch = img.shape
    assert ch == 3

    W = np.power(256, [[0],[1],[2]])

    img_id = img.dot(W).squeeze(-1) 
    values = np.unique(img_id)

    mask = np.zeros(img_id.shape)

    for i, c in enumerate(values):
        try:
            mask[img_id==c] = color2index[tuple(img[img_id==c][0])] 
        except:
            pass
    return mask


# Load matrix of indices for Willamette colormap image
X = np.loadtxt('/Users/prowe/GitHub/PENGUIN/arctic_ice_computer_science/WillametteMapInd')
mmap = np.loadtxt('/Users/prowe/GitHub/PENGUIN/arctic_ice_computer_science/WillametteMap')


# .. Set the colormap
j = 1
thismap = list(zip(mmap[j:,0], mmap[j:,1], mmap[j:,2]))
TCmap = matplotlib.colors.ListedColormap(thismap)
cm = plt.get_cmap(TCmap)

# .. Alternate way
mmap_list = []
for i in range(len(mmap)):
    mmap_list.append(list(mmap[i]))



mdir = '/Users/prowe/GitHub/PENGUIN/arctic_ice_computer_science/'

willamette_color_image = imageio.imread(mdir+'Willamette_colors2.png')
img = willamette_color_image[:,:,:3]

mask = rgb2mask(img)

mmap_list = [
    [1, 255, 255],
    [0,     0, 255],
    [0,   255, 255],
    [0,   255,   0],
    [1, 255,   0],
    [1,   0,   0]]

my_cmap = LinearSegmentedColormap.from_list('some_cmap_name', mmap_list)

# Plot the figure for the custom colormap
plt.figure(1); plt.clf()
plt.imshow(mask, cmap = my_cmap)
plt.axis('off')    # Remove axis ticks and labels
plt.axis('image')  # Set aspect ratio to obtain square pixels (mb, check matplotlib document)




'''
# Load matrix of indices for Willamette colormap image
X = np.loadtxt('/Users/prowe/GitHub/PENGUIN/arctic_ice_computer_science/WillametteMapInd')
mmap = np.loadtxt('/Users/prowe/GitHub/PENGUIN/arctic_ice_computer_science/WillametteMap')


# .. Set the colormap
j = 1
thismap = list(zip(mmap[j:,0], mmap[j:,1], mmap[j:,2]))
TCmap = matplotlib.colors.ListedColormap(thismap)
cm = plt.get_cmap(TCmap)

# .. Alternate way
mmap_list = []
for i in range(len(mmap)):
    mmap_list.append(list(mmap[i]))

my_cmap = LinearSegmentedColormap.from_list('some_cmap_name', mmap_list)

# Plot the figure for the custom colormap
plt.figure(1); plt.clf()
plt.imshow(X, interpolation='none', cmap = my_cmap)
plt.axis('off')    # Remove axis ticks and labels
plt.axis('image')  # Set aspect ratio to obtain square pixels (mb, check matplotlib document)


a = np.array([[0,1,2,3,4],[5,6,7,8,9], [10,11,12,13,14], [15,16,17,18,19],
              [20,21,22,23,24], [25,26,27,28,29], [30,31,0,0,0]])
# Plot the figure for the custom colormap
plt.figure(2); plt.clf()
plt.imshow(a, interpolation='none', cmap = my_cmap)
plt.axis('off')    # Remove axis ticks and labels
plt.axis('image')  # Set aspect ratio to obtain square pixels (mb, check matplotlib document)



mmap_list = [[0,0,0],   # black, 0
             [0,0,1],   # blue, 1
             [0,1,0],   # green, 2
             [0,1,1],   # cyan, 3
             [1,0,0],   # 4, red
             [1,0,1],   # 5, magenta
             [1,1,0],   # 6, yellow
             [.9,.9,.9],   # 7, light grey
             ]
my_cmap = LinearSegmentedColormap.from_list('some_cmap_name', mmap_list)


a = np.array([[0,1],[2,3],[4,5],[6,7]])
# Plot the figure for the custom colormap
plt.figure(2); plt.clf()
plt.imshow(a, cmap = my_cmap)
plt.axis('off')    # Remove axis ticks and labels
plt.axis('image')  # Set aspect ratio to obtain square pixels (mb, check matplotlib document)

