# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 09:47:41 2019

@author: Tom O'Reilly
"""

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib import cm
import copy


class interactivePlot(object):
    print('HELLO')
    def __init__(self,fig, ax, X, plotAxis = 2, axisLabels = None, fov = None, cmap = 'gray'):
        
#        set various event handlers
        fig.canvas.mpl_connect('scroll_event', self.onScroll)
        fig.canvas.mpl_connect('button_press_event', self.onClick)
        fig.canvas.mpl_connect('button_release_event', self.onRelease)
        fig.canvas.mpl_connect('motion_notify_event', self.onMotion)
        fig.canvas.mpl_connect('key_press_event', self.keyPress)
        
        self.fig = fig
        self.plotAxis = plotAxis
        self.ax = ax
        self.ax.set_adjustable('box')
        self.mouseClicked = None
        self.cmapRange = np.nanmax(X) - np.nanmin(X)
        self.cmapCenter = np.nanmin(X) + self.cmapRange/2
        self.tempCmapRange = self.cmapRange
        self.tempCmapCenter = self.cmapCenter
        self.X = X
        self.axisLabels = axisLabels
        self.colorMap = cmap
        
        if fov is None:
            self.fov = (1,1,1)
            self.resolution = (1,1,1)
        else:
            self.fov = fov
            self.resolution = np.divide(fov, np.shape(X))

        self.slices = np.shape(X)[plotAxis]
        self.ind = self.slices//2
        
        ax.set_title('Slice %s/%s' % (self.ind,self.slices))
        
        if self.plotAxis == 0:
            imageData = self.X[self.ind,:,:]
        elif self.plotAxis == 1:
            imageData = self.X[:,self.ind,:]
        elif self.plotAxis == 2:
            imageData = self.X[:,:,self.ind]
        else:
            print("invalid axis")
            return -1
        self.im = ax.imshow(imageData, cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        self.updateSlice()

    def keyPress(self, event):
        if event.key == " ": #change orientation on space bar press
            self.plotAxis += 1
            if self.plotAxis > 2:
                self.plotAxis = 0
            self.updateSlice()
                
    def onScroll(self, event):
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.updateSlice()
    
    def onClick(self, event):
#        print('Button pressed: %s at X: %s, Y: %s'%(event.button, event.xdata, event.ydata))
        if event.button == 3:   #Reset image on right click
            self.cmapRange = np.nanmax(self.X) - np.nanmin(self.X)
            self.cmapCenter = np.nanmin(self.X) + self.cmapRange/2
            self.tempCmapRange = self.cmapRange
            self.tempCmapCenter = self.cmapCenter
            self.im.set_clim(self.cmapCenter-self.cmapRange/2, self.cmapCenter+self.cmapRange/2)
            self.im.axes.figure.canvas.draw()
        elif event.button == 2: #change orientation on scroll wheel click
            self.plotAxis += 1
            if self.plotAxis > 2:
                self.plotAxis = 0
            self.updateSlice()
        else:
            self.mouseClicked = event.xdata, event.ydata
        
    def onRelease(self, event):
        self.mouseClicked = None
        self.cmapRange = self.tempCmapRange
        self.cmapCenter = self.tempCmapCenter
        
    def onMotion(self, event):
        if self.mouseClicked == None: return        #if mouse isn't clicked ignore movement

        dx = event.xdata - self.mouseClicked[0]
        dy = event.ydata - self.mouseClicked[1]
        
        normDx = dx/self.mouseClicked[0]
        normDy = dy/self.mouseClicked[1]
        
        self.tempCmapRange = self.cmapRange*(1+normDy)
        self.tempCmapCenter = self.cmapCenter*(1+normDx)
        
        self.im.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
        self.im.axes.figure.canvas.draw()
        
    def updateSlice(self):
        if self.plotAxis == 0:
            if self.ind >= np.size(self.X,0):
                self.ind = np.size(self.X, 0)-1
            imageData = self.X[self.ind,:,:]
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[2])
                self.ax.set_ylabel(self.axisLabels[1])
            self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
            self.ax.set_aspect(self.resolution[1]/self.resolution[2])
            self.slices = np.size(self.X,0)
        elif self.plotAxis == 1:
            if self.ind >= np.size(self.X,1):
                self.ind = np.size(self.X, 1)-1
            imageData = self.X[:,self.ind,:]
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[2])
                self.ax.set_ylabel(self.axisLabels[0])
            self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
            self.ax.set_aspect(self.resolution[0]/self.resolution[2])
            self.slices = np.size(self.X,1)
        else:
            if self.ind >= np.size(self.X,2):
                self.ind = np.size(self.X,2)-1
            imageData = self.X[:,:,self.ind]
            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[1])
                self.ax.set_ylabel(self.axisLabels[0])
            self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
            self.ax.set_aspect(self.resolution[0]/self.resolution[1])
            self.slices = np.size(self.X,2)
        # print("Plot axis: %.0f, Array dims: %.0f , %.0f" %(self.plotAxis, np.size(imageData,0), np.size(imageData,1)))
        self.im.set_data(imageData)
        self.ax.set_title('Slice %s/%s' % (self.ind,self.slices))
        self.im.axes.figure.canvas.draw()

def update(sliderVal):
    global imageData, plotIms
    
    #calulate the mask and apply to the data
    mask = (imageData > sliderVal).astype(float)
    mask[mask ==0] = np.nan
    maskedImage = mask*imageData
    
    #update the plots with the masked data
    if len(plotIms) == 1:
        plotIms[0].set_data(maskedImage)
        plotIms[0].axes.figure.canvas.draw()
    else:
        plotIms[0].set_data(maskedImage[int(np.size(maskedImage,0)/2),:,:])
        plotIms[0].axes.figure.canvas.draw()
        plotIms[1].set_data(maskedImage[:,int(np.size(maskedImage,1)/2),:])
        plotIms[1].axes.figure.canvas.draw()
        plotIms[2].set_data(maskedImage[:,:,int(np.size(maskedImage,2)/2)])
        plotIms[2].axes.figure.canvas.draw()
    
def waitforbuttonpress():
    global closed
    while plt.waitforbuttonpress(0.2) is None:#while no button is pressed
        if closed:
            return False
    return True

def handle_close(evt):
    global closed
    closed = True

def sliderPlot(image, returnMask = True):
    #define globals for sharing data between functions
    global imageData, plotIms, closed
    
    #normalise the image, removes the need to calculate range for slider
    imageData   = np.abs(image)
    # imageData   -= np.min(imageData)
    # imageData   /= np.max(imageData)
    
    numAxis     = np.size(np.shape(imageData))
    if numAxis < 2 or numAxis > 3:
        print("Invalid number of axis, data should be two or three dimensional")
        return -1

    #create list for storing image objects
    plotIms     = []
    
    maskedCmap = copy.copy(cm.get_cmap('gray'))
    maskedCmap.set_bad('red',1.)
    
    if numAxis == 2:
        fig, ax = plt.subplots(1,1)
        plotIms.append(ax.imshow(imageData, cmap = maskedCmap))
    else:
        fig, ax = plt.subplots(1,3)
        plotIms.append(ax[0].imshow(imageData[int(np.size(imageData,0)/2),:,:], cmap = maskedCmap, vmin = np.min(imageData), vmax = np.max(imageData)))
        plotIms.append(ax[1].imshow(imageData[:,int(np.size(imageData,1)/2),:], cmap = maskedCmap, vmin = np.min(imageData), vmax = np.max(imageData)))
        plotIms.append(ax[2].imshow(imageData[:,:,int(np.size(imageData,2)/2)], cmap = maskedCmap, vmin = np.min(imageData), vmax = np.max(imageData)))
    
    #move the plots to make space for the slider and button
    plt.subplots_adjust(left=0.1, bottom=0.25)
    #create axis in which to place the slider
    axSlider = plt.axes([0.1, 0.1, 0.8, 0.03])
    
    #create slider
    mask_slider = Slider(ax=axSlider, label='', valmin=np.min(imageData), valmax=np.max(imageData), valinit=np.min(imageData))
    mask_slider.on_changed(update)
    
    #set up button and wait function to pause script until plot is closed
    closed = False
    
    fig.canvas.mpl_connect('close_event', handle_close)

    
    while True:
        if not waitforbuttonpress():
            break
    if returnMask:
        return (imageData > mask_slider.val).astype(float)
    else:
        return mask_slider.val

    
