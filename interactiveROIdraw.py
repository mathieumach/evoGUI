import numpy as np
import imageProcessing as imProc
from matplotlib.widgets import Slider, Button
import skimage as ski
from datetime import datetime
from tkinter import filedialog
import os

class interactivePlot(object):   
    def __init__(self,fig, ax, X, fil, roi, plotAxis = 2, axisLabels = None, fov = None, cmap = 'gray'):
        
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
        self.roi = roi
        self.mean_roi = 0
        self.maskroi0 = np.zeros((X[0,:,:].shape))
        self.maskroi1 = np.zeros((X[:,0,:].shape))
        self.maskroi2 = np.zeros((X[:,:,0].shape))
        self.mask_ax = 'no'
        self.mask_cor = 'no'
        self.mask_sag = 'no'
        self.axisLabels = axisLabels
        self.colorMap = cmap
        self.filter_val = 0
        self.original = X
        self.fil = fil
        self.corners_plolyg_x = []
        self.corners_plolyg_y = []
        
        if fov is None:
            self.fov = (1,1,1)
            self.resolution = (1,1,1)
        else:
            self.fov = fov
            self.resolution = np.divide(fov, np.shape(X))

        self.slices = np.shape(X)[plotAxis]
        self.ind = self.slices//2
        
        ax.set_title('Slice ' + str(self.ind) + '/' + str(self.slices) + ' with ROI action --> ' + str(self.roi))
        
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
        
        #create slider
        self.axSlider = fig.add_axes([0.15, 0.03, 0.6, 0.03])
        self.slider = Slider(ax=self.axSlider, label='', valmin=0, valmax=1, valinit=0)
        
        def update(sliderVal):
            kSpace = np.fft.fftshift(np.fft.fftn((np.fft.fftshift(self.original))))
            if self.fil == "Sine bell": 
                kSpace = imProc.sineBellSquaredFilter(kSpace,sliderVal)
            elif self.fil == "Gaussian": 
                kSpace = imProc.gaussianFilter(kSpace,sliderVal)
            self.X = np.abs(np.fft.fftshift(np.fft.ifftn((np.fft.fftshift(kSpace)))))
            self.updateSlice()

        self.slider.on_changed(update)
        self.updateSlice()
        
        # Saving mask button
        def save(*args):
            if np.sum(self.maskroi0) != 0 or np.sum(self.maskroi1) != 0 or np.sum(self.maskroi2) != 0:                
                date_now = datetime.now().strftime("%Y-%m-%d")
                time_now = datetime.now().strftime("%H-%M-%S")
                
                # Will return the directory of the file
                self.filename = filedialog.askdirectory()
                file_path_label = str(self.filename)
                
                
                if self.plotAxis == 0:
                    np.save(os.path.join(file_path_label,'ROI mask #coronal #slice '+str(self.ind)+' #date '+date_now+' #time '+time_now+' #ROI value '+str(np.round(self.mean_roi,2))),self.maskroi0)
                elif self.plotAxis == 1:
                    np.save(os.path.join(file_path_label,'ROI mask #sagittal #slice '+str(self.ind)+' #date '+date_now+' #time '+time_now+' #ROI value '+str(np.round(self.mean_roi,2))),self.maskroi1)
                elif self.plotAxis == 2:
                    np.save(os.path.join(file_path_label,'ROI mask #axial #slice '+str(self.ind)+' #date '+date_now+' #time '+time_now+' #ROI value '+str(np.round(self.mean_roi,2))),self.maskroi2)

            else:
                print('Before saving a mask, it must be created')
            
        self.axButton = fig.add_axes([0.83, 0.75, 0.15, 0.1])
        self.save_button = Button(self.axButton, ('Save mask')) 
        self.save_button.on_clicked(save)

    def keyPress(self, event):
        if event.key == " ": #change orientation on space bar press
            self.plotAxis += 1
            if self.plotAxis > 2:
                self.plotAxis = 0
            self.corners_plolyg_x = []
            self.corners_plolyg_y = []
            self.maskroi0 = np.zeros((self.X[0,:,:].shape))
            self.maskroi1 = np.zeros((self.X[:,0,:].shape))
            self.maskroi2 = np.zeros((self.X[:,:,0].shape))
            self.mask_ax = 'no'
            self.mask_cor = 'no'
            self.mask_sag = 'no'
        elif event.key == "up": #change orientation on space bar press
            self.ind = (self.ind + 1) % self.slices
        elif event.key == "down": #change orientation on space bar press
            self.ind = (self.ind - 1) % self.slices
        self.updateSlice()

    def onScroll(self, event):
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.updateSlice()
    
    def onClick(self, event):
        if event.button == 3:   
            if self.roi == 'off':
                self.roi = 'on'
                self.ax.set_title('Define corners of ROI (polygon) with left click')
                self.im.axes.figure.canvas.draw()
                
            elif self.roi == 'on':
                self.roi = 'mean'
                
                if np.sum(self.corners_plolyg_x) != 0:
                
                    if self.plotAxis == 0:
                        self.mask_cor = 'yes'
                        if self.ind >= np.size(self.X,0):
                            self.ind = np.size(self.X, 0)-1
                        imageData = self.X[self.ind,:,:]
                        if self.axisLabels != None:
                            self.ax.set_xlabel(self.axisLabels[2])
                            self.ax.set_ylabel(self.axisLabels[1])
                        self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                        self.ax.set_aspect(self.resolution[1]/self.resolution[2])
                        self.slices = np.size(self.X,0)
                        if self.mask_cor == 'yes':
                            xx, yy = ski.draw.polygon(self.corners_plolyg_x, self.corners_plolyg_y)
                            self.maskroi0[xx, yy] = 100000
                            self.im.set_data(imageData + self.maskroi0)
                            self.mean_roi = np.round(np.mean(imageData * np.divide(self.maskroi0,10000), where = self.maskroi0 > 0),2)
                        elif self.mask_cor == 'no':
                            self.im.set_data(imageData)

                    elif self.plotAxis == 1:
                        self.mask_sag = 'yes'
                        if self.ind >= np.size(self.X,1):
                            self.ind = np.size(self.X,1)-1
                        imageData = self.X[:,self.ind,:]
                        if self.axisLabels != None:
                            self.ax.set_xlabel(self.axisLabels[2])
                            self.ax.set_ylabel(self.axisLabels[0])
                        self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                        self.ax.set_aspect(self.resolution[0]/self.resolution[2])
                        self.slices = np.size(self.X,1)
                        if self.mask_sag == 'yes':
                            xx, yy = ski.draw.polygon(self.corners_plolyg_x, self.corners_plolyg_y)
                            self.maskroi1[xx, yy] = 100000
                            self.im.set_data(imageData + self.maskroi1)
                            self.mean_roi = np.round(np.mean(imageData * np.divide(self.maskroi1,10000), where = self.maskroi1 > 0),2)
                        elif self.mask_sag == 'no':
                            self.im.set_data(imageData)

                    elif self.plotAxis == 2:
                        self.mask_ax = 'yes'
                        if self.ind >= np.size(self.X,2):
                            self.ind = np.size(self.X,2)-1
                        imageData = self.X[:,:,self.ind]
                        if self.axisLabels != None:
                            self.ax.set_xlabel(self.axisLabels[1])
                            self.ax.set_ylabel(self.axisLabels[0])
                        self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                        self.ax.set_aspect(self.resolution[0]/self.resolution[1])
                        self.slices = np.size(self.X,2)
                        if self.mask_ax == 'yes':
                            xx, yy = ski.draw.polygon(self.corners_plolyg_x, self.corners_plolyg_y)
                            self.maskroi2[xx, yy] = 100000
                            self.im.set_data(imageData + self.maskroi2)
                            self.mean_roi = np.round(np.mean(imageData * np.divide(self.maskroi2,10000), where = self.maskroi2 > 0),2)
                        elif self.mask_ax == 'no':
                            self.im.set_data(imageData)

                    print('Mean signal of ROI is --> ' + str(self.mean_roi))
                    self.ax.set_title('Slice ' + str(self.ind) + '/' + str(self.slices) + ' with ROI action --> ' + str(self.roi) + ' Mean of ROI --> ' + str(self.mean_roi))
                    self.im.axes.figure.canvas.draw()
                
            elif self.roi == 'mean':
                self.roi = 'off'
                self.maskroi0 = np.zeros((self.X[0,:,:].shape))
                self.maskroi1 = np.zeros((self.X[:,0,:].shape))
                self.maskroi2 = np.zeros((self.X[:,:,0].shape))
                self.corners_plolyg_x = []
                self.corners_plolyg_y = []
                self.ax.set_title('Slice ' + str(self.ind) + '/' + str(self.slices) + ' with ROI action --> ' + str(self.roi))
                self.updateSlice()
            print('Action is --> ' + str(self.roi))
            
        elif event.button == 2: #change orientation on scroll wheel click
            self.plotAxis += 1
            if self.plotAxis > 2:
                self.plotAxis = 0
            self.updateSlice()
        else:
            self.mouseClicked = event.xdata, event.ydata
            if self.roi == 'on':
                x = int(np.round(self.mouseClicked[1]))
                y = int(np.round(self.mouseClicked[0]))
                
                self.corners_plolyg_x.append(x)
                self.corners_plolyg_y.append(y)
                
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
                    self.maskroi0[x,y] = 100000
                    self.im.set_data(imageData + self.maskroi0)
                    
                elif self.plotAxis == 1:
                    if self.ind >= np.size(self.X,1):
                        self.ind = np.size(self.X,1)-1
                    imageData = self.X[:,self.ind,:]
                    if self.axisLabels != None:
                        self.ax.set_xlabel(self.axisLabels[2])
                        self.ax.set_ylabel(self.axisLabels[0])
                    self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                    self.ax.set_aspect(self.resolution[0]/self.resolution[2])
                    self.slices = np.size(self.X,1)
                    self.maskroi1[x,y] = 100000
                    self.im.set_data(imageData + self.maskroi1)
    
                elif self.plotAxis == 2:
                    if self.ind >= np.size(self.X,2):
                        self.ind = np.size(self.X,2)-1
                    imageData = self.X[:,:,self.ind]
                    if self.axisLabels != None:
                        self.ax.set_xlabel(self.axisLabels[1])
                        self.ax.set_ylabel(self.axisLabels[0])
                    self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                    self.ax.set_aspect(self.resolution[0]/self.resolution[1])
                    self.slices = np.size(self.X,2)
                    self.maskroi2[x,y] = 100000
                    self.im.set_data(imageData + self.maskroi2)
                self.im.axes.figure.canvas.draw()
        
    def onRelease(self, event):
        self.mouseClicked = None
        self.cmapRange = self.tempCmapRange
        self.cmapCenter = self.tempCmapCenter
        
    def onMotion(self, event):
        if self.mouseClicked == None: return #if mouse isn't clicked ignore movement
            
        elif self.mouseClicked != None:

            dx = event.xdata - self.mouseClicked[0]
            dy = event.ydata - self.mouseClicked[1]

            normDx = dx/self.mouseClicked[0]
            normDy = dy/self.mouseClicked[1]

            self.tempCmapRange = self.cmapRange*(1+normDy)
            self.tempCmapCenter = self.cmapCenter*(1+normDx)

            self.im.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
            self.im.axes.figure.canvas.draw()
        
    def images(self):
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
        self.im.set_data(imageData)
        self.ax.set_title('Slice ' + str(self.ind) + '/' + str(self.slices) + ' with ROI action --> ' + str(self.roi))
        self.im.axes.figure.canvas.draw()

    def updateSlice(self):

        if self.roi != 'mean':
            self.images()

        elif self.roi == 'mean':
            if np.sum(self.corners_plolyg_x) != 0:
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
                    if self.mask_cor == 'yes':
                        xx, yy = ski.draw.polygon(self.corners_plolyg_x, self.corners_plolyg_y)
                        self.maskroi0[xx, yy] = 100000
                        self.im.set_data(self.X[self.ind,:,:] + self.maskroi0)
                        self.mean_roi = np.round(np.mean(self.X[self.ind,:,:] * np.divide(self.maskroi0,10000), where = self.maskroi0 > 0),2)
                    elif self.mask_cor == 'no':
                        self.im.set_data(self.X[self.ind,:,:])

                elif self.plotAxis == 1:
                    if self.ind >= np.size(self.X,1):
                        self.ind = np.size(self.X,1)-1
                    imageData = self.X[:,self.ind,:]
                    if self.axisLabels != None:
                        self.ax.set_xlabel(self.axisLabels[2])
                        self.ax.set_ylabel(self.axisLabels[0])
                    self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                    self.ax.set_aspect(self.resolution[0]/self.resolution[2])
                    self.slices = np.size(self.X,1)
                    if self.mask_sag == 'yes':
                        xx, yy = ski.draw.polygon(self.corners_plolyg_x, self.corners_plolyg_y)
                        self.maskroi1[xx, yy] = 100000
                        self.im.set_data(self.X[:,self.ind,:] + self.maskroi1)
                        self.mean_roi = np.round(np.mean(self.X[:,self.ind,:] * np.divide(self.maskroi1,10000), where = self.maskroi1 > 0),2)
                    elif self.mask_sag == 'no':
                        self.im.set_data(self.X[:,self.ind,:])

                elif self.plotAxis == 2:
                    if self.ind >= np.size(self.X,2):
                        self.ind = np.size(self.X,2)-1
                    imageData = self.X[:,:,self.ind]
                    if self.axisLabels != None:
                        self.ax.set_xlabel(self.axisLabels[1])
                        self.ax.set_ylabel(self.axisLabels[0])
                    self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                    self.ax.set_aspect(self.resolution[0]/self.resolution[1])
                    self.slices = np.size(self.X,2)
                    if self.mask_ax == 'yes':
                        xx, yy = ski.draw.polygon(self.corners_plolyg_x, self.corners_plolyg_y)
                        self.maskroi2[xx, yy] = 100000
                        self.im.set_data(self.X[:,:,self.ind] + self.maskroi2)
                        self.mean_roi = np.round(np.mean(self.X[:,:,self.ind] * np.divide(self.maskroi2,10000), where = self.maskroi2 > 0),2)
                    elif self.mask_ax == 'no':
                        self.im.set_data(self.X[:,:,self.ind])
                print('Mean signal of ROI is --> ' + str(self.mean_roi))
                self.ax.set_title('Slice ' + str(self.ind) + '/' + str(self.slices) + ' with ROI action --> ' + str(self.roi) + ' Mean of ROI --> ' + str(self.mean_roi))
                self.im.axes.figure.canvas.draw()
            else:
                self.images()




