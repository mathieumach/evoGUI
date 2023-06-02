import numpy as np
import imageProcessing as imProc
from matplotlib.widgets import Slider, Button
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
        self.X_shape = X.shape
        self.roi = roi
        self.mean_roi = 0
        self.maskroi0 = np.zeros((X[0,:,:].shape))
        self.maskroi1 = np.zeros((X[:,0,:].shape))
        self.maskroi2 = np.zeros((X[:,:,0].shape))
        self.axisLabels = axisLabels
        self.colorMap = cmap
        self.filter_val = 0
        self.original = X
        self.fil = fil
        self.center_ROI = [0,0]
        self.hROI = 5
        self.wROI = 5
        self.mask = np.zeros((self.X[:,:,0].shape))
        
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
        
        #create filter slider
        self.axSlider = fig.add_axes([0.25, 0.05, 0.6, 0.03])
        self.slider = Slider(ax=self.axSlider, label='Filter strength', valmin=0, valmax=1, valinit=0)
        #create ROI box sliders
        ROI_h = np.concatenate([np.linspace(0, 50, 51)])
        self.axSlider_hROI = fig.add_axes([0.25, 0.03, 0.6, 0.03])
        self.slider_hROI = Slider(ax=self.axSlider_hROI, label='Heigth of ROI', valmin=1, valstep = ROI_h, valmax = ROI_h[-1], valinit=5)  
        ROI_w = np.concatenate([np.linspace(0, 50, 51)])
        self.axSlider_wROI = fig.add_axes([0.25, 0.01, 0.6, 0.03])
        self.slider_wROI = Slider(ax=self.axSlider_wROI, label='Width of ROI', valmin=1, valstep = ROI_w, valmax = ROI_w[-1], valinit=5)  
        
        def update(sliderVal):
            kSpace = np.fft.fftshift(np.fft.fftn((np.fft.fftshift(self.original))))
            if self.fil == "Sine bell": 
                kSpace = imProc.sineBellSquaredFilter(kSpace,sliderVal)
            elif self.fil == "Gaussian": 
                kSpace = imProc.gaussianFilter(kSpace,sliderVal)
            self.X = np.abs(np.fft.fftshift(np.fft.ifftn((np.fft.fftshift(kSpace)))))
            self.updateSlice()       
        
        def update_box(sliderVal):
            
            if self.roi == 'mean':
                
                A = int(self.center_ROI[1])
                B = int(self.center_ROI[0])
                # Heigth
                if self.slider_hROI.val == 0:
                    h = 5
                else:
                    h = int(self.slider_hROI.val)
                # Width
                if self.slider_wROI.val == 0:
                    w = 5
                else:
                    w = int(self.slider_wROI.val)   
                
                self.hROI = h
                self.wROI = w
                arrh = np.ones(2*h)
                arrw = np.ones(2*w)
                if self.plotAxis == 0:
                    arr2 = np.zeros((self.X_shape[1], self.X_shape[2])) # This array is only zeros with a white box
                elif self.plotAxis == 1:
                    arr2 = np.zeros((self.X_shape[2], self.X_shape[0])) 
                elif self.plotAxis == 2:
                    arr2 = np.zeros((self.X_shape[1], self.X_shape[0])) 
                    
                # ROI box
                arr2[A-h:A+h, B-w] = arrh
                arr2[A-h:A+h, B+w] = arrh
                arr2[A-h, B-w:B+w] = arrw
                arr2[A+h, B-w:B+w] = arrw
                self.mask = arr2
                
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
                    self.im.set_data(imageData + 100000*self.mask)
                    self.im.axes.figure.canvas.draw()
                    
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
                    self.im.set_data(imageData + 100000*self.mask)
                    self.im.axes.figure.canvas.draw()
                    
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
                    self.im.set_data(imageData + 100000*self.mask)
                    self.im.axes.figure.canvas.draw()
                    
                self.mean_roi = np.mean(imageData[A-h:A+h,B-w:B+w])
                print('Mean signal of ROI is --> ' + str(self.mean_roi))
                self.ax.set_title('Slice ' + str(self.ind) + '/' + str(self.slices) + ' with ROI action --> ' + str(self.roi) + ' Mean of ROI --> ' + str(np.round(self.mean_roi,2)))

        self.slider.on_changed(update)
        self.slider_hROI.on_changed(update_box)
        self.slider_wROI.on_changed(update_box)
        self.updateSlice()
        
        # Saving mask button
        def save(*args):
            if np.sum(self.mask) != 0:                
                date_now = datetime.now().strftime("%Y-%m-%d")
                time_now = datetime.now().strftime("%H-%M-%S")
                
                # Will return the directory of the file
                self.filename = filedialog.askdirectory()
                file_path_label = str(self.filename)
                
                if self.plotAxis == 0:
                    np.save(os.path.join(file_path_label,'ROI mask #coronal #slice '+str(self.ind)+' #date '+date_now+' #time '+time_now+' #ROI value '+str(np.round(self.mean_roi,2))),self.mask)
                elif self.plotAxis == 1:
                    np.save(os.path.join(file_path_label,'ROI mask #sagittal #slice '+str(self.ind)+' #date '+date_now+' #time '+time_now+' #ROI value '+str(np.round(self.mean_roi,2))),self.mask)
                elif self.plotAxis == 2:
                    np.save(os.path.join(file_path_label,'ROI mask #axial #slice '+str(self.ind)+' #date '+date_now+' #time '+time_now+' #ROI value '+str(np.round(self.mean_roi,2))),self.mask)
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
                self.mask = np.zeros((self.X[0,:,:].shape))
            if self.plotAxis == 1:
                self.mask = np.zeros((self.X[:,0,:].shape))
            if self.plotAxis == 2:
                self.mask = np.zeros((self.X[:,:,0].shape))
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
                self.roi = 'center'
                self.ax.set_title('Define center of ROI with left click')                
                self.im.axes.figure.canvas.draw()
            
            elif self.roi == 'center':
                self.roi = 'mean'
                self.im.axes.figure.canvas.draw()
                
            elif self.roi == 'mean':
                self.roi = 'off'
                self.ax.set_title('Slice ' + str(self.ind) + '/' + str(self.slices) + ' with ROI action --> ' + str(self.roi))
                
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
                    self.im.set_data(imageData)
                    
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
                    self.im.set_data(imageData)
                    
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
                    self.im.set_data(imageData)      
                
                self.im.axes.figure.canvas.draw()
            print('Current action is: ' +str(self.roi))
            
        elif event.button == 2: #change orientation on scroll wheel click
            self.plotAxis += 1
            if self.plotAxis > 2:
                self.plotAxis = 0
            self.updateSlice()
        else:
            self.mouseClicked = event.xdata, event.ydata
            if self.roi == 'center':
                self.center_ROI[0] = np.round(self.mouseClicked[0])
                self.center_ROI[1] = np.round(self.mouseClicked[1])
                
                # Define a ROI with fixe shape and draw it
                A = int(self.center_ROI[1])
                B = int(self.center_ROI[0])
                h = 5 # Heigth
                w = 5 # Width     
                arrh = np.ones(2*h)
                arrw = np.ones(2*w)
                if self.plotAxis == 0:
                    arr2 = np.zeros((self.X_shape[1], self.X_shape[2])) # This array is only zeros with a white box
                elif self.plotAxis == 1:
                    arr2 = np.zeros((self.X_shape[2], self.X_shape[0])) 
                elif self.plotAxis == 2:
                    arr2 = np.zeros((self.X_shape[1], self.X_shape[0])) 
                # ROI box
                arr2[A-h:A+h, B-w] = arrh
                arr2[A-h:A+h, B+w] = arrh
                arr2[A-h, B-w:B+w] = arrw
                arr2[A+h, B-w:B+w] = arrw
                self.mask = arr2
                
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
                    self.im.set_data(imageData + 100000*self.mask)
                    
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
                    self.im.set_data(imageData + 100000*self.mask)
    
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
                    self.im.set_data(imageData + 100000*self.mask)
                    
                self.mean_roi = np.mean(imageData[A-5:A+5,B-5:B+5])
                print('Mean signal of ROI is --> ' + str(np.round(self.mean_roi,2)))
                self.ax.set_title('Slice ' + str(self.ind) + '/' + str(self.slices) + ' with ROI action --> ' + str(self.roi) + ' Mean of ROI --> ' + str(np.round(self.mean_roi,2)))
                self.im.axes.figure.canvas.draw() 
        
    def onRelease(self, event):
        self.mouseClicked = None
        self.cmapRange = self.tempCmapRange
        self.cmapCenter = self.tempCmapCenter
        
    def onMotion(self, event):
        if self.mouseClicked == None: return        #if mouse isn't clicked ignore movement
            
        elif self.mouseClicked != None:

            dx = event.xdata - self.mouseClicked[0]
            dy = event.ydata - self.mouseClicked[1]

            normDx = dx/self.mouseClicked[0]
            normDy = dy/self.mouseClicked[1]

            self.tempCmapRange = self.cmapRange*(1+normDy)
            self.tempCmapCenter = self.cmapCenter*(1+normDx)

            self.im.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
            self.im.axes.figure.canvas.draw()
        
    def updateSlice(self):
        
        if self.roi != 'mean':
            if self.plotAxis == 0:
                if self.ind >= np.size(self.X,0):
                    self.ind = np.size(self.X,0)-1
                imageData = self.X[self.ind,:,:]
                if self.axisLabels != None:
                    self.ax.set_xlabel(self.axisLabels[2])
                    self.ax.set_ylabel(self.axisLabels[1])
                self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                self.ax.set_aspect(self.resolution[1]/self.resolution[2])
                self.slices = np.size(self.X,0)
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
        
        elif self.roi == 'mean':
            if self.plotAxis == 0:
                if self.ind >= np.size(self.X,0):
                    self.ind = np.size(self.X,0)-1
                imageData = self.X[self.ind,:,:]
                if self.axisLabels != None:
                    self.ax.set_xlabel(self.axisLabels[2])
                    self.ax.set_ylabel(self.axisLabels[1])
                self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                self.ax.set_aspect(self.resolution[1]/self.resolution[2])
                self.slices = np.size(self.X,0)
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
                      
            self.im.set_data(imageData + 100000*self.mask)
            
            A = int(self.center_ROI[1])
            B = int(self.center_ROI[0])
            h = self.hROI
            w = self.wROI
            self.mean_roi = np.mean(imageData[A-h:A+h,B-w:B+w])
            
            print('Mean signal of ROI is --> ' + str(self.mean_roi))
            self.ax.set_title('Slice ' + str(self.ind) + '/' + str(self.slices) + ' with ROI action --> ' + str(self.roi) + ' Mean of ROI --> ' + str(np.round(self.mean_roi,2)))
            self.im.axes.figure.canvas.draw()



