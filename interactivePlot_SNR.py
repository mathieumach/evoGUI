import numpy as np
import imageProcessing as imProc
from matplotlib.widgets import Slider, Button
from datetime import datetime
from tkinter import filedialog
import os

class interactivePlot(object):   
    def __init__(self,fig, ax, X, fil, Data_mat, func, plotAxis = 2, axisLabels = None, fov = None, cmap = 'gray'):
        
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
        self.filter_val = 0
        self.original = X
        self.fil = fil
        self.Data_mat = Data_mat
        self.box = 'off'
        self.std_box_defined = 'no'
        self.mean_box_defined = 'no'
        self.noise_box_center = [0,0]
        self.mean_box_center = [0,0]
        self.snr = 0
        self.mean = 0
        self.std = 0
        self.func = func
        self.length = 5
        self.mask = np.zeros((self.X[:,:,0].shape))
        
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
        
        #create sliders
        self.axSlider = fig.add_axes([0.15, 0.03, 0.6, 0.03])
        self.slider = Slider(ax=self.axSlider, label='Filter', valmin=0, valmax=1, valinit=0)
        
        def update(sliderVal):
            kSpace = np.fft.fftshift(np.fft.fftn((np.fft.fftshift(self.original))))
            if self.fil == "Sine bell": 
                kSpace = imProc.sineBellSquaredFilter(kSpace,sliderVal)
            elif self.fil == "Gaussian": 
                kSpace = imProc.gaussianFilter(kSpace,sliderVal)
            self.X = np.abs(np.fft.fftshift(np.fft.ifftn((np.fft.fftshift(kSpace)))))
            self.updateSlice()

        self.slider.on_changed(update)
            
        box_size = np.concatenate([np.linspace(0, 30, 31)])
        self.axSlider_box = fig.add_axes([0.15, 0.01, 0.6, 0.03])
        self.slider_box_size = Slider(ax=self.axSlider_box, label='box size', valmin=1, valstep = box_size, valmax = box_size[-1], valinit=5)  
        
        def update(val):    
            
            if self.box != 'okay':
                print('the std and mean boxes must be defined before computing the SNR!')
            elif self.box == 'okay':
                
                if self.std_box_defined == 'yes' and self.mean_box_defined == 'yes':
                    
                    if self.plotAxis == 0:
                        imageData = self.X[self.ind,:,:]
                        arr2 = np.zeros((self.Data_mat[1], self.Data_mat[2])) # This array is only zeros with 2 white boxes
                    elif self.plotAxis == 1:
                        imageData = self.X[:,self.ind,:]
                        arr2 = np.zeros((self.Data_mat[2], self.Data_mat[0])) 
                    elif self.plotAxis == 2:
                        imageData = self.X[:,:,self.ind]
                        arr2 = np.zeros((self.Data_mat[1], self.Data_mat[0]))
                
                    # box size 
                    if self.slider_box_size.val == 0:
                        m = 5
                    else:
                        m = self.slider_box_size.val
                    m = int(m)
                    M = 2*m
                    self.length = m
                    
                    arrstd = np.ones(M)
                    arrm = np.ones(M)
                    # std box
                    A = int(self.noise_box_center[1])
                    B = int(self.noise_box_center[0])
                    arr2[A-m:A+m, B-m] = arrstd
                    arr2[A-m:A+m, B+m] = arrstd
                    arr2[A-m, B-m:B+m] = arrstd
                    arr2[A+m, B-m:B+m] = arrstd
                    # mean box
                    a = int(self.mean_box_center[1])
                    b = int(self.mean_box_center[0])
                    arr2[a-m:a+m, b-m] = arrm
                    arr2[a-m:a+m, b+m] = arrm
                    arr2[a-m, b-m:b+m] = arrm
                    arr2[a+m, b-m:b+m] = arrm

                    self.arr2 = arr2
                    self.snr, self.mean, self.std = self.func(np.rot90(self.X[:,:,self.ind]),A-m,A+m,B-m,B+m, a-m,a+m,b-m,b+m)
                    self.fig.suptitle("SNR of image:" + str(np.round(self.snr,2))+", mean:"+str(np.round(self.mean,2))+", std:"+str(np.round(self.std,2)), fontsize=14)
                    self.ax.set_title('Slice %s/%s' % (self.ind,self.slices))       
                    
                    if self.axisLabels != None:
                        self.ax.set_xlabel(self.axisLabels[1])
                        self.ax.set_ylabel(self.axisLabels[0])
                    self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
                    self.im.set_data(imageData + 100000*self.arr2)
                    self.im.axes.figure.canvas.draw()
                    
                else:
                    print('the std and mean boxes must be defined before computing the SNR!')              
        
        self.slider_box_size.on_changed(update)        
        self.updateSlice()
        
        # Saving mask button
        def save(*args):
            if np.sum(self.arr2) != 0:                
                date_now = datetime.now().strftime("%Y-%m-%d")
                time_now = datetime.now().strftime("%H-%M-%S")

                # Will return the directory of the file
                self.filename = filedialog.askdirectory()
                file_path_label = str(self.filename)

                if self.plotAxis == 0:
                    np.save(os.path.join(file_path_label,'SNR mask #coronal #slice '+str(self.ind)+' #date '+date_now+' #time '+time_now +' #mean '+str(np.round(self.mean,2))+' #std ' +str(np.round(self.std,2))+' #SNR ' +str(np.round(self.snr,2))),self.arr2)
                elif self.plotAxis == 1:
                    np.save(os.path.join(file_path_label,'SNR mask #sagittal #slice '+str(self.ind)+' #date '+date_now+' #time '+time_now+' #mean '+str(np.round(self.mean,2))+' #std ' +str(np.round(self.std,2))+' #SNR ' +str(np.round(self.snr,2))),self.arr2)
                elif self.plotAxis == 2:
                    np.save(os.path.join(file_path_label,'SNR mask #axial #slice '+str(self.ind)+' #date '+date_now+' #time '+time_now+' #mean '+str(np.round(self.mean,2))+' #std ' +str(np.round(self.std,2))+' #SNR ' +str(np.round(self.snr,2))),self.arr2)
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
        elif event.key == "up": #change orientation on space bar press
            self.ind = (self.ind + 1) % self.slices
        elif event.key == "down": #change orientation on space bar press
            self.ind = (self.ind - 1) % self.slices
        elif event.key == "right":
            if self.filter_val < 1:
                self.filter_val += 1/100
                kSpace = np.fft.fftshift(np.fft.fftn((np.fft.fftshift(self.original))))
                
                if self.fil == "Sine bell": 
                    kSpace = imProc.sineBellSquaredFilter(kSpace,self.filter_val)
                elif self.fil == "Gaussian": 
                    kSpace = imProc.gaussianFilter(kSpace,self.filter_val)
                    
                self.X = np.abs(np.fft.fftshift(np.fft.ifftn((np.fft.fftshift(kSpace)))))
        elif event.key == "left":
            if self.filter_val > 0.01:
                self.filter_val -= 1/100
                kSpace = np.fft.fftshift(np.fft.fftn((np.fft.fftshift(self.original))))
                
                if self.fil == "Sine bell": 
                    kSpace = imProc.sineBellSquaredFilter(kSpace,self.filter_val)
                elif self.fil == "Gaussian": 
                    kSpace = imProc.gaussianFilter(kSpace,self.filter_val)
                
                self.X = np.abs(np.fft.fftshift(np.fft.ifftn((np.fft.fftshift(kSpace)))))    
        self.updateSlice()


    def onScroll(self, event):
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.updateSlice()
    
    def onClick(self, event):
        if event.button == 3:   #Reset image on right click
            if self.box == 'off':
                self.box = 'std'
                self.ax.set_title('Select center of noise box (for standard deviation)')
                self.im.axes.figure.canvas.draw()
                
            elif self.box == 'std':
                self.box = 'mean'
                self.ax.set_title('Select center of mean box')
                self.im.axes.figure.canvas.draw()
                
            elif self.box == 'mean':
                self.box = 'okay'
                self.updateSlice()
                self.im.axes.figure.canvas.draw()
                
            elif self.box == 'okay':
                self.box = 'off'
                self.std_box_defined = 'no'
                self.mean_box_defined = 'no'
                self.updateSlice()
                
            print('Action is ' + str(self.box))
            
        else:
            self.mouseClicked = event.xdata, event.ydata
            if self.box == 'std':
                self.std_box_defined = 'yes'
                self.noise_box_center[0] = np.round(self.mouseClicked[0])
                self.noise_box_center[1] = np.round(self.mouseClicked[1])
            elif self.box == 'mean':
                self.mean_box_defined = 'yes'
                self.mean_box_center[0] = np.round(self.mouseClicked[0])
                self.mean_box_center[1] = np.round(self.mouseClicked[1])
                
                m = 5
                A = int(self.noise_box_center[1])
                B = int(self.noise_box_center[0])
                a = int(self.mean_box_center[1])
                b = int(self.mean_box_center[0])
                
                self.snr, self.mean, self.std = self.func(np.rot90(self.X[:,:,self.ind]),A-m,A+m,B-m,B+m, a-m,a+m,b-m,b+m)
                self.fig.suptitle("SNR of image:" + str(np.round(self.snr,2))+", mean:"+str(np.round(self.mean,2))+", std:"+str(np.round(self.std,2)), fontsize=14)
            self.updateSlice()
        
        
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
            arr2 = np.zeros((self.Data_mat[1], self.Data_mat[2])) # This array is only zeros with 2 white boxes
        elif self.plotAxis == 1:
            arr2 = np.zeros((self.Data_mat[2], self.Data_mat[0])) 
        elif self.plotAxis == 2:
            arr2 = np.zeros((self.Data_mat[1], self.Data_mat[0])) 
        
        if self.box == 'off': # Will not draw the boxes

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
            self.ax.set_title('Slice %s/%s' % (self.ind,self.slices))
            self.fig.suptitle('')
            self.im.axes.figure.canvas.draw()    
      
        elif self.box != 'okay':             
            if self.std_box_defined == 'no' and self.mean_box_defined == 'no':    
            
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
                self.ax.set_title('Slice %s/%s' % (self.ind,self.slices))
                self.im.axes.figure.canvas.draw()   

            elif self.std_box_defined == 'yes' and self.mean_box_defined == 'no':
            
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
                
                # Creation of the boxes with center defined by user but start size is 10x10
                M = 10; m = 5
                arrstd = np.ones(M)               
                # std box
                A = int(self.noise_box_center[1]); B = int(self.noise_box_center[0])
                arr2[A-m:A+m, B-m] = arrstd; arr2[A-m:A+m, B+m] = arrstd
                arr2[A-m, B-m:B+m] = arrstd; arr2[A+m, B-m:B+m] = arrstd
                self.arr2 = arr2
                self.im.set_data(imageData + 100000*self.arr2)                
                self.im.axes.figure.canvas.draw()    
                
            elif self.std_box_defined == 'yes' and self.mean_box_defined == 'yes': 
            
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
                           
                # Creation of the boxes with center defined by user but start size is 10x10
                M = 10; m = 5
                arrstd = np.ones(M); arrm = np.ones(M)
                arr2 = np.zeros((self.Data_mat[1], self.Data_mat[0])) # This array is only zeros with 2 white boxes
                # std box
                A = int(self.noise_box_center[1]); B = int(self.noise_box_center[0])
                arr2[A-m:A+m, B-m] = arrstd; arr2[A-m:A+m, B+m] = arrstd
                arr2[A-m, B-m:B+m] = arrstd; arr2[A+m, B-m:B+m] = arrstd
                # mean box
                a = int(self.mean_box_center[1]); b = int(self.mean_box_center[0])
                arr2[a-m:a+m, b-m] = arrm; arr2[a-m:a+m, b+m] = arrm
                arr2[a-m, b-m:b+m] = arrm; arr2[a+m, b-m:b+m] = arrm
                self.arr2 = arr2
                self.im.set_data(imageData + 100000*self.arr2)
                self.ax.set_title('Slice %s/%s' % (self.ind,self.slices))
                self.im.axes.figure.canvas.draw()
              
        elif self.box == 'okay': # Will draw both boxes
            
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
                
                
            m = self.slider_box_size.val
            m = int(m)
            M = 2*m
            self.length = m

            arrstd = np.ones(M)
            arrm = np.ones(M)
            # std box
            A = int(self.noise_box_center[1])
            B = int(self.noise_box_center[0])
            arr2[A-m:A+m, B-m] = arrstd
            arr2[A-m:A+m, B+m] = arrstd
            arr2[A-m, B-m:B+m] = arrstd
            arr2[A+m, B-m:B+m] = arrstd
            # mean box
            a = int(self.mean_box_center[1])
            b = int(self.mean_box_center[0])
            arr2[a-m:a+m, b-m] = arrm
            arr2[a-m:a+m, b+m] = arrm
            arr2[a-m, b-m:b+m] = arrm
            arr2[a+m, b-m:b+m] = arrm

            self.arr2 = arr2
            self.snr, self.mean, self.std = self.func(np.rot90(imageData),A-m,A+m,B-m,B+m, a-m,a+m,b-m,b+m)
            self.fig.suptitle("SNR of image:" + str(np.round(self.snr,2))+", mean:"+str(np.round(self.mean,2))+", std:"+str(np.round(self.std,2)), fontsize=14)      

            if self.axisLabels != None:
                self.ax.set_xlabel(self.axisLabels[1])
                self.ax.set_ylabel(self.axisLabels[0])
            self.im.set_extent((-0.5,np.size(imageData,1) - 0.5,np.size(imageData,0) - 0.5,-0.5))
            self.im.set_data(imageData + 100000*self.arr2)
            self.ax.set_title('SNR will is defined based on the boxes, Slice %s/%s' % (self.ind,self.slices))
            self.im.axes.figure.canvas.draw()

