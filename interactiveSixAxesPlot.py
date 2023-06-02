import numpy as np
import imageProcessing as imProc
from matplotlib.widgets import Slider, Button, RadioButtons

class interactivePlot(object):   
    def __init__(self,fig, ax, X, K, fil, plotAxis = 2, axisLabels = None, fov = None, cmap = 'gray'):
        
#        set various event handlers
        fig.canvas.mpl_connect('scroll_event', self.onScroll)
        fig.canvas.mpl_connect('button_press_event', self.onClick)
        fig.canvas.mpl_connect('button_release_event', self.onRelease)
        fig.canvas.mpl_connect('motion_notify_event', self.onMotion)
        fig.canvas.mpl_connect('key_press_event', self.keyPress)
        
        self.fig = fig
        self.plotAxis = plotAxis
        self.ax = ax
        self.mouseClicked = None
        self.cmapRange = np.nanmax(X) - np.nanmin(X)
        self.cmapCenter = np.nanmin(X) + self.cmapRange/2
        self.tempCmapRange = self.cmapRange
        self.tempCmapCenter = self.cmapCenter
        self.X = X
        self.K = K
        self.axisLabels = axisLabels
        self.colorMap = cmap
        self.originalim = X
        self.originalk = K
        self.originallog = np.log(K)
        self.fil = fil
        self.index = 0
        
        if fov is None:
            self.fov = (1,1,1)
            self.resolution = (1,1,1)
        else:
            self.fov = fov
            self.resolution = np.divide(fov, np.shape(X))
            
        self.shape = np.shape(X)
        self.slices = np.shape(X)[plotAxis]
        self.ind = self.slices//2

        self.x_ind = self.ind
        self.y_ind = self.ind
        self.z_ind = self.ind
        
        imageDataX = self.X[self.ind,:,:]
        imageDataY = self.X[:,self.ind,:]
        imageDataZ = self.X[:,:,self.ind]
        imageK_DataX = self.K[self.ind,:,:]
        imageK_DataY = self.K[:,self.ind,:]
        imageK_DataZ = self.K[:,:,self.ind]
        
        self.imx = ax[0,0].imshow(imageDataX, cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        self.ax[0,0].set_title('Slice cor %s/%s' % (self.x_ind,self.shape[0]))
        
        self.imy = ax[0,1].imshow(imageDataY, cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        self.ax[0,1].set_title('Slice sag %s/%s' % (self.y_ind,self.shape[1]))
        
        self.imz = ax[0,2].imshow(imageDataZ, cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        self.ax[0,2].set_title('Slice ax %s/%s' % (self.z_ind,self.shape[2]))
        
        self.imxK = ax[1,0].imshow(imageK_DataX, cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        self.ax[1,0].set_title('Slice cor %s/%s' % (self.x_ind,self.shape[0]))
        
        self.imyK = ax[1,1].imshow(imageK_DataY, cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        self.ax[1,1].set_title('Slice sag %s/%s' % (self.y_ind,self.shape[1]))
        
        self.imzK = ax[1,2].imshow(imageK_DataZ, cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)
        self.ax[1,2].set_title('Slice ax %s/%s' % (self.z_ind,self.shape[2]))
        
        # Global title
        self.fig.suptitle("Action will be used on coronal axe", fontsize=14)
        
        #create slider
        self.axSlider = fig.add_axes([0.15, 0.03, 0.6, 0.03])
        self.slider = Slider(ax=self.axSlider, label='', valmin=0, valmax=1, valinit=0)
        
        def update(sliderVal):
            kSpace = np.fft.fftshift(np.fft.fftn((np.fft.fftshift(self.originalim))))
            if self.fil == "Sine bell": 
                kSpace = imProc.sineBellSquaredFilter(kSpace,sliderVal)
            elif self.fil == "Gaussian": 
                kSpace = imProc.gaussianFilter(kSpace,sliderVal)
            self.X = np.abs(np.fft.fftshift(np.fft.ifftn((np.fft.fftshift(kSpace)))))

            if self.log == 'raw Kspace data':
                if self.fil == "Sine bell": 
                    kSpace = imProc.sineBellSquaredFilter(self.originalk,sliderVal)
                elif self.fil == "Gaussian": 
                    kSpace = imProc.gaussianFilter(self.originalk,sliderVal)
                self.K = kSpace
            elif self.log == 'log(Kspace data)':
                if self.fil == "Sine bell": 
                    kSpace = imProc.sineBellSquaredFilter(self.originallog,sliderVal)
                elif self.fil == "Gaussian": 
                    kSpace = imProc.gaussianFilter(self.originallog,sliderVal)
                self.K = kSpace
                
            self.updateAllSlice()
        self.slider.on_changed(update)
        
        #create radiobutton
        self.axRadio = fig.add_axes([0.85, 0.8, 0.15, 0.1])
        self.radio = RadioButtons(self.axRadio, ('no cross line', 'cross line'))
        self.plane = 'no cross line'
        
        def planefunc(label):
            self.plane = label
            if self.plane == 'cross line':
                
                imageDataX = self.X[self.x_ind,:,:]
                s = self.X[self.x_ind,:,:].shape
                arr1 = np.ones(s[0])
                arr2 = np.ones(s[1])
                arr3 = np.zeros((s)) # This array is only zeros with a vertical and horizontal white lines
                arr3[int(self.shape[1]/2), :] = arr2
                arr3[:, int(self.shape[2]/2)] = arr1
                if self.x_ind >= np.size(self.X,0):
                    self.x_ind = np.size(self.X, 0)-1
                self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
                self.slices = np.size(self.X,0)
                self.imx.set_data(imageDataX + 10000*arr3)
                self.ax[0,0].set_title('Slice cor %s/%s' % (self.x_ind,self.slices))
                self.imx.axes.figure.canvas.draw()
                
                imageDataY = self.X[:,self.y_ind,:]
                s = self.X[:,self.y_ind,:].shape
                arr1 = np.ones(s[0])
                arr2 = np.ones(s[1])
                arr3 = np.zeros((s)) # This array is only zeros with a vertical and horizontal white lines
                arr3[int(self.shape[0]/2), :] = arr2
                arr3[:, int(self.shape[2]/2)] = arr1
                if self.y_ind >= np.size(self.X,0):
                    self.y_ind = np.size(self.X, 0)-1
                self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
                self.slices = np.size(self.X,0)
                self.imy.set_data(imageDataY + 10000*arr3)
                self.ax[0,1].set_title('Slice sag %s/%s' % (self.y_ind,self.slices))
                self.imy.axes.figure.canvas.draw()
                
                imageDataZ = self.X[:,:,self.z_ind]
                s = self.X[:,:,self.z_ind].shape
                arr1 = np.ones(s[0])
                arr2 = np.ones(s[1])
                arr3 = np.zeros((s)) # This array is only zeros with a vertical and horizontal white lines
                arr3[int(self.shape[0]/2), :] = arr2
                arr3[:, int(self.shape[2]/2)] = arr1
                if self.z_ind >= np.size(self.X,0):
                    self.z_ind = np.size(self.X, 0)-1
                self.imz.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
                self.slices = np.size(self.X,0)
                self.imz.set_data(imageDataZ + 10000*arr3)
                self.ax[0,2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
                self.imz.axes.figure.canvas.draw()    
                
            elif self.plane == 'no cross line':
                imageDataX = self.X[self.x_ind,:,:]

                if self.x_ind >= np.size(self.X,0):
                    self.x_ind = np.size(self.X, 0)-1
                self.imxK.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
                self.slices = np.size(self.X,0)
                self.imxK.set_data(imageDataX)
                self.ax[0,0].set_title('Slice cor %s/%s' % (self.x_ind,self.slices))
                self.imxK.axes.figure.canvas.draw()
                
                imageDataY = self.X[:,self.y_ind,:]
                if self.y_ind >= np.size(self.X,0):
                    self.y_ind = np.size(self.X, 0)-1
                self.imyK.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
                self.slices = np.size(self.X,0)
                self.imyK.set_data(imageDataY)
                self.ax[0,1].set_title('Slice sag %s/%s' % (self.y_ind,self.slices))
                self.imyK.axes.figure.canvas.draw()
                
                imageDataZ = self.X[:,:,self.z_ind]
                if self.z_ind >= np.size(self.X,0):
                    self.z_ind = np.size(self.X, 0)-1
                self.imzK.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
                self.slices = np.size(self.X,0)
                self.imzK.set_data(imageDataZ)
                self.ax[0,2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
                self.imzK.axes.figure.canvas.draw()   
            
        self.radio.on_clicked(planefunc)
        
        self.updateSlice()
        
        #create radiobutton
        self.axRadio = fig.add_axes([0.85, 0.1, 0.15, 0.1])
        self.radiolog = RadioButtons(self.axRadio, ('raw Kspace data', 'log(Kspace data)'))
        self.log = 'raw Kspace data'
        
        def data_transform(label):
            self.log = label
            if label == 'log(Kspace data)':
                self.K = np.log(self.K)
            elif label == 'raw Kspace data':
                self.K = np.exp(self.K)
            self.updateAllSlice()
            
        self.radiolog.on_clicked(data_transform) 
            
    def keyPress(self, event):
        if event.key == "up": 
            self.ind = (self.ind + 1) % self.slices
        elif event.key == "down": 
            self.ind = (self.ind - 1) % self.slices  
        elif event.key == "right":
            if self.index < 2 :
                self.index += 1
                if self.index == 1:
                    self.fig.suptitle("Action will be used on sagittal axe", fontsize=14)
                elif self.index == 2:
                    self.fig.suptitle("Action will be used on axial axe", fontsize=14)
            elif self.index == 2:
                self.index = 0
                self.fig.suptitle("Action will be used on coronal axe", fontsize=14)
        self.updateSlice()

    def onScroll(self, event):
        if event.button == 'up':
            if self.index == 0:
                self.x_ind = (self.x_ind + 1) % self.slices
            elif self.index == 1:
                self.y_ind = (self.y_ind + 1) % self.slices
            elif self.index == 2:
                self.z_ind = (self.z_ind + 1) % self.slices
        else:
            if self.index == 0:
                self.x_ind = (self.x_ind - 1) % self.slices
            elif self.index == 1:
                self.y_ind = (self.y_ind - 1) % self.slices
            elif self.index == 2:
                self.z_ind = (self.z_ind - 1) % self.slices        
        self.updateSlice()

    def onClick(self, event):
#        print('Button pressed: %s at X: %s, Y: %s'%(event.button, event.xdata, event.ydata))
        if event.button == 3:   #Reset image on right click
            if self.index == 0:
                self.y_ind = round(event.ydata)
                self.z_ind = round(event.xdata)
            elif self.index == 1:
                self.x_ind = round(event.ydata)
                self.z_ind = round(event.xdata)
            elif self.index == 2:
                self.x_ind = round(event.ydata)
                self.y_ind = round(event.xdata)
            self.updateAllSlice()
            
            
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
        if self.mouseClicked == None: #return        #if mouse isn't clicked ignore movement
            if self.plane == 'cross line':
                if event.xdata != None:
                    if event.ydata != None:

                        if self.index == 0:  
                            #IM space
                            imageDataX = self.X[self.x_ind,:,:]
                            s = self.X[self.x_ind,:,:].shape
                            arr1 = np.ones(s[0])
                            arr2 = np.ones(s[1])
                            arr3 = np.zeros((s)) # This array is only zeros with a vertical and horizontal white lines
                            arr3[int(event.ydata), :] = arr2
                            arr3[:, int(event.xdata)] = arr1

                            if self.x_ind >= np.size(self.X,0):
                                self.x_ind = np.size(self.X, 0)-1
                            self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
                            self.slices = np.size(self.X,0)
                            self.imx.set_data(imageDataX + 1000000*arr3)
                            self.ax[0,0].set_title('Slice cor %s/%s' % (self.x_ind,self.slices))
                            self.imx.axes.figure.canvas.draw()
                            
                            self.y_ind = int(event.ydata)
                            imageDataY = self.X[:,self.y_ind,:]
                            if self.y_ind >= np.size(self.X,0):
                                self.y_ind = np.size(self.X, 0)-1
                            self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
                            self.slices = np.size(self.X,0)
                            self.imy.set_data(imageDataY)
                            self.ax[0,1].set_title('Slice sag %s/%s' % (self.y_ind,self.slices))
                            self.imy.axes.figure.canvas.draw()
                                                        
                            #IM space
                            self.z_ind = int(event.xdata)
                            imageDataZ = self.X[:,:,self.z_ind]
                            if self.z_ind >= np.size(self.X,0):
                                self.z_ind = np.size(self.X, 0)-1
                            self.imz.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
                            self.slices = np.size(self.X,0)
                            self.imz.set_data(imageDataZ)
                            self.ax[0,2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
                            self.imz.axes.figure.canvas.draw()
                            
                        elif self.index == 1:
                            #IM space
                            imageDataY = self.X[:,self.y_ind,:]
                            s = self.X[:,self.y_ind,:].shape
                            arr1 = np.ones(s[0])
                            arr2 = np.ones(s[1])
                            arr3 = np.zeros((s)) # This array is only zeros with a vertical and horizontal white lines
                            arr3[int(event.ydata), :] = arr2
                            arr3[:, int(event.xdata)] = arr1

                            if self.y_ind >= np.size(self.X,0):
                                self.y_ind = np.size(self.X, 0)-1
                            self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
                            self.slices = np.size(self.X,0)
                            self.imy.set_data(imageDataY + 1000000*arr3)
                            self.ax[0,1].set_title('Slice sag %s/%s' % (self.y_ind,self.slices))
                            self.imy.axes.figure.canvas.draw()
                            
                            #IM space
                            self.x_ind = int(event.ydata)
                            imageDataX = self.X[self.x_ind,:,:]
                            if self.x_ind >= np.size(self.X,0):
                                self.x_ind = np.size(self.X, 0)-1
                            self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
                            self.slices = np.size(self.X,0)
                            self.imx.set_data(imageDataX)
                            self.ax[0,0].set_title('Slice cor %s/%s' % (self.x_ind,self.slices))
                            self.imx.axes.figure.canvas.draw()
                            
                            #IM space
                            self.z_ind = int(event.xdata)
                            imageDataZ = self.X[:,:,self.z_ind]
                            if self.z_ind >= np.size(self.X,0):
                                self.z_ind = np.size(self.X, 0)-1
                            self.imz.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
                            self.slices = np.size(self.X,0)
                            self.imz.set_data(imageDataZ)
                            self.ax[0,2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
                            self.imz.axes.figure.canvas.draw()
                            
                        elif self.index == 2:
                            #IM space
                            imageDataZ = self.X[:,:,self.z_ind]
                            s = self.X[:,:,self.z_ind].shape
                            arr1 = np.ones(s[0])
                            arr2 = np.ones(s[1])
                            arr3 = np.zeros((s)) # This array is only zeros with a vertical and horizontal white lines
                            arr3[int(event.ydata), :] = arr2
                            arr3[:, int(event.xdata)] = arr1

                            if self.z_ind >= np.size(self.X,0):
                                self.z_ind = np.size(self.X, 0)-1
                            self.imz.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
                            self.slices = np.size(self.X,0)
                            self.imz.set_data(imageDataZ + 1000000*arr3)
                            self.ax[0,2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
                            self.imz.axes.figure.canvas.draw()
                            
                            #IM space
                            self.x_ind = int(event.ydata)
                            imageDataX = self.X[self.x_ind,:,:]
                            if self.x_ind >= np.size(self.X,0):
                                self.x_ind = np.size(self.X, 0)-1
                            self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
                            self.slices = np.size(self.X,0)
                            self.imx.set_data(imageDataX)
                            self.ax[0,0].set_title('Slice cor %s/%s' % (self.x_ind,self.slices))
                            self.imx.axes.figure.canvas.draw()
                            
                            #IM space
                            self.y_ind = int(event.xdata)
                            imageDataY = self.X[:,self.y_ind,:]
                            if self.y_ind >= np.size(self.X,0):
                                self.y_ind = np.size(self.X, 0)-1
                            self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
                            self.slices = np.size(self.X,0)
                            self.imy.set_data(imageDataY)
                            self.ax[0,1].set_title('Slice sag %s/%s' % (self.y_ind,self.slices))
                            self.imy.axes.figure.canvas.draw()
                            
        elif self.mouseClicked != None:

            dx = event.xdata - self.mouseClicked[0]
            dy = event.ydata - self.mouseClicked[1]

            normDx = dx/self.mouseClicked[0]
            normDy = dy/self.mouseClicked[1]

            self.tempCmapRange = self.cmapRange*(1+normDy)
            self.tempCmapCenter = self.cmapCenter*(1+normDx)

            if self.index == 0:
                self.imx.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
                self.imx.axes.figure.canvas.draw()
                self.imxK.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
                self.imxK.axes.figure.canvas.draw()
            elif self.index == 1:
                self.imy.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
                self.imy.axes.figure.canvas.draw()
                self.imyK.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
                self.imyK.axes.figure.canvas.draw()
            elif self.index == 2:
                self.imz.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
                self.imz.axes.figure.canvas.draw()
                self.imzK.set_clim(self.tempCmapCenter-self.tempCmapRange/2, self.tempCmapCenter+self.tempCmapRange/2)
                self.imzK.axes.figure.canvas.draw()
        
    def updateSlice(self):
        imageDataX = self.X[self.x_ind,:,:]
        imageDataY = self.X[:,self.y_ind,:]
        imageDataZ = self.X[:,:,self.z_ind]
        
        imageK_DataX = self.K[self.x_ind,:,:]
        imageK_DataY = self.K[:,self.y_ind,:]
        imageK_DataZ = self.K[:,:,self.z_ind]
        
        if self.index == 0:
            if self.x_ind >= np.size(self.X,0):
                self.x_ind = np.size(self.X, 0)-1
            imageDataX = self.X[self.x_ind,:,:]
            self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
            self.slices = np.size(self.X,0)
            self.imx.set_data(imageDataX)
            self.ax[0,0].set_title('Slice cor %s/%s' % (self.x_ind,self.slices))
            self.imx.axes.figure.canvas.draw()
            
            if self.x_ind >= np.size(self.X,0):
                self.x_ind = np.size(self.X, 0)-1
            imageK_DataX = self.K[self.x_ind,:,:]
            self.imxK.set_extent((-0.5,np.size(imageK_DataX,1) - 0.5,np.size(imageK_DataX,0) - 0.5,-0.5))
            self.slices = np.size(self.X,0)
            self.imxK.set_data(imageK_DataX)
            self.ax[1,0].set_title('Slice cor %s/%s' % (self.x_ind,self.slices))
            self.imxK.axes.figure.canvas.draw()
            
        elif self.index == 1:
            if self.y_ind >= np.size(self.X,1):
                self.y_ind = np.size(self.X, 1)-1
            imageDataY = self.X[:,self.y_ind,:]
            self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
            self.slices = np.size(self.X,1)
            self.imy.set_data(imageDataY)
            self.ax[0,1].set_title('Slice sag %s/%s' % (self.y_ind,self.slices))
            self.imy.axes.figure.canvas.draw()
            
            if self.y_ind >= np.size(self.X,1):
                self.y_ind = np.size(self.X, 1)-1
            imageK_DataY = self.K[:,self.y_ind,:]
            self.imyK.set_extent((-0.5,np.size(imageK_DataY,1) - 0.5,np.size(imageK_DataY,0) - 0.5,-0.5))
            self.slices = np.size(self.X,1)
            self.imyK.set_data(imageK_DataY)
            self.ax[1,1].set_title('Slice sag %s/%s' % (self.y_ind,self.slices))
            self.imyK.axes.figure.canvas.draw()
            
        elif self.index == 2:
            if self.z_ind >= np.size(self.X,2):
                self.z_ind = np.size(self.X,2)-1
            imageDataZ = self.X[:,:,self.z_ind]
            self.imz.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
            self.slices = np.size(self.X,2)
            self.imz.set_data(imageDataZ)        
            self.ax[0,2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
            self.imz.axes.figure.canvas.draw()
            
            if self.z_ind >= np.size(self.X,2):
                self.z_ind = np.size(self.X,2)-1
            imageK_DataZ = self.K[:,:,self.z_ind]
            self.imzK.set_extent((-0.5,np.size(imageK_DataZ,1) - 0.5,np.size(imageK_DataZ,0) - 0.5,-0.5))
            self.slices = np.size(self.X,2)
            self.imzK.set_data(imageK_DataZ)        
            self.ax[1,2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
            self.imzK.axes.figure.canvas.draw()
            
    def updateAllSlice(self):
        imageDataX = self.X[self.x_ind,:,:]
        imageDataY = self.X[:,self.y_ind,:]
        imageDataZ = self.X[:,:,self.z_ind]
        
        if self.x_ind >= np.size(self.X,0):
            self.x_ind = np.size(self.X, 0)-1
        imageDataX = self.X[self.x_ind,:,:]
        self.imx.set_extent((-0.5,np.size(imageDataX,1) - 0.5,np.size(imageDataX,0) - 0.5,-0.5))
        self.slices = np.size(self.X,0)
        self.imx.set_data(imageDataX)
        self.ax[0,0].set_title('Slice cor %s/%s' % (self.x_ind,self.slices))
        self.imx.axes.figure.canvas.draw()

        if self.y_ind >= np.size(self.X,1):
            self.y_ind = np.size(self.X, 1)-1
        imageDataY = self.X[:,self.y_ind,:]
        self.imy.set_extent((-0.5,np.size(imageDataY,1) - 0.5,np.size(imageDataY,0) - 0.5,-0.5))
        self.slices = np.size(self.X,1)
        self.imy.set_data(imageDataY)
        self.ax[0,1].set_title('Slice sag %s/%s' % (self.y_ind,self.slices))
        self.imy.axes.figure.canvas.draw()

        if self.z_ind >= np.size(self.X,2):
            self.z_ind = np.size(self.X,2)-1
        imageDataZ = self.X[:,:,self.z_ind]
        self.imz.set_extent((-0.5,np.size(imageDataZ,1) - 0.5,np.size(imageDataZ,0) - 0.5,-0.5))
        self.slices = np.size(self.X,2)
        self.imz.set_data(imageDataZ)        
        self.ax[0,2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
        self.imz.axes.figure.canvas.draw()
        
        ## Kspace
        
        imageK_DataX = self.K[self.x_ind,:,:]
        imageK_DataY = self.K[:,self.y_ind,:]
        imageK_DataZ = self.K[:,:,self.z_ind]
        
        if self.x_ind >= np.size(self.X,0):
            self.x_ind = np.size(self.X, 0)-1
        imageK_DataX = self.K[self.x_ind,:,:]
        self.imxK.set_extent((-0.5,np.size(imageK_DataX,1) - 0.5,np.size(imageK_DataX,0) - 0.5,-0.5))
        self.slices = np.size(self.X,0)
        self.cmapRange = np.nanmax(self.K) - np.nanmin(self.K)
        self.cmapCenter = np.nanmin(self.K) + self.cmapRange/2
        self.tempCmapRange = self.cmapRange
        self.tempCmapCenter = self.cmapCenter
        self.imxK = self.ax[1,0].imshow(imageK_DataX, cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)        
        self.ax[1,0].set_title('Slice cor %s/%s' % (self.x_ind,self.slices))
        self.imxK.axes.figure.canvas.draw()

        if self.y_ind >= np.size(self.X,1):
            self.y_ind = np.size(self.X, 1)-1
        imageK_DataY = self.K[:,self.y_ind,:]
        self.imyK.set_extent((-0.5,np.size(imageK_DataY,1) - 0.5,np.size(imageK_DataY,0) - 0.5,-0.5))
        self.slices = np.size(self.X,1)
        self.cmapRange = np.nanmax(self.K) - np.nanmin(self.K)
        self.cmapCenter = np.nanmin(self.K) + self.cmapRange/2
        self.tempCmapRange = self.cmapRange
        self.tempCmapCenter = self.cmapCenter
        self.imyK = self.ax[1,1].imshow(imageK_DataY, cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)  
        self.ax[1,1].set_title('Slice sag %s/%s' % (self.y_ind,self.slices))
        self.imyK.axes.figure.canvas.draw()

        if self.z_ind >= np.size(self.X,2):
            self.z_ind = np.size(self.X,2)-1
        imageK_DataZ = self.K[:,:,self.z_ind]
        self.imzK.set_extent((-0.5,np.size(imageK_DataZ,1) - 0.5,np.size(imageK_DataZ,0) - 0.5,-0.5))
        self.slices = np.size(self.X,2)
        self.cmapRange = np.nanmax(self.K) - np.nanmin(self.K)
        self.cmapCenter = np.nanmin(self.K) + self.cmapRange/2
        self.tempCmapRange = self.cmapRange
        self.tempCmapCenter = self.cmapCenter
        self.imzK = self.ax[1,2].imshow(imageK_DataZ, cmap = self.colorMap,vmin = self.cmapCenter-self.cmapRange/2, vmax= self.cmapCenter+self.cmapRange/2)  
        self.ax[1,2].set_title('Slice ax %s/%s' % (self.z_ind,self.slices))
        self.imzK.axes.figure.canvas.draw()


