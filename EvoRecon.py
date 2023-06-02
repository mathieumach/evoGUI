import numpy as np 
import matplotlib.pyplot as plt
import interactive3Dplot as plt3D
import EvoDataProcessing as evoProc
import imageProcessing as imProc
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()
dataFolder = filedialog.askopenfilename(filetypes=[("MRI Files", ".MRD"),
                                    ("All Files", ".*")],
                                    title="Select MR Data File",
                                    multiple=False)

scanParams = evoProc.readPar(dataFolder)

print(dataFolder)

'''Read k-space data'''
kSpace = evoProc.readKSpace(dataFolder, scanParams = scanParams, correctOversampling = True)

'''Noise correction'''
# kSpace              = imProc.noiseCorrection(kSpace, scanParams, dataFolder)

'''apply filter to k-space data'''
kSpace = imProc.sineBellSquaredFilter(kSpace, filterStrength = 0.5) # strength = 0 no filter, 1 = max

'''Distortion correction'''
reconImage = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(kSpace)))
if int(scanParams['NO_VIEWS_2'][1]) >1:
    reconImage = np.fft.fftshift(reconImage, axes=2)

###################################        Plots           #######################################################
if 'NO_VIEWS_2' in scanParams:
    if int(scanParams['NO_VIEWS_2'][1]) >1:
        fig, ax = plt.subplots()
        fig3D = plt3D.interactivePlot(fig, ax, np.log1p(np.abs(kSpace)), plotAxis = 2, fov = scanParams["FOV"], axisLabels = None)       # add your FOV
        
        fig, ax = plt.subplots()
        fig3D = plt3D.interactivePlot(fig, ax, np.abs(reconImage), plotAxis = 2, fov = None, axisLabels = None)
        fig, ax = plt.subplots(2,10)
        for i in range(10):
            ax[0][i].imshow(np.abs(reconImage[:,:,i]), cmap = 'gray')
            ax[0][i].set_title(i)
            ax[0][i].axis('off')
            ax[1][i].imshow(np.log1p(np.abs(kSpace[:,:,i])), cmap = 'gray')
            ax[1][i].axis('off')
    else:
        fig, ax = plt.subplots()
        plt.imshow(np.abs(reconImage), cmap = 'gray')
        plt.title("magimage")
        
        fig, ax = plt.subplots()
        plt.imshow(np.log1p(np.abs(kSpace)), cmap = 'gray')
        plt.title("K space") 

else:
    fig, ax = plt.subplots(1,3)
    ax[0].imshow(np.angle(kSpace))
    ax[0].set_title("Phase")
    ax[1].imshow(np.abs(reconImage), cmap = 'gray')
    ax[1].set_title("magimage")
    ax[2].imshow(np.log1p(np.abs(kSpace), cmap = 'gray'))
    ax[2].set_title("K space")

plt.show()