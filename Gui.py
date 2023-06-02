import numpy as np 
import os
import matplotlib.pyplot as plt
import interactive3Dplot as plt3D
import interactiveKspace3Dplot as pltK3D
import interactiveROIbox as pltroibox
import interactiveROIdraw as pltroidraw
import interactiveAxesPlot as pltaxes
import interactiveKspaceAxesPlot as pltKaxes
import interactiveSixAxesPlot as pltsixaxes
import interactivePlot_SNR as pltsnr
import imageProcessing as imProc
import EvoDataProcessing as evoProc
#import distortionCorrection as distCorr
from tkinter import *
from tkinter import filedialog
from PIL import ImageTk, Image
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from mpl_interactions import ipyplot as iplt

root = Tk()
root.title("Low field reconstruction GUI")

#/////////////////////// FUNCTIONS ///////////////////////#

##### Functions for the tutorial #####    
# Forward function for the next tutorial image
def forward(image_num):
    """
    Will change the slide in the tutorial to the next one
    
    Input:  image_num --> the number of the slide
    Output: none
    """
    # To be able to update the buttons, need global buttons
    global Image_label
    global forward_button
    global back_button

    Image_label.grid_forget()                                    # To remove the original image in my_label when we click forward
    Image_label = Label(Tutorial, image=image_list[image_num-1]) # Going from 1st to 2nd image means going from 0 to 1, we passed 2 so need a -1

    # Updates, so that forward doesn't have the '2' as input
    forward_button = Button(Tutorial, text=">>", command=lambda: forward(image_num+1)) # Want the next image so +1
    back_button = Button(Tutorial, text="<<", command=lambda: back(image_num-1))       # Want the previous image so -1

    if image_num == len(image_list):
        forward_button = Button(Tutorial, text=">>", state=DISABLED)

    Image_label.grid(row=0,column=0,columnspan=3)
    back_button.grid(row=1,column=0)
    forward_button.grid(row=1,column=2)
    
    # We want the status to be updated is we press forward button
    status = Label(Tutorial, text="Image " + str(image_num)  + " of " + str(len(image_list)), bd=1, relief=SUNKEN, anchor=E)
    status.grid(row=2,column=0, columnspan=3, sticky=W+E)

# Backward function for the next tutorial image
def back(image_num):
    """
    Will change the slide in the tutorial to the previous one
    
    Input:  image_num --> the number of the slide
    Output: none
    """
    # To be able to update the buttons, need global buttons
    global Image_label
    global forward_button
    global back_button

    Image_label.grid_forget()
    Image_label = Label(Tutorial, image=image_list[image_num-1])    # Going from 1st to 2nd image means going from 0 to 1, we passed 2 so need a -1

    # Updates
    forward_button = Button(Tutorial, text=">>", command=lambda: forward(image_num+1)) # Want the next image so +1
    back_button = Button(Tutorial, text="<<", command=lambda: back(image_num-1))       # Want the previous image so -1

    if image_num == 1:
        back_button = Button(Tutorial, text="<<", state=DISABLED)

    Image_label.grid(row=0,column=0,columnspan=3)
    back_button.grid(row=1,column=0)
    forward_button.grid(row=1,column=2)
    
    # We want the status to be updated is we press back button
    status = Label(Tutorial, text="Image " + str(image_num)  + " of " + str(len(image_list)), bd=1, relief=SUNKEN, anchor=E)
    status.grid(row=2,column=0, columnspan=3, sticky=W+E)

# Function starting the tutorial
def tutorial():
    """
    Will open a new window with the first slide of the tutorial and enable the user to change slides or exit the tutorial
    
    Input: none 
    Output: none 
    """
    global Tutorial
    global Image_label
    global image_list
    Tutorial = Toplevel() # New window
    Tutorial.title("Tutorial")
    
    # The images have to come from the correct folder
    image1 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide1.png"))
    image2 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide2.png"))
    image3 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide3.png"))
    image4 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide4.png"))
    image5 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide5.png"))
    image6 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide6.png"))
    image7 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide7.png"))
    image8 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide8.png"))
    image9 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide9.png"))
    image10 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide10.png"))
    image11 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide11.png"))
    image12 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide12.png"))
    image13 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide13.png"))
    image14 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide14.png"))
    image15 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide15.png"))
    image16 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide16.png"))
    image17 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide17.png"))
    image18 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide18.png"))
    image19 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide19.png"))
    image20 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide20.png"))
    image21 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide21.png"))
    image22 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide22.png"))
    image23 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide23.png"))
    image24 = ImageTk.PhotoImage(Image.open(r"Tutorial images\slides\Slide24.png"))
    
    # Creating a list of all images
    image_list = [image1, image2, image3, image4, image5, image6, image7, image8, image9, image10, image11, image12, image13, image14, image15, image16, image17, image18, image19, image20, image21, image22, image23, image24] 
    
    # The status, bd for border, relief for sunnken style, anchor E -> so it is fixed at buttom right
    status = Label(Tutorial, text="Image 1 of " + str(len(image_list)), bd=1, relief=SUNKEN, anchor=E)
    # sticky for stretching, here from west to east
    status.grid(row=2,column=0, columnspan=3, sticky=W+E)
    
    Image_label = Label(Tutorial, image=image1) # Show the first image
    Image_label.grid(row=0,column=0,columnspan=3)
    
    # Buttons
    back_button = Button(Tutorial, text="<<", command=back, state = DISABLED) # We don't want to be able to go back when the first image arrives
    exit_button = Button(Tutorial, text="Exit tutorial", command=Tutorial.destroy)
    forward_button = Button(Tutorial, text=">>", command=lambda: forward(2)) # We pass '2' to go to the next image
    back_button.grid(row=1,column=0)
    exit_button.grid(row=1,column=1)
    forward_button.grid(row=1,column=2)
    
# Function computing SNR based on 2 boxes defined by the parameters (box "s" is in a noisy region of the image (to compute the std of noise) and box "m" is in a image region where their is signal (to compute the mean of the signal))
def snr_homemade(im,s1,s2,s3,s4,m1,m2,m3,m4):
    """
    Computes the SNR of an image
    
    Input: im --> image
           s1,s2 --> top and bottum of "noise box" (height)
           s3,s4 --> left and rigth of "noise box" (width)
           m1,m2 --> top and bottum of "signal box" (height)
           m3,m4 --> left and rigth of "signal box" (width)
           
    Output: SNR 
    """
    m = np.mean(im[m1:m2,m3:m4])
    std = np.std(im[s1:s2,s3:s4])
    
    if std == 0:
        snr = 1000
    else:
        snr = m/std    
    return snr, m, std
    

#### Function showing the user's selection of parameters
def show():
    global var_scan_num
    global locations
    global file_path_label
    
    print('////////////////////SELECTION/////////////////////////')
    print('Filter selected -- > ' + str(fil.get()))
    print('Filter strength -- > ' + str(fil_strength.get()))
    print('B0 correction --> ' + str(B0.get()))
    print('Gradient correction --> ' + str(G.get()))
    print('Noise correction --> ' + str(N.get()))
    print('Zero filling (yes or no) -- > ' + str(Z.get()))
    
    if Z.get() == "1":
        print('zero-filling dim 1 --> ' + str(z1_entry.get()))
        print('zero-filling dim 2 --> ' + str(z2_entry.get()))
        print('zero-filling dim 3 --> ' + str(z3_entry.get()))
    
    print('Directory to file --> ' + str(file_path_label))
        
    if var_scan_num.get() != "None":
        x = var_scan_num.get()
        x_elements = var_scan_num.get().split(" ")
        scan_num = x_elements[len(x_elements)-1]   # Scan number
        datafol = x_elements[len(x_elements)-6]   
        path = file_path_label + '/' + datafol     # Path to folder containing that scan number
        print('Scan number -- > ' + str(scan_num))
        print('Path to datafolder -- > ' + str(path))
    else:
        print('Scan number -- > not selected yet')
        print('Path to datafolder -- > not selected yet')
    
    print('//////////////////////////////////////////////////////')
     
##### Scan selection 
def open():
    global file_path_label
    global var_scan_num
    global parFile
    global scan_num
    global path
    global drop
    
    # Will return the directory of the file
    root.filename = filedialog.askdirectory()
    file_path_label = str(root.filename)
               
    parFile = []
    numbers = []
    locations = []
    options = []
    index = 0
    for (Root,dirs,files) in os.walk(file_path_label):
        for file in files:
            if file.endswith('.mrd'):
                parFile.append(os.path.join(Root, file))
                separate = parFile[index].split("\\")
                numbers.append(separate[len(separate)-2])   # Stores the number of the scans that have a "data.3d" file
                locations.append(separate[len(separate)-3]) # Stores the file of the scans that have a "data.3d" file
                options.append(' Folder --> ' + locations[index])
                index += 1
    
    if index > 0: # A dropdown menu will appear only if their is at least one 'data.3d' file
        var_scan_num.set("None") # Default value
        drop = OptionMenu(frame1, var_scan_num, *options) # Using a list
        drop.grid(row = 8, column = 2, columnspan = 4)    
    elif index == 0:
        var_scan_num.set("None") # Default value
        options = ['no ".mrd" files in this folder']
        drop = OptionMenu(frame1, var_scan_num, *options) 
        drop.grid(row = 8, column = 2, columnspan = 4)
        
##### Applying the reconstruction code
def recon():
    global file_path_label
    global scanParams
    global reconImage
    global var_scan_num
    global B0
    global G
    global N
    global reconImage
    global kSpace
    
    # To acces the selection
    x_elements = var_scan_num.get().split(" ")  
    datafol = x_elements[len(x_elements)-1]   
    path = file_path_label + '/' + datafol     # Path to folder containing that scan number

    dataFolder = path
    
    scanParams = evoProc.readPar(dataFolder)
    
    '''Read k-space data'''
    kSpace = evoProc.readKSpace(dataFolder, scanParams = scanParams, correctOversampling = True)
    
    '''Noise correction'''
    if N.get() == "yes": 
        kSpace              = imProc.noiseCorrection(kSpace, scanParams, dataFolder)
    
    '''apply filter to k-space data'''
    if str(fil.get()) == "Sine bell": 
        kSpace          = imProc.sineBellSquaredFilter(kSpace, filterStrength = fil_strength.get()/100) # strength = 0 no filter, 1 = max
    elif str(fil.get()) == "Gaussian": 
        kSpace          = imProc.gaussianFilter(kSpace, filterStrength = fil_strength.get()/100)
    
    '''Distortion correction'''
    reconImage          = np.fft.fftshift(np.fft.fftn((np.fft.fftshift(kSpace))))
    
    ''' B0 correction '''
    if B0.get() == "yes":
        reconImage          = distCorr.b0Correction(kSpace, scanParams, dataFolder, shOrder = 15)
    
    ''' Gradient correction '''
    if G.get() == "yes":
        reconImage          = distCorr.gradUnwarp(reconImage, scanParams, dataFolder)

    if Z.get() == "1":
        '''post processing''' #data should revert back to k-space for this
        kSpace = np.fft.fftshift(np.fft.ifftn((np.fft.fftshift(reconImage))))

        '''Zero fill data'''
        a = int(z1_entry.get()); b = int(z2_entry.get()); c = int(z3_entry.get())
        kSpace          = imProc.zeroFill(kSpace, (a,b,c))

        ''' Shift image'''
        # shiftDistance = 0.025 #shift distance in meters
        # shiftAxis = 2 # FE = 0, phase 1 = 1, phase 2 = 2
        # kSpace = imProc.shiftImage(kSpace, scanParams, shiftDistance,shiftAxis)

        reconImage       = np.fft.fftshift(np.fft.fftn((np.fft.fftshift(kSpace))))

    kSpace = np.fft.fftshift(np.fft.ifftn((np.fft.fftshift(reconImage))))
    
    recon_label = Label(root, text = "Reconstruction done").grid(row = 3, column = 1)

    
def display():
    global reconImage
    global kSpace 
    
    method = im.get()
    
    if method == "One axe (IM space)":
        fig, ax = plt.subplots()
        fig3D = plt3D.interactivePlot(fig, ax, np.abs(reconImage), str(fil.get()), plotAxis = 2, fov = scanParams["FOV"])
        plt.show()
    elif method == "Three axes (IM space)":
        fig, ax = plt.subplots(1,3)
        fig3D = pltaxes.interactivePlot(fig, ax, np.abs(reconImage), str(fil.get()), plotAxis = 2, fov = scanParams["FOV"])
        plt.show()
    elif method == "One axe (K space)":
        fig, ax = plt.subplots()
        fig3D = pltK3D.interactivePlot(fig, ax, np.abs(kSpace), str(fil.get()), plotAxis = 2, fov = scanParams["FOV"])
        plt.show()
    elif method == "Three axes (K space)":
        fig, ax = plt.subplots(1,3)
        fig3D = pltKaxes.interactivePlot(fig, ax, np.abs(kSpace), str(fil.get()), plotAxis = 2, fov = scanParams["FOV"])
        plt.show()
    elif method == "Six axes (IM + K space)":
        fig, ax = plt.subplots(2,3)
        fig3D = pltsixaxes.interactivePlot(fig, ax, np.abs(reconImage), np.abs(kSpace), str(fil.get()), plotAxis = 2, fov = scanParams["FOV"])
        plt.show()
    elif method == "ROI box":
        roi = 'off'
        fig, ax = plt.subplots()
        fig3D = pltroibox.interactivePlot(fig, ax, np.abs(reconImage), str(fil.get()), roi, plotAxis = 2, fov = scanParams["FOV"])
        plt.show()
    elif method == "ROI draw":
        roi = 'off'
        fig, ax = plt.subplots()
        fig3D = pltroidraw.interactivePlot(fig, ax, np.abs(reconImage), str(fil.get()), roi, plotAxis = 2, fov = scanParams["FOV"])
        plt.show()        
        
def display_all():
    global dis
    global reconImage
    
    shape = reconImage.shape
    roots = np.floor(np.sqrt(shape))
    
    if dis.get() == 'axial':
        root = int(roots[2])
        if root**2 == shape[2]: # perfect square
            fig, ax = plt.subplots(root,root)
            index = 0
            for i in range(root):
                for j in range(root):
                    ax[i, j].imshow(np.abs(reconImage[:,:,index]),cmap='gray')
                    ax[i, j].axis('off'); #ax[i, j].set_title(str(index))
                    index += 1      
            
        elif root**2 != shape[2]: # not perfect square
            fig, ax = plt.subplots(root+1,root+1)
            index = 0
            for i in range(root+1):
                for j in range(root+1):
                    if index < shape[2]:
                        ax[i, j].imshow(np.abs(reconImage[:,:,index]),cmap='gray')
                        ax[i, j].axis('off'); #ax[i, j].set_title(str(index))
                        index += 1   
                    else:
                        ax[i, j].axis('off')
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0.015)
        plt.show()   
                
    elif dis.get() == 'sagittal':      
        root = int(roots[1])      
        if root**2 == shape[1]: # perfect square
            fig, ax = plt.subplots(root,root)
            index = 0
            for i in range(root):
                for j in range(root):
                    ax[i, j].imshow(np.abs(reconImage[:,index,:]),cmap='gray')
                    ax[i, j].axis('off'); #ax[i, j].set_title(str(index))
                    index += 1
            
        elif root**2 != shape[1]:
            fig, ax = plt.subplots(root+1,root+1)
            index = 0
            for i in range(root+1):
                for j in range(root+1):
                    if index < shape[1]:
                        ax[i, j].imshow(np.abs(reconImage[:,index,:]),cmap='gray')
                        ax[i, j].axis('off'); #ax[i, j].set_title(str(index))
                        index += 1
                    else:
                        ax[i, j].axis('off')
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0.015)
        plt.show()
        
    elif dis.get() == 'coronal':
        root = int(roots[0])
        if root**2 == shape[0]: # perfect square
            fig, ax = plt.subplots(root,root)
            index = 0
            for i in range(root):
                for j in range(root):
                    ax[i, j].imshow(np.abs(reconImage[index,:,:]),cmap='gray')
                    ax[i, j].axis('off'); #ax[i, j].set_title(str(index))
                    index += 1
            
        elif root**2 != shape[0]:
            fig, ax = plt.subplots(root+1,root+1)
            index = 0
            for i in range(root+1):
                for j in range(root+1):
                    if index < shape[0]:
                        ax[i, j].imshow(np.abs(reconImage[index,:,:]),cmap='gray')
                        ax[i, j].axis('off'); #ax[i, j].set_title(str(index))
                        index += 1
                    else:
                        ax[i, j].axis('off')
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0.015)
        plt.show()
          
def display_selected():
    global axe
    global Start_entry
    global End_entry
    global reconImage
    global kSpace 
    
    shape = reconImage.shape
    
    S = int(Start_entry.get())
    E = int(End_entry.get())
    
    if S<E:
        
        num = E-S    
        root = int(np.floor(np.sqrt(num)))

        if axe.get() == 'axial':
            if root**2 == num: # perfect square
                fig, ax = plt.subplots(root,root)
                index = S-1
                for i in range(root):
                    for j in range(root):
                        ax[i, j].imshow(np.abs(reconImage[:,:,index]),cmap='gray')
                        ax[i, j].axis('off'); #ax[i, j].set_title(str(index))
                        index += 1      
            elif root**2 != num: # not perfect square
                fig, ax = plt.subplots(root+1,root+1)
                index = S-1
                for i in range(root+1):
                    for j in range(root+1):
                        if index < E-1:
                            ax[i, j].imshow(np.abs(reconImage[:,:,index]),cmap='gray')
                            ax[i, j].axis('off'); #ax[i, j].set_title(str(index))
                            index += 1   
                        else:
                            ax[i, j].axis('off')
            plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0.015)
            plt.show()  

        elif axe.get() == 'sagittal':
            if root**2 == num: # perfect square
                fig, ax = plt.subplots(root,root)
                index = S-1
                for i in range(root):
                    for j in range(root):
                        ax[i, j].imshow(np.abs(reconImage[:,index,:]),cmap='gray')
                        ax[i, j].axis('off'); #ax[i, j].set_title(str(index))
                        index += 1      
            elif root**2 != num: # not perfect square
                fig, ax = plt.subplots(root+1,root+1)
                index = S-1
                for i in range(root+1):
                    for j in range(root+1):
                        if index < E-1:
                            ax[i, j].imshow(np.abs(reconImage[:,index,:]),cmap='gray')
                            ax[i, j].axis('off'); #ax[i, j].set_title(str(index))
                            index += 1   
                        else:
                            ax[i, j].axis('off')
            plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0.015)
            plt.show()

        elif axe.get() == 'coronal':
            if root**2 == num: # perfect square
                fig, ax = plt.subplots(root,root)
                index = S-1
                for i in range(root):
                    for j in range(root):
                        ax[i, j].imshow(np.abs(reconImage[index,:,:]),cmap='gray')
                        ax[i, j].axis('off'); #ax[i, j].set_title(str(index))
                        index += 1      
            elif root**2 != num: # not perfect square
                fig, ax = plt.subplots(root+1,root+1)
                index = S-1
                for i in range(root+1):
                    for j in range(root+1):
                        if index < E-1:
                            ax[i, j].imshow(np.abs(reconImage[index,:,:]),cmap='gray')
                            ax[i, j].axis('off'); #ax[i, j].set_title(str(index))
                            index += 1   
                        else:
                            ax[i, j].axis('off')
            plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0.015)
            plt.show()

#/////////////////////// FRAMES ///////////////////////#
global fil
global fil_strength
global file_path_label
global scanParams
global reconImage
global scan_num
global var_scan_num
global path
global B0
global G
global N

file_path_label = 'None'

############### Frame regarding the post-processing ###############
frame1 = LabelFrame(root, text = "Post-processing")
frame1.grid(row = 0, column = 0, rowspan = 4)

# Labels of the parameters
fliter_selection_label = Label(frame1, text = "Filter selection ").grid(row = 1, column = 0)
fliter_strength_label = Label(frame1, text = "Filter strength ").grid(row = 1, column = 2, columnspan = 3)
zero_filling_label = Label(frame1, text = "Zero filling").grid(row = 5, column = 0)
z1_label = Label(frame1, text = "Zero filling dimension ").grid(row = 7, column = 0)
z2_label = Label(frame1, text = "x").grid(row = 7, column = 2)
z3_label = Label(frame1, text = "x").grid(row = 7, column = 4)
noise_correction_label = Label(frame1, text = "Noise correction").grid(row = 5, column = 2, columnspan = 3)
file_selection_label = Label(frame1, text = "File selection").grid(row = 8, column = 0)
B0_selection_label = Label(frame1, text = "B0 correction").grid(row = 3, column = 0)
Gradient_selection_label = Label(frame1, text = "Gradient correction").grid(row = 3, column = 2, columnspan = 3)

#Filter selection 
fil = StringVar()
fil.set("Sine bell") # To set what is inside before somebody click on it
Radiobutton(frame1, text="Sine bell", variable=fil, value="Sine bell").grid(row = 1, column = 1)
Radiobutton(frame1, text="Gaussian", variable=fil, value="Gaussian").grid(row = 2, column = 1)

#Filter strength
fil_strength = Scale(frame1, from_=0, to=100, orient=HORIZONTAL)
fil_strength.grid(row = 1, column = 5)

#Zero filling
Z = StringVar()
Z.set("0") # To set what is inside before somebody click on it
Radiobutton(frame1, text="Yes", variable=Z, value="1").grid(row = 5, column = 1)
Radiobutton(frame1, text="No", variable=Z, value="0").grid(row = 6, column = 1)

z1_entry = Entry(frame1); z1_entry.grid(row = 7, column = 1); z1_entry.insert(0, '0')
z2_entry = Entry(frame1); z2_entry.grid(row = 7, column = 3); z2_entry.insert(0, '0')
z3_entry = Entry(frame1); z3_entry.grid(row = 7, column = 5); z3_entry.insert(0, '0')

#File selection
button = Button(frame1, text="Open file", command=open).grid(row = 8, column = 1)

var_scan_num = StringVar()
var_scan_num.set("None") # Default value
drop = OptionMenu(frame1, var_scan_num, *['None']) # Using a list
drop.grid(row = 8, column = 2, columnspan = 4)

#B0 correction selection
B0 = StringVar()
B0.set("no") # To set what is inside before somebody click on it
Radiobutton(frame1, text="yes", variable=B0, value="yes").grid(row = 3, column = 1)
Radiobutton(frame1, text="no", variable=B0, value="no").grid(row = 4, column = 1)

#Gradient correction selection
G = StringVar()
G.set("no") # To set what is inside before somebody click on it
Radiobutton(frame1, text="yes", variable=G, value="yes").grid(row = 3, column = 5)
Radiobutton(frame1, text="no", variable=G, value="no").grid(row = 4, column = 5)

#Noise correction
N = StringVar()
N.set("no") # To set what is inside before somebody click on it
Radiobutton(frame1, text="Yes", variable=N, value="yes").grid(row = 5, column = 5)
Radiobutton(frame1, text="No", variable=N, value="no").grid(row = 6, column = 5)


############### Frame regarding the display options ###############
frame3 = LabelFrame(root, text = "Display options")
frame3.grid(row = 4, column = 0, columnspan = 2)

#Plot selection
frame4 = LabelFrame(frame3, text = "3D data visualization")
frame4.grid(row = 0, column = 0)

plot_selection_label = Label(frame4, text = "Plot selection").grid(row = 0, column = 0)

# Dropdown menu to select visualization method   
options = [
    "One axe (IM space)",
    "Three axes (IM space)",
    "One axe (K space)",
    "Three axes (K space)",
    "Six axes (IM + K space)",
    "ROI box",
    "ROI draw"
]
im = StringVar()
im.set("No selection")                             # Default value, or could use; options[0]
image_visu_drop = OptionMenu(frame4, im, *options) # Using a list, NEEDS a star in front
image_visu_drop.grid(row = 0, column = 1)

plot_button = Button(frame4, text="Display plot", command=display).grid(row = 0, column = 2)

#All images display
frame5 = LabelFrame(frame3, text = "Display all data")
frame5.grid(row = 1, column = 0)

Display_all_label = Label(frame5, text = "Axe of images").grid(row = 0, column = 0)
dis = StringVar()
dis.set("axial") # To set what is inside before somebody click on it
Radiobutton(frame5, text="Axial", variable=dis, value="axial").grid(row = 0, column = 1)
Radiobutton(frame5, text="Coronal", variable=dis, value="coronal").grid(row = 1, column = 1)
Radiobutton(frame5, text="Sagittal", variable=dis, value="sagittal").grid(row = 2, column = 1)

display_button = Button(frame5, text="Display all images", command=display_all).grid(row = 0, column = 2)

#Selected images display
frame6 = LabelFrame(frame3, text = "Display selected data")
frame6.grid(row = 2, column = 0)

Display_selected_label = Label(frame6, text = "Select images").grid(row = 0, column = 0)
Display_selected_label = Label(frame6, text = "Start").grid(row = 0, column = 1)
Display_selected_label = Label(frame6, text = "End").grid(row = 0, column = 2)

Start_entry = Entry(frame6); Start_entry.grid(row = 1, column = 1); Start_entry.insert(0, '0')
End_entry = Entry(frame6); End_entry.grid(row = 1, column = 2); End_entry.insert(0, '0')

Display_all_label = Label(frame6, text = "Axe of images").grid(row = 2, column = 0)
axe = StringVar()
axe.set("axial") # To set what is inside before somebody click on it
Radiobutton(frame6, text="Axial", variable=axe, value="axial").grid(row = 2, column = 1)
Radiobutton(frame6, text="Coronal", variable=axe, value="coronal").grid(row = 2, column = 2)
Radiobutton(frame6, text="Sagittal", variable=axe, value="sagittal").grid(row = 2, column = 3)

display_button = Button(frame6, text="Display selected images", command=display_selected).grid(row = 3, column = 1)

#Frame regarding the SNR
frame7 = LabelFrame(frame3, text = "SNR")
frame7.grid(row = 0, column = 1)

def snr():
    fig, ax = plt.subplots() # figsize=(4, 4)
    figSNR = pltsnr.interactivePlot(fig, ax, np.abs(reconImage), str(fil.get()), np.abs(reconImage).shape, snr_homemade)        
    plt.show() 

SNR_button = Button(frame7, text="Visualize SNR", command=snr).grid(row = 0, column = 0)

############### Root ###############
# Tutorial button
Tuto_button = Button(root, text = "Show tutorial", command = tutorial).grid(row = 0, column = 1)
# Showing the parameters selected button
mybutton = Button(root, text="Show selection", command=show).grid(row = 1, column = 1)
# Apply reconstruction of the data
recon_button = Button(root, text="Reconstruct", command=recon).grid(row = 2, column = 1)
# Label telling the user that the reconstruction is finished
recon_label = Label(root, text = "Not yet reconstructed").grid(row = 3, column = 1)

root.mainloop()

