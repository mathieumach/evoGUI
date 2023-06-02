# -*- coding: utf-8 -*-
"""
Created on Tue Jan  1 10:49:50 2023

@author: Javad Parsa

Version: 2.0
"""

import numpy as np
import re
import os
import pathlib 
import interactive3Dplot as plt3D

def readPar(file_path):
    with open(file_path, 'rb') as parFile:
        parameters = {}
        Parlist=[]
        lines = parFile.read().split(b':')
        lines.reverse()
        for line in lines:
            line = line.rstrip().decode("utf-8", "ignore")
            line = line.replace(',', '').replace('\n', '').replace('\r', '').replace('\0', '').replace('VAR ', '').replace('VAR_ARRAY ', '')
            parData = line.split(' ')
            if len(parData) >= 2:
                parameters[parData[0]] = parData[1:]
            Parlist.append(line)
            if re.search('PPL', line) :
                break
        Parlist.reverse()
        
        if 'NO_VIEWS_2' in parameters:
            try:
                parameters["FOV"] = [float(parameters['FOV'][0]), float(parameters['fov_sl'][0]), float(parameters['SLICE_THICKNESS'][2])]
            except KeyError:
                print("Error reading field of view, could not find all/some of the parameters: FOV, fov_sl, SLICE_THICKNESS")
                print("Defaulting to 256, 256, 256 mm field of view")
                parameters["FOV"] = [256, 256, 256]
            except ValueError:
                print("Error converting FOV, fov_sl, or SLICE_THICKNESS to float")
                print("Defaulting to 256, 256, 256 mm field of view")
                parameters["FOV"] = [256, 256, 256]

        filename = os.path.splitext(pathlib.Path(file_path).name)[0]
        with open(os.path.join(pathlib.Path(file_path).parent, filename + '_Para.txt'), 'w') as f:
            for line in Parlist:
                f.write("%s\n" % line) 
        parFile.close()

    return parameters

def readKSpace(file_path, scanParams = {}, correctOversampling = True, forceOversamplingCorrection = False):
    fid = open(file_path, 'rb')
    nmrsamples = np.double(np.fromfile(fid, np.int32, 1))
    nmrviews = np.double(np.fromfile(fid, np.int32, 1))
    nmrviews2 = np.fromfile(fid, np.int32, 1)
    nmrslices = np.fromfile(fid, np.int32, 1)
    fid.seek(2, 1)
    nmrdatatype = np.fromfile(fid, np.int16, 1)
    fid.seek(152, 0)
    nmrechoes = np.fromfile(fid, np.int32, 1)
    nmrexperiments = np.fromfile(fid, np.int32, 1)
    nfids = np.double(nmrviews*nmrviews2*nmrslices*nmrechoes*nmrexperiments)

    ####   read data   ###

    fid.seek(512, 0)
    AR = np.fromfile(fid, np.float32, int(nmrsamples*2*nfids),'')
    fid.seek(int(-8*nmrsamples*nfids + 4), 1)
    AI = np.fromfile(fid, np.float32, int(nmrsamples*2*nfids),'')
    imagMatrix = AI[::2]
    imagMatrix = np.reshape(imagMatrix, (int(nfids), int(nmrsamples)))
    imagMatrix = np.double(imagMatrix)
    realMatrix = AR[::2]
    realMatrix = np.reshape(realMatrix, (int(nfids), int(nmrsamples)))
    realMatrix = np.double(realMatrix)
    
    nmrdata = realMatrix + imagMatrix*1j
    nmrdata = nmrdata.T 

    if 'NO_VIEWS_2' in scanParams:
        ###   Trajectory table   ####
        no_views = int(nmrviews)
        array_count =0
        echo_table = np.zeros(no_views)
        while array_count < no_views/2:
            if array_count == 0:
                echo_table[array_count] = array_count
            
            if array_count >= 1 :
                echo_table[2*(array_count-1)+1] = array_count
                echo_table[2*(array_count-1)+2] = -array_count
            echo_table[int(no_views)-1] = (no_views/2)
            array_count+=1
            
        current_view=0
        C=[]
        views_per_seg = int(scanParams['VIEWS_PER_SEGMENT'][1])
        while current_view<no_views:
            echo_cnt = 0
            while echo_cnt < views_per_seg:
                gp_mul = echo_table[int(echo_cnt+(current_view/views_per_seg)+(echo_cnt*(no_views/views_per_seg -1)))]
                C.append(gp_mul)
                echo_cnt+=1
            current_view +=views_per_seg
        table=[]
        i=0
        while i <no_views:
            table.append(int(C[i]+ no_views/2))
            i +=1

        ###  reordering Kspace data      ###
        ETL = int(scanParams['VIEWS_PER_SEGMENT'][1])
        if int(nmrviews2)>>1 :
            MessedUpdata = np.zeros((int(nmrsamples), int(nmrviews), int(nmrviews2)), dtype =complex)
            shot = 0
            for kk in range(0, int(nmrviews/ETL)):
                for ii in range(0, int(nmrviews2)):
                    shot_window = slice(shot*ETL, (shot+1)*ETL)
                    phase1_window = slice(ETL*kk, ETL*(kk+1))
                    MessedUpdata[:, phase1_window, ii] = nmrdata[:, shot_window]
                    shot += 1

            sortdata = np.zeros((int(nmrsamples), int(nmrviews), int(nmrviews2)), dtype =complex)
            for ii in range(int(nmrviews2)):
                temp1 = np.squeeze(MessedUpdata[:,:,ii])
                for kk in range(len(table)):
                    sortdata[table[kk]-1,:,ii] = temp1[:,kk]
        else:
            sortdata = np.zeros((int(nmrsamples), int(nmrviews)), dtype =complex)
            for kk in range(len(table)):
                sortdata[table[kk]-1,:] = nmrdata[:,kk]
    else:
        sortdata = nmrdata.T
    return sortdata
