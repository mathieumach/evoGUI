# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:33:25 2021

@author: Low field
"""
import numpy as np
from pathlib import Path

def gaussianFilter(rawData, p1_param = 1/2, p2_param = 1/6):
    inputShape = np.shape(rawData)
    filterMat = 1
    for dimSize in inputShape:
        N = dimSize
        p1 = N*p1_param
        p2 = N*p2_param
        filterVec = np.exp(-(np.square(np.arange(N) - p1)/(p2**2)))
        filterMat = np.multiply.outer(filterMat, filterVec)
    return np.multiply(rawData, filterMat)

def sineBellSquaredFilter(complexData,  filterStrength = 1):
    inputShape = np.shape(complexData)
    filterMat = 1
    for dimSize in inputShape:
        N = dimSize
        p1 = N/2
        axis = np.linspace(-N/2, N/2, N)
        filterVec = 1-filterStrength*np.square(np.cos(0.5*np.pi*(axis-p1)/(N-p1)))
        filterMat = np.multiply.outer(filterMat, filterVec)
    return np.multiply(complexData, filterMat)


def shiftImage(kSpace,scanParams, distance, axis):
    """ Shift image by adding a time dependent phase to one of the components,
        Distance defines how far in mm, axis the axis along which it should be moved"""
    if axis == 0:
        fov = float(scanParams['FOVread'])
        nrPts = int(scanParams['nrPnts']) 
    elif axis == 1: 
        fov = float(scanParams['FOVphase1'])
        nrPts = int(scanParams['nPhase1']) 
    elif axis == 2:
        fov = float(scanParams['FOVphase2'])
        nrPts = int(scanParams['nPhase2']) 
    else:
        raise("Invalid axis")
    
    hzPerMM = 1e3*float(scanParams['bandwidth'])/float(scanParams['FOVread'])
    freqShift = distance * hzPerMM
    acqTime = nrPts/float(scanParams['bandwidth'])
    timeScale = np.linspace(-acqTime/2, acqTime/2, nrPts)
    phaseTerm = np.exp(-1j*2*np.pi*timeScale*freqShift)
    
    if axis == 0:
        shiftedKspace = np.multiply(kSpace, phaseTerm[:,np.newaxis, np.newaxis])
    elif axis == 1: 
        shiftedKspace = np.multiply(kSpace, phaseTerm[np.newaxis,:, np.newaxis])
    elif axis == 2:
        shiftedKspace = np.multiply(kSpace, phaseTerm[np.newaxis, np.newaxis,:])
    return shiftedKspace

def noiseCorrection(kSpace, scanParams, dataFolder):
    try:
        '''store descriptor in parent folder to be accessed by other reconstruction methods'''
        dataPath            = Path(dataFolder)
        mainDir             = str(dataPath.parent)
        noiseDict           = np.load(r"%s/NoiseData.npy"%(mainDir), allow_pickle = True).item()
    except:
        print("Could not load noise scan data")
        return kSpace
    
    freqOffset = float(scanParams["b1Freq"][:-1]) - float(noiseDict["center_freq"][:-1])
    
    print("Freq difference: %.0f Hz"%(freqOffset*1e6))
    
    imageBandwidth  = np.linspace(-float(scanParams["reconBW"])/2,float(scanParams["reconBW"])/2, np.size(kSpace, 0))
    
    polyFunc        = np.poly1d(noiseDict["noise_fit"])
    noiseFit        = polyFunc(imageBandwidth+freqOffset)
    
    image           = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(kSpace)))

    if np.size(np.shape(kSpace)) == 2 or np.shape(kSpace)[2] == 1:
        image           /= noiseFit[:,np.newaxis]
    else:
        image           /= noiseFit[:,np.newaxis,np.newaxis]
    
    return np.fft.fftshift(np.fft.ifftn(np.fft.fftshift(image)))

def zeroFill(data, zeroFillDimensions):
    if np.size(np.shape(data)) is not np.size(zeroFillDimensions):
        print("Dimensions of input array must match dimensions of output array")
        return -1
        
    centerPoint = np.array(np.floor(np.divide(zeroFillDimensions,2)), dtype = int)
    
    dataSizeMin = np.array(np.floor(np.divide(np.shape(data),2)), dtype = int)
    dataSizeMax = np.array(np.ceil(np.divide(np.shape(data),2)), dtype = int)
    
    zeroFilledData = np.zeros(zeroFillDimensions, dtype = np.complex)
    zeroFilledData[centerPoint[0] - dataSizeMin[0]: centerPoint[0] + dataSizeMax[0],\
                   centerPoint[1] - dataSizeMin[1]: centerPoint[1] + dataSizeMax[1],\
                   centerPoint[2] - dataSizeMin[2]: centerPoint[2] + dataSizeMax[2]] = data
    return zeroFilledData
    
    