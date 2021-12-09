font = {'family' : 'sans-serif',
        'color' : 'k',
        'style' : 'normal', 
        'variant' : 'normal',
        'weight' : 'normal',
        'size' : 'medium'}


import tensorflow as tf
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)

    except RuntimeError as e:
        print(e)

import keras
import os
import pickle
import numpy as np
import statsmodels

# import matplotlib
# import matplotlib.pyplot as plt
# from matplotlib import cm

import seaborn as sns
import math

import pandas as pd

import scipy

from matplotlib.colors import LogNorm, LinearSegmentedColormap

import matplotlib.patches as mpatches
import matplotlib.lines as mlines

from sklearn.cluster import KMeans, MeanShift
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import KernelDensity



def getCoord(x, y):

    coordXY = []

    for trackX, trackY in zip(x, y):
        coxy = []

        for j in range(len(trackX) - 1):
            coxy.append([trackX[j], trackY[j], trackX[j + 1], trackY[j + 1]])
        coordXY.append(coxy)

    return coordXY

def getMeanMSDs(x, y, maxOrder, shift):

    D= []
    for trackX, trackY in zip(x, y):
        dTrack = []
        for order in range(1, maxOrder + 1):
            dOrder = []
            for j in range(len(trackX) - order):
                dOrder.append(np.sqrt((trackX[j + order] - trackX[j]) ** 2 + (trackY[j + order] - trackY[j]) ** 2))
            dTrack.append(dOrder)
        D.append(dTrack)

    dVec = []
    for track in D:
        vector = []
        for order in range(1, maxOrder + 1):
            v = []
            dO = track[order - 1]

            for i in range(len(track[0])):
                if order % 2 == 0:
                    iStart = (i - shift - ((order - 1) // 2))
                    iEnd = iStart + (2 * shift) - 1
                else:
                    iStart = (i - shift - (order // 2))
                    iEnd = iStart + (2 * shift)

                if iStart < 0:
                    iStart = 0

                if iEnd > (len(track[0]) - 1):
                    iEnd = len(track[0]) - 1

                meanPoint = np.mean(dO[iStart:(iEnd + 1)])
                v.append(meanPoint)
            vector.append(v)
        dVec.append(vector)

    return dVec


def getFeatVec(d, x, y, addFeat, maxOrder, shift):

    if 'meanMSD' in addFeat:
        dMeanVec = getMeanMSDs(x, y, maxOrder, shift)

    if 'xy' in addFeat:
        xy = getCoord(x, y)

    featVec = []
    for i in range(len(d)):
        featTrack = []
        for j in range(len(d[i])):
            featPoint = []
            featPoint.append(d[i][j])

            if 'meanMSD' in addFeat:
                for order in range(maxOrder):
                    featPoint.append(dMeanVec[i][order][j])

            if 'xy' in addFeat:
                featPoint.append(xy[i][j][0])
                featPoint.append(xy[i][j][1])
                featPoint.append(xy[i][j][2])
                featPoint.append(xy[i][j][3])

            featTrack.append(featPoint)
        featVec.append(featTrack)

    return featVec

def getDist(x, y, deltaT = 1):

    distances = []

    for trackX, trackY in zip(x, y):
        distances.append(np.sqrt((np.array(trackX)[deltaT:] - np.array(trackX)[:-deltaT]) ** 2 + \
                                 (np.array(trackY)[deltaT:] - np.array(trackY)[:-deltaT]) ** 2))

    return distances


def loadRealData(filename):
    xCoord = []
    yCoord = []
    trInd = 0
    pos = 0

    while trInd <= filename[-1, 3]:
        if trInd == filename[-1, 3]:
            nextt = len(filename)
        else:
            nextt = list(filename[:,3]).index(trInd + 1)

        if nextt - pos > 3:
            frames = filename[pos:nextt, 0]
            xx = filename[pos:nextt, 1]
            yy = filename[pos:nextt, 2]

            missing = list(set(np.arange(frames[0], frames[-1] + 1)) - set(frames))
            if len(missing) != 0:
                missing = sorted(missing)
                xxnew = xx
                yynew = yy
                for i in missing:
                    xxnew = np.insert(xxnew, int(i - frames[0]), np.interp(i, frames, xx))
                    yynew = np.insert(yynew, int(i - frames[0]), np.interp(i, frames, yy))
                xCoord.append(xxnew)
                yCoord.append(yynew)
            else:
                xCoord.append(xx)
                yCoord.append(yy)

        pos = nextt
        trInd += 1

    return xCoord, yCoord

def getMSDandMSS(x, y, numPmsd, numPmss, p, b = 'unknown'):

    # Compute MSD
    dts = np.arange(1, numPmsd + 1)
    timemean = dts - dts.mean()
    sumtsq = sum(timemean ** 2)

    MSD = []
    for dt in dts:
        ddt = getDist(x, y, deltaT = dt)
        flatD = np.asarray([di ** 2 for Ds in ddt for di in Ds])
        MSD.append(flatD.mean())

    if b == 'zero':
        diff = np.array(MSD).mean() / dts.mean() * pixSize ** 2 / t / 4
    else:
        diff = sum(timemean * (MSD - np.array(MSD).mean())) / sumtsq * pixSize ** 2 / t / 4

    # Compute MSS
    dts = np.arange(1, numPmss + 1)
    MSS = []

    for pi in p:
        logM = []
        logT = []
        for dt in dts:
            ddt = getDist(x, y, deltaT = dt)
            flatD = np.asarray([di for Ds in ddt for di in Ds])
            logM.append(np.log10((flatD ** pi).mean()))
            logT.append(np.log10(dt))
        logTmean = np.array(logT) - np.array(logT).mean()
        sumLogTsq = sum(logTmean ** 2)
        logMmean = np.array(logM) - np.array(logM).mean()
        MSS.append(sum(logTmean * logMmean) / sumLogTsq)

    # Compute Smss and intercept
    pmean = p - p.mean()
    sumpsq = sum(pmean ** 2)
    MSSmean = np.array(MSS) - np.array(MSS).mean()
    Smss = sum(pmean * MSSmean) / sumpsq
    intercept = np.array(MSS).mean() - Smss * p.mean()

    return diff, MSS, Smss, intercept

def getMSDandMSS3D(x, y, z, numPmsd, numPmss, p, b = 'unknown'):

    # Compute MSD
    dts = np.arange(1, numPmsd + 1)
    timemean = dts - dts.mean()
    sumtsq = sum(timemean ** 2)

    MSD = []
    for dt in dts:
        ddt = getDist3D(x, y, z, deltaT = dt)
        flatD = np.asarray([di ** 2 for Ds in ddt for di in Ds])
        MSD.append(flatD.mean())

    if b == 'zero':
        diff = np.array(MSD).mean() / dts.mean() * pixSize ** 2 / t /6
    else:
        diff = sum(timemean * (MSD - np.array(MSD).mean())) / sumtsq * pixSize ** 2 / t / 6

    # Compute MSS
    dts = np.arange(1, numPmss + 1)
    MSS = []

    for pi in p:
        logM = []
        logT = []
        for dt in dts:
            ddt = getDist3D(x, y, deltaT = dt)
            flatD = np.asarray([di for Ds in ddt for di in Ds])
            logM.append(np.log10((flatD ** pi).mean()))
            logT.append(np.log10(dt))
        logTmean = np.array(logT) - np.array(logT).mean()
        sumLogTsq = sum(logTmean ** 2)
        logMmean = np.array(logM) - np.array(logM).mean()
        MSS.append(sum(logTmean * logMmean) / sumLogTsq)

    # Compute Smss and intercept
    pmean = p - p.mean()
    sumpsq = sum(pmean ** 2)
    MSSmean = np.array(MSS) - np.array(MSS).mean()
    Smss = sum(pmean * MSSmean) / sumpsq
    intercept = np.array(MSS).mean() - Smss * p.mean()

    return diff, MSS, Smss, intercept


def getMSDandMSSandC(x, y, numPmsd, numPmss, p, b = 'unknown'):

    # Compute MSD
    dts = np.arange(1, numPmsd + 1)
    timemean = dts - dts.mean()
    sumtsq = sum(timemean ** 2)

    MSD = []
    for dt in dts:
        ddt = getDist(x, y, deltaT = dt)
        flatD = np.asarray([di ** 2 for Ds in ddt for di in Ds])
        MSD.append(flatD.mean())

    if b == 'zero':
        diff = np.array(MSD).mean() / dts.mean() * pixSize ** 2 / t / 4
    else:
        diff = sum(timemean * (MSD - np.array(MSD).mean())) / sumtsq * pixSize ** 2 / t / 4

    # Compute MSS
    dts = np.arange(1, numPmss + 1)
    MSS = []
    C = []
    cD = []

    for pi in p:
        logM = []
        logT = []
        for dt in dts:
            ddt = getDist(x, y, deltaT = dt)
            flatD = np.asarray([di for Ds in ddt for di in Ds])
            logM.append(np.log((flatD ** pi).mean()))
            logT.append(np.log(dt))
        logTmean = np.array(logT) - np.array(logT).mean()
        sumLogTsq = sum(logTmean ** 2)
        logMmean = np.array(logM) - np.array(logM).mean()
        gamma =  sum(logTmean * logMmean) / sumLogTsq
        MSS.append(gamma)
        C.append(np.array(logM).mean() - gamma * np.array(logT).mean())
        cD.append((np.exp(np.array(logM).mean()) / 2 ** (pi / 2) / scipy.special.gamma(1 + pi / 2)) ** (2 / pi) \
                  / 2 / np.exp(np.array(logT).mean()))

    # Compute Smss and intercept
    pmean = p - p.mean()
    sumpsq = sum(pmean ** 2)
    MSSmean = np.array(MSS) - np.array(MSS).mean()
    Smss = sum(pmean * MSSmean) / sumpsq
    intercept = np.array(MSS).mean() - Smss * p.mean()

    return diff, MSS, C, cD, Smss, intercept

def getTrackPieces(x, y, allStates):

    piecesX0 = []
    piecesY0 = []
    piecesX1 = []
    piecesY1 = []
    piecesX2 = []
    piecesY2 = []

    for trX, trY, trS in zip(x, y, allStates):
        sw = 0
        pos = 0
        swList = []
        while list(trS)[sw:].count(trS[sw]) != len(trS[sw:]):
            sw = next(ind for ind in range(sw + 1, len(trS)) if trS[ind] != trS[sw])
            if trS[sw - 1] == 0:
                piecesX0.append(trX[pos:sw])
                piecesY0.append(trY[pos:sw])
            elif trS[sw - 1] == 1:
                piecesX1.append(trX[pos:sw])
                piecesY1.append(trY[pos:sw])
            else:
                piecesX2.append(trX[pos:sw])
                piecesY2.append(trY[pos:sw])
            pos = sw
        else:
            if trS[pos] == 0:
                piecesX0.append(trX[pos:])
                piecesY0.append(trY[pos:])
            elif trS[pos] == 1:
                piecesX1.append(trX[pos:])
                piecesY1.append(trY[pos:])
            else:
                piecesX2.append(trX[pos:])
                piecesY2.append(trY[pos:])

    return piecesX0, piecesY0, piecesX1, piecesY1, piecesX2, piecesY2

def getTrackPiecesForInfo(x, y, allStates):

    piecesX0 = []
    piecesY0 = []
    piecesX1 = []
    piecesY1 = []
    piecesX2 = []
    piecesY2 = []
    trnum = 0
    trnums0 = []
    trnums1 = []
    trnums2 = []

    for trX, trY, trS in zip(x, y, allStates):
        sw = 0
        pos = 0
        swList = []
        n = 1
        while list(trS)[sw:].count(trS[sw]) != len(trS[sw:]):
            sw = next(ind for ind in range(sw + 1, len(trS)) if trS[ind] != trS[sw])
            if trS[sw - 1] == 0:
                piecesX0.append(trX[pos:sw])
                piecesY0.append(trY[pos:sw])
                trnums0.append(trnum + n/10)
            elif trS[sw - 1] == 1:
                piecesX1.append(trX[pos:sw])
                piecesY1.append(trY[pos:sw])
                trnums1.append(trnum + n/10)
            else:
                piecesX2.append(trX[pos:sw])
                piecesY2.append(trY[pos:sw])
                trnums2.append(trnum + n/10)
            pos = sw
            n += 1
        else:
            if trS[pos] == 0:
                piecesX0.append(trX[pos:])
                piecesY0.append(trY[pos:])
                trnums0.append(trnum + n/10)
            elif trS[pos] == 1:
                piecesX1.append(trX[pos:])
                piecesY1.append(trY[pos:])
                trnums1.append(trnum + n/10)
            else:
                piecesX2.append(trX[pos:])
                piecesY2.append(trY[pos:])
                trnums2.append(trnum + n/10)
        trnum += 1

    return piecesX0, piecesY0, piecesX1, piecesY1, piecesX2, piecesY2, trnums0, trnums1, trnums2

def patchmatch(match, pattern):
    indices = []
    for pat in pattern:
        found = False
        patlen = len(pat)
        for i in range(len(match)):
            rollAr = np.roll(match, -i)
            if rollAr[:patlen].tolist() == pat:
                ind = np.arange(i, i + patlen)
                found = True
                break
        if found == False:
            print('No match for %s' %(pat))
        else:
            indices.append(ind.tolist())

    return indices
    
def releaseMemory():
    from numba import cuda
    cuda.select_device(0)
    cuda.close()
    
