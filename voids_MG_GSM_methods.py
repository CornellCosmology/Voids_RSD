#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 11:56:47 2022

@author: christopherwilson
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 10:22:17 2022

@author: christopherwilson
"""
import numpy as np
import time
import copy
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import scipy as sp
import scipy.linalg
from scipy import integrate
from scipy import special 
from scipy import misc
from scipy import special 
#from mpmath import *
from scipy import interpolate
from numpy import polynomial 
from scipy.interpolate import interp1d
#from forceMethods import *
from scipy.misc import derivative as diff


from numpy import nan

import seaborn as sns
import cmath as cmath 
import math
from pylab import *
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from scipy.special import binom as binomial


z=0.48251
a=1./(1.+z)

omega_M=0.3089
omega_L=1.-omega_M
H=100.*np.sqrt(omega_L+omega_M*(1./(a**3)))

sqrt2pi=np.sqrt(2*np.pi)


"##############################################################################"
voidDir="/Users/christopherwilson/Desktop/Research/GLAM/125/catalogs/"   
"#### REPLACE THIS WITH YOUR LOCAL PATH TO THE CATALOG DIRECTORY ####"




newCalcDir="/Users/christopherwilson/Desktop/Research/GLAM/125/newCalculations/"
"#### Directory where calculations may be saved to speed up some methods. Not required that you change or use this ####"




"#########################################################################"
class catalogClass:
    
    "#########################################################################"
    def __init__(self,voidInfos=None,partsRandMu=None,densVelsSigmas=None,velZArrays=None,scaling="",centerType="CC",gravType="GR"):
    
    
        self.scaling=scaling
        self.centerType=centerType
        self.gravType=gravType
        loadMuProfs=True
        loadDensVelsSigmas=True
        loadVelZArrays=True
    
        if partsRandMu==None:
            loadMuProfs=False
            # print("NOT loading redshift space R,mu parts profile")
            
        if densVelsSigmas==None:
            loadDensVelsSigmas=False
            print("NOT loading densities, velocities, sigmas")
    
        if velZArrays==None:
            loadVelZArrays=False
            
        
        if loadVelZArrays==True:
            print("loading velZArrays")
    
    
    
        if scaling=="Reff":
            xpoints=(np.linspace(0,3,31)[1:] + np.linspace(0,3,31)[:-1])/2
            self.xpoints=xpoints
            self.muPoints=(np.linspace(0,1,51)[1:] + np.linspace(0,1,51)[:-1])/2
        if scaling=="R24":
            xpoints=(np.linspace(0,5,51)[1:] + np.linspace(0,5,51)[:-1])/2*24.
            self.xpoints=xpoints
            self.muPoints=(np.linspace(0,1,51)[1:] + np.linspace(0,1,51)[:-1])/2
        
        if voidInfos==None:
            voidList=None
        

        if voidInfos != None:
        
            q=0
            
            for i in range(len(voidInfos)):
                for j in range(len(voidInfos[i])):
                    q+=1
            voidList=[[] for i in range(q)] # DONT MAKE THIS A NUMPY ARRAY
            qq=0
            
            A=None
            B=[None,None,None]
            C=None
            
            try:
                w=len(velZArrays)
            except:
                w=0
            "w tells us how many realizations for velZArrays have been loaded in, assuming we always load in from real 1 first"
            
            self.maxRealForVelArray=w


            for i in range(len(voidInfos)):
                for j in range(len(voidInfos[i])):
                    if loadVelZArrays==True:
                        if i<w:
                            C=velZArrays[i][j]
                        if i>=w:
                            C=None
                    
                    
                    
                    
                    if loadMuProfs==True:
                        A=partsRandMu[i][j]
                    if loadDensVelsSigmas==True:
                        B=densVelsSigmas[i][j]
                        
                    voidList[qq]=voidClass(voidInfos[i][j],partsRandMu=A,densVelsSigmas=B,velZArray=C,scaling=scaling,xpoints=xpoints)
                    qq+=1
                
            self.voids=voidList
            self.numVoids=len(voidList)
    "#########################################################################"   
    
    "#########################################################################"
    def append(self,voidsToAppend):
        voidsTotal=np.array( ( list( self.voids ) ) + ( list( voidsToAppend ) ) )
        self.voids=voidsTotal
        self.numVoids=len(self.voids)
    "#########################################################################"
    
    "#########################################################################"
    def filterVoids(self,rMin,rMax,voidTypes=["R","S"],realsToUse=range(1,101)):
        # [v for void in self.voids if rMin<=v.radius<=rMax and v.voidType in voidTypes]
        self.voidsToUse=[v for v in self.voids if rMin<=v.radius<=rMax and v.voidType in voidTypes and v.realNum in realsToUse]
    "#########################################################################"    
    
    "#########################################################################"
    def calcQuadMoment(self,rMin,rMax,voidTypes=["R","S"],realsToUse=range(1,101)):
        
        self.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes,realsToUse=realsToUse)
        
        partsArray=np.zeros(np.shape(self.voids[0].partsRandMu))
        nMuBins=len(partsArray[0]) # should be 50 
        
        muBinEdges=np.linspace(0,1,nMuBins+1)
        dMu=muBinEdges[1]-muBinEdges[0]
        muPoints=(muBinEdges[0:-1]+muBinEdges[1:])/2
        
        
        
        totalShellVols=np.zeros(len(self.voids[0].shellVols))
        nBar=0
        for v in self.voidsToUse:
            partsArray+=v.partsRandMu
            totalShellVols+=v.shellVols
            nBar+=v.nBar
        nBar=nBar/len(self.voidsToUse)
        densArray=np.array([partsArray[i]/(totalShellVols[i]/nMuBins)/nBar - 1. for i in range(len(partsArray))])
        self.densArray=densArray
        quadMoment=np.array([np.sum([densArray[i][j]*5/2*(3*muPoints[j]**2-1)*dMu for j in range(len(partsArray[i]))]) for i in range(len(partsArray))])
        
        # quad=np.zeros(30)
        # for i in range(len(quad)):
        #     quad[i]=np.sum([partsArray[i][j]*5/2*(3*muPoints[j]**2-1) for j in range(len(partsArray[i]))])/nBar/totalShellVols[i]
        self.quadMoment=quadMoment
        self.rMinQuad=rMin
        self.rMaxQuad=rMax
        
        
        
        
        # self.quad=quad
        return np.array([self.xpoints,quadMoment])
    "#########################################################################"
    
    "#########################################################################"    
    def jackKnifeQuadMoment(self,rMin=0,rMax=100,voidTypes=["R","S"],numSamples=1000,realsToUse=range(1,101)):
        """
        Using the jackknife method, computes the mean quadrupole, and the standard uncertainty on the mean quadrupole 
        according to all available data provided. To get the standard uncertainty on the mean of the void quadrupole 
        from a (1 Mpc/h)^3 volume of space, the standard uncertainty on the mean according to a volume of X (Mpc/h)^3
        must be scaled by \sqrt(X)
        """
        
        
        
        self.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes,realsToUse=realsToUse)
        partsArraySamples=np.zeros((numSamples,np.shape(self.voids[0].partsRandMu)[0],np.shape(self.voids[0].partsRandMu)[1]))
        shellVolArraySamples=np.zeros((numSamples,len(self.voids[0].shellVols)))
        nBarSamples=np.zeros(numSamples)
        numVoidsSamples=np.zeros(numSamples)
    
        
        
        partsArrayTotal=np.zeros(np.shape(self.voids[0].partsRandMu))
        shellVolTotal=np.zeros(len(self.voids[0].shellVols))
        
        
        
        nMuBins=len(partsArrayTotal[0]) # should be 50 
        
        muBinEdges=np.linspace(0,1,nMuBins+1)
        dMu=muBinEdges[1]-muBinEdges[0]
        muPoints=(muBinEdges[0:-1]+muBinEdges[1:])/2
        
        
        
        

        for i in range(numSamples):
            voidsInSample=[self.voidsToUse[j] for j in range(len(self.voidsToUse)) if j%numSamples == i]
            for v in voidsInSample:
                partsArraySamples[i] += v.partsRandMu
                shellVolArraySamples[i] += v.shellVols
                nBarSamples[i] += v.nBar
                numVoidsSamples[i] += 1
                    
        
        
        nBarTotal=np.sum(nBarSamples)/np.sum(numVoidsSamples)
        
        for i in range(numSamples):
            partsArrayTotal += partsArraySamples[i]
            shellVolTotal += shellVolArraySamples[i]
            

        moments=np.array([np.zeros(len(self.xpoints)) for i in range(numSamples)])
        for i in range(numSamples):
            
            nBarTemp=(np.sum(nBarSamples) - nBarSamples[i])/(np.sum(numVoidsSamples)-numVoidsSamples[i])
            partsTemp = partsArrayTotal-partsArraySamples[i]
            shellVolsTemp = shellVolTotal-shellVolArraySamples[i]            
            
            densArrayTemp=np.array([partsTemp[k]/(shellVolsTemp[k]/nMuBins)/nBarTemp - 1. for k in range(len(partsTemp))])
            
            moments[i]=np.array([np.sum([densArrayTemp[k][j]*5/2*(3*muPoints[j]**2-1)*dMu for j in range(len(partsTemp[k]))]) for k in range(len(partsTemp))])
            
        densArrayTotal=np.array([partsArrayTotal[k]/(shellVolTotal[k]/nMuBins)/nBarTotal - 1. for k in range(len(partsArrayTotal))])
        meanMoment=np.array([np.sum([densArrayTotal[k][j]*5/2*(3*muPoints[j]**2-1)*dMu for j in range(len(partsArrayTotal[k]))]) for k in range(len(partsArrayTotal))])
                    
        varSq=np.zeros(len(moments[0]))
        for i in range(len(moments)):
            varSq+=(len(moments)-1)/len(moments)*(moments[i]-meanMoment)**2
            
        var=np.sqrt(varSq)
        
        
        
        self.rMinQuad=rMin
        self.rMaxQuad=rMax
        
        
        
        self.quadErrors=var
        self.quadMoment=meanMoment
        
        self.densArrayTotal=densArrayTotal
        
        
        return np.array([self.xpoints,meanMoment,var])
    "#########################################################################"
    
    "#########################################################################"    
    def jackKnifeXi4Moment(self,rMin=0,rMax=100,voidTypes=["R","S"],numSamples=1000,realsToUse=range(1,101)):
        
        self.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes,realsToUse=realsToUse)
        partsArraySamples=np.zeros((numSamples,np.shape(self.voids[0].partsRandMu)[0],np.shape(self.voids[0].partsRandMu)[1]))
        shellVolArraySamples=np.zeros((numSamples,len(self.voids[0].shellVols)))
        nBarSamples=np.zeros(numSamples)
        numVoidsSamples=np.zeros(numSamples)
    
        
        partsArrayTotal=np.zeros(np.shape(self.voids[0].partsRandMu))
        shellVolTotal=np.zeros(len(self.voids[0].shellVols))
        
        
        
        nMuBins=len(partsArrayTotal[0]) # should be 50 
        
        muBinEdges=np.linspace(0,1,nMuBins+1)
        dMu=muBinEdges[1]-muBinEdges[0]
        muPoints=(muBinEdges[0:-1]+muBinEdges[1:])/2
        
        

        for i in range(numSamples):
            voidsInSample=[self.voidsToUse[j] for j in range(len(self.voidsToUse)) if j%numSamples == i]
            for v in voidsInSample:
                partsArraySamples[i] += v.partsRandMu
                shellVolArraySamples[i] += v.shellVols
                nBarSamples[i] += v.nBar
                numVoidsSamples[i] += 1
                    
        
        
        nBarTotal=np.sum(nBarSamples)/np.sum(numVoidsSamples)
        
        for i in range(numSamples):
            partsArrayTotal += partsArraySamples[i]
            shellVolTotal += shellVolArraySamples[i]
            

        moments=np.array([np.zeros(len(self.xpoints)) for i in range(numSamples)])
        for i in range(numSamples):
            
            nBarTemp=(np.sum(nBarSamples) - nBarSamples[i])/(np.sum(numVoidsSamples)-numVoidsSamples[i])
            partsTemp = partsArrayTotal-partsArraySamples[i]
            shellVolsTemp = shellVolTotal-shellVolArraySamples[i]            
            
            densArrayTemp=np.array([partsTemp[k]/(shellVolsTemp[k]/nMuBins)/nBarTemp - 1. for k in range(len(partsTemp))])
            
            moments[i]=np.array([np.sum([densArrayTemp[k][j]*(2*(4) + 1)/2*(1/8)*(2)*(35*muPoints[j]**4 - 30*muPoints[j]**2 + 3)*dMu for j in range(len(partsTemp[k]))]) for k in range(len(partsTemp))])
            
            # meanMoment=np.array([np.sum([densArrayTotal[k][j]*5/2*(3*muPoints[j]**2-1)*dMu for j in range(len(partsArrayTotal[k]))]) for k in range(len(partsArrayTotal))])
            
            
            
            
            
            
            
        densArrayTotal=np.array([partsArrayTotal[k]/(shellVolTotal[k]/nMuBins)/nBarTotal - 1. for k in range(len(partsArrayTotal))])
        meanMoment=np.array([np.sum([densArrayTotal[k][j]*(2*(4) + 1)/2*(1/8)*(2)*(35*muPoints[j]**4 - 30*muPoints[j]**2 + 3)*dMu for j in range(len(partsArrayTotal[k]))]) for k in range(len(partsArrayTotal))])
                    
    

        
        
        varSq=np.zeros(len(moments[0]))
        for i in range(len(moments)):
            varSq+=(len(moments)-1)/len(moments)*(moments[i]-meanMoment)**2
            
        var=np.sqrt(varSq)
        
        
        
        self.rMinXi4=rMin
        self.rMaxXi4=rMax
        
        
        
        self.Xi4Errors=var
        self.Xi4Moment=meanMoment
        
        self.densArrayTotal=densArrayTotal
        
        
        return np.array([self.xpoints,meanMoment,var])
    "#########################################################################"
    
    "#########################################################################"    
    def jackKnifeMonopoleMoment(self,rMin=0,rMax=100,voidTypes=["R","S"],numSamples=1000,realsToUse=range(1,101)):
        
        self.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes,realsToUse=realsToUse)
        partsArraySamples=np.zeros((numSamples,np.shape(self.voids[0].partsRandMu)[0],np.shape(self.voids[0].partsRandMu)[1]))
        shellVolArraySamples=np.zeros((numSamples,len(self.voids[0].shellVols)))
        nBarSamples=np.zeros(numSamples)
        numVoidsSamples=np.zeros(numSamples)
    
        
        
        partsArrayTotal=np.zeros(np.shape(self.voids[0].partsRandMu))
        shellVolTotal=np.zeros(len(self.voids[0].shellVols))
        
        
        
        nMuBins=len(partsArrayTotal[0]) # should be 50 
        
        muBinEdges=np.linspace(0,1,nMuBins+1)
        dMu=muBinEdges[1]-muBinEdges[0]
        muPoints=(muBinEdges[0:-1]+muBinEdges[1:])/2
        
        
        
        
        
        for i in range(numSamples):
            voidsInSample=[self.voidsToUse[j] for j in range(len(self.voidsToUse)) if j%numSamples == i]
            for v in voidsInSample:
                partsArraySamples[i] += v.partsRandMu
                shellVolArraySamples[i] += v.shellVols
                nBarSamples[i] += v.nBar
                numVoidsSamples[i] += 1
                    
        
        
        nBarTotal=np.sum(nBarSamples)/np.sum(numVoidsSamples)
        
        for i in range(numSamples):
            partsArrayTotal += partsArraySamples[i]
            shellVolTotal += shellVolArraySamples[i]
            

        moments=np.array([np.zeros(len(self.xpoints)) for i in range(numSamples)])
        for i in range(numSamples):
            
            nBarTemp=(np.sum(nBarSamples) - nBarSamples[i])/(np.sum(numVoidsSamples)-numVoidsSamples[i])
            partsTemp = partsArrayTotal-partsArraySamples[i]
            shellVolsTemp = shellVolTotal-shellVolArraySamples[i]            
            
            densArrayTemp=np.array([partsTemp[k]/(shellVolsTemp[k]/nMuBins)/nBarTemp - 1. for k in range(len(partsTemp))])
            
            moments[i]=np.array([np.sum([densArrayTemp[k][j]*1*dMu for j in range(len(partsTemp[k]))]) for k in range(len(partsTemp))])
            
        densArrayTotal=np.array([partsArrayTotal[k]/(shellVolTotal[k]/nMuBins)/nBarTotal - 1. for k in range(len(partsArrayTotal))])
        meanMoment=np.array([np.sum([densArrayTotal[k][j]*1*dMu for j in range(len(partsArrayTotal[k]))]) for k in range(len(partsArrayTotal))])
                    
        varSq=np.zeros(len(moments[0]))
        for i in range(len(moments)):
            varSq+=(len(moments)-1)/len(moments)*(moments[i]-meanMoment)**2
            
        var=np.sqrt(varSq)
        
        
        
        self.rMinMono=rMin
        self.rMaxMono=rMax
        
        
        
        self.monoErrors=var
        self.monoMoment=meanMoment
        
        
        
        return np.array([self.xpoints,meanMoment,var])
    "#########################################################################"
       
    "#########################################################################"
    def calcVelocityFromBeta(self,beta,v=0,rMin=0,rMax=100,voidTypes=["R","S"],plotting=False,shiftFac=1):
        dx=(self.xpoints[1]-self.xpoints[0])/2
        

        Delta=self.calcIntegratedDensityProfile(rMin=rMin,rMax=rMax,voidTypes=["R","S"])
        
        if self.scaling=="Reff":
            weight=-1
        if self.scaling=="R24":
            weight=0
        
        
        
        betaVelInterp=sp.interpolate.CubicSpline((Delta[0] + dx*(1.+shiftFac)),-1/3*a*H*beta*Delta[1]*(Delta[0]+dx),bc_type = 'natural')
        
        betaVel=np.array([self.xpoints,betaVelInterp(self.xpoints)])
        self.betaVel=betaVel
        
        if plotting==True:
            trueVel=self.jackKnifeVelocityProfile(rMin=rMin, rMax=rMax,weight=weight)
            
            plt.errorbar(trueVel[0],trueVel[1],trueVel[2])
            plt.plot(betaVel[0],betaVel[1])
            
            plt.legend(["True Velocity","Linear Theory"])
            plt.title(str(int(rMin)) + " - " + str(int(rMax)) + ", " + self.scaling + ", " + self.centerType +", beta = " + str(beta))
            plt.show()
            plt.clf()
        
        
        
        return betaVel
    "#########################################################################"
    
    "#########################################################################"
    def quadMomentAvgOverReals(self,rMin=0,rMax=100,voidTypes=["R","S"],realsToUse=range(1,101)):
        
        quadsArray=np.zeros((len(realsToUse),len(self.xpoints)))
        quadErrors=np.zeros(len(self.xpoints))
        q=0
        for real in realsToUse:
            quadsArray[q]=self.calcQuadMoment(rMin=rMin,rMax=rMax,voidTypes=voidTypes,realsToUse=[real])[1]
            q+=1
            
        quadMaster=self.calcQuadMoment(rMin=rMin,rMax=rMax,voidTypes=voidTypes,realsToUse=realsToUse)[1]
        for i in range(len(self.xpoints)):
            quadErrors[i]=np.std([quadsArray[j][i] for j in range(len(quadsArray))])
        
        
        
        self.rMinQuad=rMin
        self.rMaxQuad=rMax
        
        
        
        self.quadStd=quadErrors
        self.quadMoment=quadMaster
        
        return np.array([self.xpoints,quadMaster,quadErrors])
    "#########################################################################"
                    
    "#########################################################################"
    def jackKnifeDensityProfile(self,rMin,rMax,voidTypes=["R","S"],numSamples=100):
        
        self.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
        partsArraySamples=np.zeros((numSamples,len(self.voids[0].partsProf)))
        shellVolArraySamples=np.zeros((numSamples,len(self.voids[0].shellVols)))
        nBarSamples=np.zeros(numSamples)
        numVoidsSamples=np.zeros(numSamples)
    
        
        
        partsArrayTotal=np.zeros(len(self.voids[0].partsProf))
        shellVolTotal=np.zeros(len(self.voids[0].shellVols))
        
        
        for i in range(numSamples):
            voidsInSample=[self.voidsToUse[j] for j in range(len(self.voidsToUse)) if j%numSamples == i]
            for v in voidsInSample:
                partsArraySamples[i] += v.partsProf
                shellVolArraySamples[i] += v.shellVols
                nBarSamples[i] += v.nBar
                numVoidsSamples[i] += 1
                    
        
        
        nBarTotal=np.sum(nBarSamples)/np.sum(numVoidsSamples)
        
        for i in range(numSamples):
            partsArrayTotal += partsArraySamples[i]
            shellVolTotal += shellVolArraySamples[i]
            

        densArray=np.array([np.zeros(len(self.xpoints)) for i in range(numSamples)])
        for i in range(numSamples):
            
            nBarTemp=(np.sum(nBarSamples) - nBarSamples[i])/(np.sum(numVoidsSamples)-numVoidsSamples[i])
            partsTemp = partsArrayTotal-partsArraySamples[i]
            shellVolsTemp = shellVolTotal-shellVolArraySamples[i]
            densArrayTemp=np.array([partsTemp[k]/(shellVolsTemp[k])/nBarTemp - 1. for k in range(len(partsTemp))])
            
            # moments[i]=np.array([np.sum([densArrayTemp[k][j]*5/2*(3*muPoints[j]**2-1)*dMu for j in range(len(partsTemp[k]))]) for k in range(len(partsTemp))])
            densArray[i]=densArrayTemp
            
            
        densArrayTotal=np.array([partsArrayTotal[k]/(shellVolTotal[k])/nBarTotal - 1. for k in range(len(partsArrayTotal))])
        # meanMoment=np.array([np.sum([densArrayTotal[k][j]*5/2*(3*muPoints[j]**2-1)*dMu for j in range(len(partsArrayTotal[k]))]) for k in range(len(partsArrayTotal))])
        # meanMoment=densArrayTotal
        # moments=densArray                
    
        varSq=np.zeros(len(densArray[0]))
        for i in range(len(densArray)):
            varSq+=(len(densArray)-1)/len(densArray)*(densArray[i]-densArrayTotal)**2
            
        var=np.sqrt(varSq)
        
        self.rMinDens=rMin
        self.rMaxDens=rMax
        
        
        
        
        self.densErrors=var
        self.densProf=densArrayTotal
        
        return np.array([self.xpoints,densArrayTotal,var])
    "#########################################################################"
    
    "#########################################################################"
    def calcIntegratedDensityProfile(self,rMin,rMax,voidTypes=["R","S"]):
        
        self.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
        partsArrayTotal=np.zeros(len(self.voids[0].partsProf))
        shellVolTotal=np.zeros(len(self.voids[0].shellVols))
        nBarTotal=0.
        
        for v in self.voidsToUse:

            shellVolTotal += v.shellVols
            nBarTotal += v.nBar
            partsArrayTotal += v.partsProf
                            
        nBarTotal=nBarTotal/len(self.voidsToUse)          
        intDensArrayTotal=np.array([np.sum(partsArrayTotal[:(k+1)])/(np.sum(shellVolTotal[:(k+1)])*nBarTotal) - 1. for k in range(len(partsArrayTotal))])        
        self.Delta=np.array([self.xpoints,intDensArrayTotal])
        
        self.rMinDelta=rMin
        self.rMaxDelta=rMax
        
        
        return np.array([self.xpoints,intDensArrayTotal])
    "#########################################################################"
    
    "#########################################################################"
    def jackKnifeVelocityProfile(self,rMin,rMax,voidTypes=["R","S"],numSamples=10,weight=0):
        
        self.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
        partsArraySamples=np.zeros((numSamples,len(self.voids[0].partsProf)))
        shellVolArraySamples=np.zeros((numSamples,len(self.voids[0].shellVols)))
        nBarSamples=np.zeros(numSamples)
        numVoidsSamples=np.zeros(numSamples)
    
        partsArrayTotal=np.zeros(len(self.voids[0].partsProf))
        shellVolTotal=np.zeros(len(self.voids[0].shellVols))
        
        weightedVelsSamples=np.zeros((numSamples,len(self.xpoints)))
        weightedVelsTotal=np.zeros(len(self.xpoints))
        
        
        for i in range(numSamples):
            voidsInSample=[self.voidsToUse[j] for j in range(len(self.voidsToUse)) if j%numSamples == i]
            for v in voidsInSample:
                partsArraySamples[i] += v.partsProf
                nBarSamples[i] += v.nBar
                numVoidsSamples[i] += 1
                weightedVelsSamples[i] += np.array([v.velProf[i]*v.partsProf[i] for i in range(len(v.velProf))])*((v.radius)**weight)
                


        for i in range(numSamples):
            partsArrayTotal += partsArraySamples[i]
            shellVolTotal += shellVolArraySamples[i]
            weightedVelsTotal += weightedVelsSamples[i]
            

 
        weightedVels=np.array([np.zeros(len(self.xpoints)) for i in range(numSamples)])
        
        for i in range(numSamples):
            
            nBarTemp=(np.sum(nBarSamples) - nBarSamples[i])/(np.sum(numVoidsSamples)-numVoidsSamples[i])
            partsTemp = partsArrayTotal-partsArraySamples[i]
  
        
            weightedVelsTemp=np.zeros(len(self.xpoints))
            for j in range(len(self.xpoints)):
                if partsTemp[j] != 0:
                    weightedVelsTemp[j] = (weightedVelsTotal[j] - weightedVelsSamples[i][j])/partsTemp[j]
                    
            # moments[i]=np.array([np.sum([densArrayTemp[k][j]*5/2*(3*muPoints[j]**2-1)*dMu for j in range(len(partsTemp[k]))]) for k in range(len(partsTemp))])
            weightedVels[i]=weightedVelsTemp
            
            
        # densArrayTotal=np.array([partsArrayTotal[k]/(shellVolTotal[k])/nBarTotal - 1. for k in range(len(partsArrayTotal))])
        for i in range(len(partsArrayTotal)):
            if partsArrayTotal[i] != 0:
                weightedVelsTotal[i]=weightedVelsTotal[i]/partsArrayTotal[i] 
        
        # meanMoment=np.array([np.sum([densArrayTotal[k][j]*5/2*(3*muPoints[j]**2-1)*dMu for j in range(len(partsArrayTotal[k]))]) for k in range(len(partsArrayTotal))])
        # meanMoment=weightedVelsTotal
        # moments=weightedVels                
    
        varSq=np.zeros(len(weightedVels[0]))
        for i in range(len(weightedVels)):
            varSq+=(len(weightedVels)-1)/len(weightedVels)*(weightedVels[i]-weightedVelsTotal)**2
            
        var=np.sqrt(varSq)
        
        
        self.rMinVels=rMin
        self.rMaxVels=rMax
        self.velWeight=weight
        
        
        self.velErrors=var
        self.velProf=weightedVelsTotal
        
        return np.array([self.xpoints,weightedVelsTotal,var])
    "#########################################################################"

    "#########################################################################"
    def calcSigmaV_r(self,rMin,rMax,voidTypes=["R","S"],weight=0,startReal=1,endReal=100):
        self.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes,realsToUse=range(startReal,endReal+1))
                
        #calc mean velocity in each bin
        velProf=self.jackKnifeVelocityProfile(rMin=rMin, rMax=rMax,voidTypes=voidTypes,numSamples=10,weight=weight)[1]
        #above automatically updates voidsToUse
        numerator=np.zeros(len(self.xpoints))
        totalParts=np.zeros(len(self.xpoints))
        sigmaV_r=np.zeros(len(self.xpoints))
         
        for v in self.voidsToUse:
            parts=v.partsProf
            sigmas=v.sigmaProf*(v.radius**weight)
            voidAvgVel=v.velProf*(v.radius**weight)
             
            numerator += np.array([parts[i]*(sigmas[i]**2 + (velProf[i]-voidAvgVel[i])**2) for i in range(len(sigmas))])
            totalParts += parts
         
        for i in range(len(self.xpoints)):
            if totalParts[i] != 0:
                sigmaV_r[i] = numerator[i] / totalParts[i] 
                
        
        sigmaV_r =np.sqrt(sigmaV_r)
        self.sigmaV_r=sigmaV_r
        
        
        self.sigmaWeight=weight
        self.rMinSigma=rMin
        self.rMaxSigma=rMax
        
        
        
        
        return sigmaV_r
    "#########################################################################"
    
    "#########################################################################"
    def jackKnifeVoidSizeFunction(self,rMin=0,rMax=80,nBins=80,voidTypes=["R","S"],numSamples=10,density=True):
        self.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
        
        
        histSamplesArray=np.zeros((numSamples,nBins))
        rads=[v.radius for v in self.voidsToUse]
        hist=np.histogram(rads,bins=nBins,range=(rMin,rMax))
        
        masterHist=hist[0]
        binEdges=hist[1]
        binWidth=binEdges[1]-binEdges[0]
        binCenters=(binEdges[1:]+binEdges[:-1])/2
        
        
        for i in range(numSamples):
            voidsInSample=[self.voidsToUse[j] for j in range(len(self.voidsToUse)) if j%numSamples == i]
            radsInSample=[v.radius for v in voidsInSample]
            histSamplesArray[i]=masterHist-np.histogram(radsInSample,bins=nBins,range=(rMin,rMax))[0]
            
            "Now we need to convert HistSamplesArray into a density "
            normFac=np.sum(histSamplesArray[i])*binWidth
            histSamplesArray[i]=histSamplesArray[i]/normFac
            
        
    
    
        
        masterNormFac=np.sum(masterHist)*binWidth
        masterHistNormed=masterHist/masterNormFac
    
        varSq=np.zeros(len(histSamplesArray[0]))
        for i in range(len(histSamplesArray)):
            varSq+=(len(histSamplesArray)-1)/len(histSamplesArray)*(histSamplesArray[i]-masterHistNormed)**2
            
        var=np.sqrt(varSq)
        
        self.VSF=np.array([binCenters,masterHistNormed,var])
        
        return np.array([binCenters,masterHistNormed,var])
    "#########################################################################"
    
    "#########################################################################"
    def calcVoidSizeFunction(self,rMin=0,rMax=80,nBins=80,voidTypes=["R","S"],realsToUse=range(1,101),density=True):
        
        # self.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
        numSamples=len(realsToUse)
        histSamplesArray=np.zeros((numSamples,nBins))
        numSamples=len(realsToUse)

        
        i=0
        for real in realsToUse:
            self.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes,realsToUse=[real])
            # histSamplesArray=np.zeros((numSamples,nBins))
            rads=[v.radius for v in self.voidsToUse]
            print(len(rads))
            # hist=np.histogram(rads,bins=nBins,range=(rMin,rMax))
            histSamplesArray[i],binEdges=np.histogram(rads,bins=nBins,range=(rMin,rMax))
            # binEdges=range(1,80)
            binWidth=binEdges[1]-binEdges[0]
            histSamplesArray[i]=histSamplesArray[i]/1024**3/binWidth
            
            
            
            
            i+=1
        binCenters=(binEdges[1:] + binEdges[:-1])/2
        # return histSamplesArray,binEdges
    
        
        meanHist=[np.mean([histSamplesArray[i][j] for i in range(len(realsToUse))]) for j in range(len(histSamplesArray[0]))]
        std=[np.std([histSamplesArray[i][j] for i in range(len(realsToUse))]) for j in range(len(histSamplesArray[0]))]
        
        return binCenters,meanHist,std,histSamplesArray
    "#########################################################################"   
    
    "#########################################################################"
    def jackKnifeVSF_VolNorm(self,rMin=0,rMax=80,nBins=80,voidTypes=["R","S"],realsToUse=range(1,101),density=False):
        self.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
        
        numSamples=len(realsToUse)
        
        histSamplesArray=np.zeros((numSamples,nBins))
        rads=[v.radius for v in self.voidsToUse]
        hist=np.histogram(rads,bins=nBins,range=(rMin,rMax))
        
        masterHist=hist[0]
        binEdges=hist[1]
        binWidth=binEdges[1]-binEdges[0]
        binCenters=(binEdges[1:]+binEdges[:-1])/2
        
        masterVol=1024**3*len(realsToUse)
        
        
        for i in realsToUse:
            voidsInSample=[self.voidsToUse[j] for j in range(len(self.voidsToUse)) if self.voidsToUse[j].realNum==i] 
            
            
            radsInSample=[v.radius for v in voidsInSample]
            histSamplesArray[i-1]=masterHist-np.histogram(radsInSample,bins=nBins,range=(rMin,rMax))[0]
            
            
            
            histSamplesArray[i-1]=histSamplesArray[i-1]/(1024**3*(len(realsToUse)-1))

    
    
        
        # masterNormFac=np.sum(masterHist)*binWidth
        masterHistNormed=masterHist/masterVol
    
        
        cumHistArray=np.zeros((len(realsToUse),nBins))
        cumHistDeriv=np.zeros((len(realsToUse),nBins))
        
        for i in range(len(cumHistArray)):
            cumHistArray[i]=np.array([np.sum(histSamplesArray[i][:(j+1)]) for j in range(len(histSamplesArray[i]))])
            for j in range(1,nBins):
                cumHistDeriv[i][j]=(cumHistArray[i][j]-cumHistArray[i][j-1])/binWidth
            cumHistDeriv[i][0]=cumHistArray[i][0]/binWidth
            
                            
            
    
        masterHistNormed=masterHistNormed/binWidth
    
    
        histSamplesArray=cumHistDeriv
    
    
    
    
    
        varSq=np.zeros(len(histSamplesArray[0]))
        for i in range(len(histSamplesArray)):
            varSq+=(len(histSamplesArray)-1)/len(histSamplesArray)*(histSamplesArray[i]-masterHistNormed)**2
            
        var=np.sqrt(varSq)
        
        self.VSF=np.array([binCenters,masterHistNormed,var])
        
        return np.array([binCenters,masterHistNormed,var])
    "#########################################################################"   
    
    "#########################################################################"
    def calcBetaAndSigmaHessian(self,startIndex,endIndex,betaAndSigma=[0,0],rMin=35,rMax=100,voidTypes=["R","S"],load=False,newCalcDir="/Users/christopherwilson/Desktop/Research/GLAM/125/newCalculations/",XYZ = False,shiftFac=1,includeMono=False):
    

        "Canonical choices thus far"
        if self.scaling=="Reff":
            # endIndex=28
            scaleByReff=True
        if self.scaling=="R24":
            # endIndex=45
            scaleByReff=False
    
        catType=self.gravType + "_" + self.centerType + "_" + self.scaling
        voidTypeText="allTypes"
        

                

        print("betaAndSigma = " + str(betaAndSigma))
 


        if load==True:
            if includeMono==False:
                if XYZ==False:
                    catType=self.gravType + "_" + self.centerType + "_" + self.scaling
                    
                    Hessian=np.load(newCalcDir + "Hessian" + str(int(shiftFac)) + "_" + catType + "_" + voidTypeText + "_" + str(rMin) + "_" + str(rMax) + "_" + str(startIndex) + "to" + str(endIndex) + ".npy")
                    
                    
                    self.Hessian=Hessian
        
                    self.startIndexHessian=startIndex
                    self.endIndexHessian=endIndex
        
                    
                    return Hessian
                if XYZ == True:
                    catType=self.gravType + "_" + self.centerType + "_" + self.scaling
                    
                    Hessian=np.load(newCalcDir + "Hessian" + str(int(shiftFac)) + "XYZ_" + catType + "_" + voidTypeText + "_" + str(rMin) + "_" + str(rMax) + "_" + str(startIndex) + "to" + str(endIndex) + ".npy")
                    
                    
                    self.Hessian=Hessian
        
                    self.startIndexHessian=startIndex
                    self.endIndexHessian=endIndex
        
                    
                    return Hessian
                
            if includeMono==True:
                            
                if XYZ==False:
                    catType=self.gravType + "_" + self.centerType + "_" + self.scaling
                    
                    Hessian=np.load(newCalcDir + "Hessian" + str(int(shiftFac)) + "_" + "inclMono_" + catType  + "_" + voidTypeText + "_" + str(rMin) + "_" + str(rMax) + "_" + str(startIndex) + "to" + str(endIndex) + ".npy")
                    
                    
                    self.Hessian=Hessian
        
                    self.startIndexHessian=startIndex
                    self.endIndexHessian=endIndex
        
                    
                    return Hessian
                if XYZ == True:
                    catType=self.gravType + "_" + self.centerType + "_" + self.scaling
                    
                    Hessian=np.load(newCalcDir + "Hessian" + str(int(shiftFac)) + "XYZ_" + "inclMono_" + catType  + "_" + voidTypeText + "_" + str(rMin) + "_" + str(rMax) + "_" + str(startIndex) + "to" + str(endIndex) + ".npy")
                    
                    
                    self.Hessian=Hessian
        
                    self.startIndexHessian=startIndex
                    self.endIndexHessian=endIndex
        
                    
                    return Hessian


        "Construct Quad Data, if not already constructed, with JackKnife errors" 
        
        Q=0
        try:
            if (self.rMinQuad==rMin and self.rMaxQuad==rMax):
                print("Loading Quadrupole")
                quadMoment=np.array([self.xpoints,self.quadMoment,self.quadErrors])
                Q=1
            if Q==0:
                print("quadrupole doesnt match required params, calculating now")
                # velProf=self.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
                quadMoment=self.jackKnifeQuadMoment(rMin=rMin,rMax=rMax,voidTypes=voidTypes,numSamples=1000)
                
        except:
            print("quadrupole not pre-calculated, calculating now")
            quadMoment=self.jackKnifeQuadMoment(rMin=rMin,rMax=rMax,voidTypes=voidTypes,numSamples=1000)



        if includeMono==True:
            Q=0
            try:
                if (self.rMinMono==rMin and self.rMaxMono==rMax):
                    print("Loading Quadrupole")
                    monoMoment=np.array([self.xpoints,self.monoMoment,self.monoErrors])
                    Q=1
                if Q==0:
                    print("quadrupole doesnt match required params, calculating now")
                    # velProf=self.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
                    monoMoment=self.jackKnifeMonopoleMoment(rMin=rMin,rMax=rMax,voidTypes=voidTypes,numSamples=1000)
                    
            except:
                print("quadrupole not pre-calculated, calculating now")
                monoMoment=self.jackKnifeMonopoleMoment(rMin=rMin,rMax=rMax,voidTypes=voidTypes,numSamples=1000)
            


        "##################" 
        def betaSigmaFunc(betaAndSigma,quadData=quadMoment,rMin=rMin,rMax=rMax,a=1/(1+0.48251),H=130.29171508616594,voidType=["R","S"],scaleByReff=scaleByReff,startIndex=startIndex,endIndex=endIndex,nMuBins=50,includeMono=includeMono):
            beta=betaAndSigma[0]
            sigma=betaAndSigma[1]
            sigmas=sigma*np.ones(len(self.xpoints))
            lambdas=sigma*np.zeros(len(self.xpoints))
            kappas=sigma*np.zeros(len(self.xpoints))
            
    
            quadThy,deltaSandMu=calcQuadEdgeworthStreaming(cat=self,sigmas=sigmas,lambdas=lambdas,kappas=kappas,voidTypes=voidTypes,catType="halo",rMin=rMin,rMax=rMax,modelVelocity=True,beta=beta,a=a,H=H,Integrated=True,nMuBins=nMuBins,muDependance=False,degree=0,velFac=1,startIndex=startIndex,endIndex=endIndex,shiftFac=shiftFac)
            plt.plot(quadThy[0],quadThy[1])
            plt.errorbar(quadThy[0],quadData[1],quadData[2])
            plt.title("beta = " + str(beta) + ", sigma = " + str(sigma))
            plt.show()
            plt.clf()
    
        
            chiSquared=calcChiSquaredForModel(modelData=quadThy[1],realData=quadData[1],dataErrors=quadData[2],startIndex=startIndex,endIndex=endIndex)
            
            if includeMono==True:
                monoThy=np.mean(deltaSandMu,axis=1)-1
                chiSquared2=calcChiSquaredForModel(modelData=monoThy,realData=monoMoment[1],dataErrors=monoMoment[2],startIndex=startIndex,endIndex=endIndex)
                chiSquared = chiSquared + chiSquared2
                
                
                
        
            print("chiSquared = " + str(chiSquared))
            print("beta is " + str(betaAndSigma[0]))
            print("sigma is " + str(betaAndSigma[1]))
            
            
            
            return chiSquared
        "##################" 
        
        

        
        
        "##################" 
        def calcChiSquaredForModel(modelData,realData,dataErrors,startIndex=startIndex,endIndex=endIndex):
            "4, 30 or 5, 45"
            if (len(modelData) != len(realData)) and (len(dataErrors) != len(realData)):
                print("input length mismatch")
                return -1
            
            print("start index is " + str(startIndex))
            print("end index is " + str(endIndex))
            chiSquared=0
            for i in range(startIndex,endIndex):
            # for i in range(4,30):
                chiSquared += (modelData[i] - realData[i])**2/(dataErrors[i]**2)
            return chiSquared
        "##################" 
    
    
    

  
        beta=betaAndSigma[0]
        sigma=betaAndSigma[1]
        
        dBeta=0.02
        dSigma=0.1
        
        if self.scaling=="R24":
            dSigma=10.
        
        
        
        # F_beta_sigma=optimizeBetaAndSigma([beta,sigma],cat=cat,quadData=quadData,rMin=rMin,rMax=rMax,a=a,H=H,voidType=voidType,scaleByReff=scaleByReff,startIndex=startIndex,endIndex=endIndex,nMuBins=nMuBins)
        F_beta_sigma=betaSigmaFunc([beta,sigma])
        print("d.o.f = " + str(endIndex-startIndex))
        print("USE THIS CHI SQUARED")
        "###"
        F_betaPlus_sigma=betaSigmaFunc([beta + dBeta,sigma])#,cat=cat,quadData=quadData,rMin=rMin,rMax=rMax,a=a,H=H,voidType=voidType,scaleByReff=scaleByReff,startIndex=startIndex,endIndex=endIndex,nMuBins=nMuBins)
        F_betaMinus_sigma=betaSigmaFunc([beta - dBeta,sigma])#,cat=cat,quadData=quadData,rMin=rMin,rMax=rMax,a=a,H=H,voidType=voidType,scaleByReff=scaleByReff,startIndex=startIndex,endIndex=endIndex,nMuBins=nMuBins)
        F_beta_sigmaPlus=betaSigmaFunc([beta,sigma+dSigma])#,cat=cat,quadData=quadData,rMin=rMin,rMax=rMax,a=a,H=H,voidType=voidType,scaleByReff=scaleByReff,startIndex=startIndex,endIndex=endIndex,nMuBins=nMuBins)
        F_beta_sigmaMinus=betaSigmaFunc([beta,sigma-dSigma])#,cat=cat,quadData=quadData,rMin=rMin,rMax=rMax,a=a,H=H,voidType=voidType,scaleByReff=scaleByReff,startIndex=startIndex,endIndex=endIndex,nMuBins=nMuBins)
        "###"
        
        "###"
        F_betaPlus_sigmaPlus=betaSigmaFunc([beta + dBeta,sigma + dSigma])#,cat=cat,quadData=quadData,rMin=rMin,rMax=rMax,a=a,H=H,voidType=voidType,scaleByReff=scaleByReff,startIndex=startIndex,endIndex=endIndex,nMuBins=nMuBins)
        F_betaPlus_sigmaMinus=betaSigmaFunc([beta + dBeta,sigma - dSigma])#,cat=cat,quadData=quadData,rMin=rMin,rMax=rMax,a=a,H=H,voidType=voidType,scaleByReff=scaleByReff,startIndex=startIndex,endIndex=endIndex,nMuBins=nMuBins)
        F_betaMinus_sigmaPlus=betaSigmaFunc([beta - dBeta,sigma + dSigma])#,cat=cat,quadData=quadData,rMin=rMin,rMax=rMax,a=a,H=H,voidType=voidType,scaleByReff=scaleByReff,startIndex=startIndex,endIndex=endIndex,nMuBins=nMuBins)
        F_betaMinus_sigmaMinus=betaSigmaFunc([beta - dBeta,sigma - dSigma])#,cat=cat,quadData=quadData,rMin=rMin,rMax=rMax,a=a,H=H,voidType=voidType,scaleByReff=scaleByReff,startIndex=startIndex,endIndex=endIndex,nMuBins=nMuBins)
        "###"
    
    
    
        # print(dBeta)
        # print(dSigma)
        # print(F_beta_Sigma)
        # print(F_betaPlus_sigma )
    
        # dBeta=0.015
        # dSigma=0.05
    
    
        d2beta_F=1/(dBeta**2)*(-2*F_beta_sigma + F_betaPlus_sigma + F_betaMinus_sigma)
        d2sigma_F=1/(dSigma**2)*(-2*F_beta_sigma + F_beta_sigmaPlus + F_beta_sigmaMinus)
        dbeta_dsigma_F=1/(4*dBeta*dSigma)*(F_betaPlus_sigmaPlus + F_betaMinus_sigmaMinus - F_betaPlus_sigmaMinus -  F_betaMinus_sigmaPlus)
        
        Hessian=np.array([[d2beta_F, dbeta_dsigma_F],[dbeta_dsigma_F,d2sigma_F]])
        
        
        catType=self.gravType + "_" + self.centerType + "_" + self.scaling
        
        self.Hessian=Hessian

        self.startIndexHessian=startIndex
        self.endIndexHessian=endIndex

        
        return Hessian
    "#########################################################################"
    
    "#########################################################################"        
    def calcConfidenceEllipses(self,startIndex,endIndex,rMin=35,rMax=100,voidTypes=["R","S"],calcDir="/Users/christopherwilson/Desktop/Research/GLAM/125/newCalculations/",load=False,loadHessian=False,save=False,errors="new",XYZ=True,shiftFac=1,tol=0.05,includeMono=False):
        
        

        
        #canonical choices for endIndex thus far 
        if self.scaling=="Reff":
            # endIndex=28
            scaleByReff=True
        if self.scaling=="R24":
            # endIndex=45
            scaleByReff=False
        
        
        
        catType=self.gravType + "_" + self.centerType + "_" + self.scaling
        voidTypeText="allTypes"    
        
        
        if includeMono==False:
            if load==False:
                optParams=self.optimizeConstantBetaAndSigma(startIndex=startIndex,endIndex=endIndex,rMin=rMin,rMax=rMax,voidTypes=voidTypes,shiftFac=shiftFac,tol=tol,includeMono=False)
            if load==True:
                if XYZ==False:
                    optParams=np.load(calcDir + "betaAndSigma" + str(int(shiftFac)) + "_" + catType + "_" + voidTypeText + "_" + str(rMin) + "_" + str(rMax) + "_" + str(startIndex) + "to" + str(endIndex) + ".npy")
                if XYZ==True:
                    optParams=np.load(calcDir + "betaAndSigma" + str(int(shiftFac)) + "XYZ_" + catType + "_" + voidTypeText + "_" + str(rMin) + "_" + str(rMax) + "_" + str(startIndex) + "to" + str(endIndex) + ".npy")
            
        if includeMono==True:
            if load==False:
                optParams=self.optimizeConstantBetaAndSigma(startIndex=startIndex,endIndex=endIndex,rMin=rMin,rMax=rMax,voidTypes=voidTypes,shiftFac=shiftFac,tol=tol,includeMono=True)
        
            if load==True:
                if XYZ==False:
                    optParams=np.load(calcDir + "betaAndSigma" + str(int(shiftFac)) + "_" + "inclMono_" + catType + "_" + voidTypeText + "_" + str(rMin) + "_" + str(rMax) + "_" + str(startIndex) + "to" + str(endIndex) + ".npy")
                if XYZ==True:
                    optParams=np.load(calcDir + "betaAndSigma" + str(int(shiftFac)) + "XYZ_" + "inclMono_" + catType + "_" + voidTypeText + "_" + str(rMin) + "_" + str(rMax) + "_" + str(startIndex) + "to" + str(endIndex) + ".npy")
            
        
        
        
        
        print("successfully optimized for beta and sigma")
        

        
        Hessian=self.calcBetaAndSigmaHessian(startIndex=startIndex,endIndex=endIndex,betaAndSigma=optParams,rMin=rMin,rMax=rMax,voidTypes=voidTypes,load=loadHessian,XYZ = XYZ,shiftFac=shiftFac,newCalcDir=calcDir,includeMono=includeMono)
    
    
    
        if save==True:
            if XYZ==False:
                if includeMono==False:
                    np.save(newCalcDir + "Hessian" + str(int(shiftFac)) + "_" + catType + "_" + voidTypeText + "_" + str(int(rMin)) + "_" + str(int(rMax)) + "_" + str(int(startIndex)) + "to" + str(int(endIndex)),Hessian)
                    np.save(newCalcDir + "betaAndSigma" + str(int(shiftFac)) + "_" + catType + "_" + voidTypeText + "_" + str(int(rMin)) + "_" + str(int(rMax)) + "_" + str(int(startIndex)) + "to" + str(int(endIndex)),optParams)
                    
                if includeMono==True:
                    np.save(newCalcDir + "Hessian" + str(int(shiftFac)) + "_" + "inclMono_" + catType + "_" + voidTypeText + "_" + str(int(rMin)) + "_" + str(int(rMax)) + "_" + str(int(startIndex)) + "to" + str(int(endIndex)),Hessian)
                    np.save(newCalcDir + "betaAndSigma" + str(int(shiftFac)) + "_" + "inclMono_" + catType + "_" + voidTypeText + "_" + str(int(rMin)) + "_" + str(int(rMax)) + "_" + str(int(startIndex)) + "to" + str(int(endIndex)),optParams)
                    
                    
            if XYZ==True:
                if includeMono==False:
                    np.save(newCalcDir + "Hessian" + str(int(shiftFac)) + "XYZ_" + "inclMono_" + catType + "_" + voidTypeText + "_" + str(int(rMin)) + "_" + str(int(rMax)) + "_" + str(int(startIndex)) + "to" + str(int(endIndex)),Hessian)
                    np.save(newCalcDir + "betaAndSigma"  + str(int(shiftFac)) + "XYZ_" + "inclMono_" + catType + "_" + voidTypeText + "_" + str(int(rMin)) + "_" + str(int(rMax)) + "_" + str(int(startIndex)) + "to" + str(int(endIndex)),optParams)
                if includeMono==True:
                    np.save(newCalcDir + "Hessian" + str(int(shiftFac)) + "XYZ_" + "inclMono_" + catType + "_" + voidTypeText + "_" + str(int(rMin)) + "_" + str(int(rMax)) + "_" + str(int(startIndex)) + "to" + str(int(endIndex)),Hessian)
                    np.save(newCalcDir + "betaAndSigma"  + str(int(shiftFac)) + "XYZ_" + "inclMono_" + catType + "_" + voidTypeText + "_" + str(int(rMin)) + "_" + str(int(rMax)) + "_" + str(int(startIndex)) + "to" + str(int(endIndex)),optParams)
    
    
        if errors=="new":
            Hessian = Hessian / 100

        if XYZ==True:
            Hessian = Hessian / 3


        "Confidence Ellipses"
            
        alpha=0.5*Hessian
        

        
        
        # optParams=self.betaAndSigma
        

        deltaChi2=2.3 #2.3 for a 68% confidence ellipse
        deltaChi2_95=1
        
        
        lambdas=np.linalg.eig(alpha)[0]
        weights=np.sqrt(lambdas)
        vecs=np.transpose(np.linalg.eig(alpha)[1])
        
        
        def deltaA(z,deltaChi2,weights,vecs):
            "z must lie on unit circle"
            dA=np.zeros(len(vecs))
            for i in range(len(weights)):
                dA+=np.sqrt(deltaChi2)*1./weights[i]*(z@vecs[i])*vecs[i]
            return dA
        
    
        def quadraticForm(alphaMatrix,dA):
            "returns x @ alphaMatric @ x"
            return dA @ (alphaMatrix@dA)
            
            
        zList=np.array([[np.cos(theta),np.sin(theta)] for theta in np.linspace(0,2*np.pi,1000)])
        dAList=[deltaA(z,deltaChi2=deltaChi2,weights=weights,vecs=vecs) for z in zList]
        dAList_95=[deltaA(z,deltaChi2=deltaChi2_95,weights=weights,vecs=vecs) for z in zList]
                
        
        dAList_NC=np.array([dAList[i]+optParams for i in range(len(dAList))])
        dAList_95_NC=np.array([dAList_95[i]+optParams for i in range(len(dAList))])
        # plt.plot([dAList_NC[i][0] for i in range(len(dAList_NC))],[dAList_NC[i][1] for i in range(len(dAList_NC))])
        # plt.plot(optParams[0],optParams[1],marker="*",color="C0")
        
        dAList_NC=np.transpose(dAList_NC)
        dAList_95_NC=np.transpose(dAList_95_NC)
        
        
        catType=self.gravType + "_" + self.centerType + "_" + self.scaling
        voidTypeText="allTypes"    
        

            
        self.dAList=dAList_NC
        self.dAList_95=dAList_95_NC
        return dAList_NC,optParams
    
    "#########################################################################" 
    "#########################################################################" 
    
    "#########################################################################"
    def optimizeConstantBetaAndSigma(self,startIndex,endIndex,guess=None,rMin=0,rMax=100,voidTypes=["R","S"],shiftFac=1,tol=0.05,includeMono=False):
        "endIndex is fixed depending on the scaling of the void "
        
        self.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
        
        #canonical choices for endIndex thus far
        
        
        if self.scaling=="Reff":
            betaAndSigma=[0.35,7.]        
            # endIndex=28
            scaleByReff=True
            
            
        if self.scaling=="R24":
            betaAndSigma=[0.35,335.]
            # endIndex=45
            scaleByReff=False
            
        "Construct Quad Data, if not already constructed, with JackKnife errors"    
        Q=0
        try:
            if (self.rMinQuad==rMin and self.rMaxQuad==rMax):
                print("Loading Quadrupole")
                quadMoment=np.array([self.xpoints,self.quadMoment,self.quadErrors])
                Q=1
            if Q==0:
                print("quadrupole doesnt match required params, calculating now")
                # velProf=self.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
                quadMoment=self.jackKnifeQuadMoment(rMin=rMin,rMax=rMax,voidTypes=voidTypes,numSamples=1000)
                
        except:
            print("quadrupole not pre-calculated, calculating now")
            quadMoment=self.jackKnifeQuadMoment(rMin=rMin,rMax=rMax,voidTypes=voidTypes,numSamples=1000)
        

        if includeMono==True:
            Q=0
            try:
                if (self.rMinMono==rMin and self.rMaxMono==rMax):
                    print("Loading Quadrupole")
                    monoMoment=np.array([self.xpoints,self.monoMoment,self.monoErrors])
                    Q=1
                if Q==0:
                    print("quadrupole doesnt match required params, calculating now")
                    # velProf=self.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
                    monoMoment=self.jackKnifeMonopoleMoment(rMin=rMin,rMax=rMax,voidTypes=voidTypes,numSamples=1000)
                    
            except:
                print("quadrupole not pre-calculated, calculating now")
                monoMoment=self.jackKnifeMonopoleMoment(rMin=rMin,rMax=rMax,voidTypes=voidTypes,numSamples=1000)
            


        "##################" 
        def calcChiSquaredForModel(modelData,realData,dataErrors,startIndex=startIndex,endIndex=endIndex):
            "4, 30 or 5, 45"
            if (len(modelData) != len(realData)) and (len(dataErrors) != len(realData)):
                print("input length mismatch")
                return -1
            
            print("start index is " + str(startIndex))
            print("end index is " + str(endIndex))
            chiSquared=0
            for i in range(startIndex,endIndex):
            # for i in range(4,30):
                chiSquared += (modelData[i] - realData[i])**2/(dataErrors[i]**2)
            return chiSquared
        "##################" 




        "##################" 
        def betaSigmaFunc(betaAndSigma,quadData=quadMoment,rMin=rMin,rMax=rMax,a=1/(1+0.48251),H=130.29171508616594,voidType=["R","S"],scaleByReff=scaleByReff,startIndex=startIndex,endIndex=endIndex,nMuBins=50,includeMono=includeMono):
            beta=betaAndSigma[0]
            sigma=betaAndSigma[1]
            sigmas=sigma*np.ones(len(self.xpoints))
            lambdas=sigma*np.zeros(len(self.xpoints))
            kappas=sigma*np.zeros(len(self.xpoints))
            
    
            quadThy,deltaSandMu=calcQuadEdgeworthStreaming(cat=self,sigmas=sigmas,lambdas=lambdas,kappas=kappas,voidTypes=voidTypes,catType="halo",rMin=rMin,rMax=rMax,modelVelocity=True,beta=beta,a=a,H=H,Integrated=True,nMuBins=nMuBins,muDependance=False,degree=0,velFac=1,startIndex=startIndex,endIndex=endIndex,shiftFac=shiftFac)
            # plt.plot(quadThy[0],quadThy[1])
            # plt.errorbar(quadThy[0],quadData[1],quadData[2])
            # plt.title("beta = " + str(beta) + ", sigma = " + str(sigma))
            # plt.show()
            # plt.clf()
    
        
            chiSquared=calcChiSquaredForModel(modelData=quadThy[1],realData=quadData[1],dataErrors=quadData[2],startIndex=startIndex,endIndex=endIndex)
            
            if includeMono==True:
                monoThy=np.mean(deltaSandMu,axis=1)-1
                chiSquared2=calcChiSquaredForModel(modelData=monoThy,realData=monoMoment[1],dataErrors=monoMoment[2],startIndex=startIndex,endIndex=endIndex)
                chiSquared = chiSquared + chiSquared2
                # plt.plot(quadThy[0],monoThy)
                # plt.errorbar(quadThy[0],monoMoment[1],monoMoment[2])
                # plt.title("beta = " + str(beta) + ", sigma = " + str(sigma))
                # plt.show()
                # plt.clf()
                
                
                
                
        
            print("chiSquared = " + str(chiSquared))
            print("beta is " + str(betaAndSigma[0]))
            print("sigma is " + str(betaAndSigma[1]))
            
            
            
            return chiSquared
        "##################" 
            

                
        optimObject=scipy.optimize.minimize(betaSigmaFunc,x0=betaAndSigma, method=None, jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=tol, callback=None, options=None)
        betaAndSigmaOpt=optimObject.x
        
        self.optimObject=optimObject
        self.betaAndSigmaOpt=betaAndSigmaOpt
        print("Successfully terminated optimization!")

        return betaAndSigmaOpt
    "#########################################################################"
    
    "#########################################################################"
    def calcFullGSM(self,rMin=25,rMax=100,voidTypes=["R","S"],sigmaFac=1.0,velFac=1.0,numBinFac=1,sigmas=None):
        
        
        if self.scaling=="Reff":
            weight=-1
        if self.scaling=="R24":
            weight=0
        
        
        
        "if sigmas is given as 'None' above, then we calculate sigmas from the catalogs. If anything else is given as sigmas, that array will be used as sigmas instead"
        if sigmas is None:
            

            Q=0
            try:
                if self.rMinSigma==rMin and self.rMaxSigma==rMax and self.sigmaWeight==weight:
                    print("loading sigma profile")
                    sigmas=self.sigmaV_r*sigmaFac
                    
                    Q=1
                if Q==0:
                    print("Sigma profile doesnt match required params, calculating now")
                    sigmas = self.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)*sigmaFac
                    

                    
            except:
                print("sigma Prof not pre-calculated, calculating now")
                sigmas = self.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)*sigmaFac    
    
    


        "0's in sigmas is a problem, have to set them equal to something small but nonzero"
        
        if weight==0:
            for i in range(len(sigmas)):
                if sigmas[i]==0:
                    sigmas[i]=1.
        if weight==-1:
            for i in range(len(sigmas)):
                if sigmas[i]==0:
                    sigmas[i]=0.2


        quadThy,deltaSandMuThy=calcQuadEdgeworthStreaming(self,
                               sigmas=sigmas,
                               lambdas=np.zeros(len(self.xpoints)),
                               kappas=np.zeros(len(self.xpoints)),
                               voidTypes=voidTypes,
                               catType="halo",
                               rMin=rMin,
                               rMax=rMax,
                               modelVelocity=False,
                               calcVelocityForPlot=False,
                               beta=1,
                               a=1/(1+0.48251),
                               H=130.29171508616594,
                               Integrated=True,
                               nMuBins=50,
                               numBinFac=numBinFac,
                               muDependance=False,
                               degree=0,
                               velFac=velFac,
                               startIndex=0,
                               endIndex=round(len(self.xpoints)*numBinFac))
        
        self.quadThy=quadThy
        self.deltaSandMuThy=deltaSandMuThy
        
        return quadThy
    "#########################################################################"
    
    "#########################################################################"
    def calcNoIntGSM(self,rMin=25,rMax=100,voidTypes=["R","S"],velFac=1.0,numBinFac=1):
        
        
        if self.scaling=="Reff":
            weight=-1
        if self.scaling=="R24":
            weight=0

        sigmas=np.zeros(30)
    

        quadThy,deltaSandMuThy=calcQuadEdgeworthStreaming(self,
                               sigmas=sigmas,
                               lambdas=np.zeros(len(self.xpoints)),
                               kappas=np.zeros(len(self.xpoints)),
                               voidTypes=voidTypes,
                               catType="halo",
                               rMin=rMin,
                               rMax=rMax,
                               modelVelocity=False,
                               calcVelocityForPlot=False,
                               beta=1,
                               a=1/(1+0.48251),
                               H=130.29171508616594,
                               Integrated=False,
                               nMuBins=50,
                               numBinFac=numBinFac,
                               muDependance=False,
                               degree=0,
                               velFac=velFac,
                               startIndex=0,
                               endIndex=round(len(self.xpoints)*numBinFac))
        
        self.quadThy=quadThy
        self.deltaSandMuThy=deltaSandMuThy
        
        return quadThy
    "#########################################################################"
    
    "#########################################################################"
    def calcConstGSM(self,betaAndSigma,rMin=25,rMax=100,voidTypes=["R","S"],numBinFac=1):
        
    

        beta=betaAndSigma[0]
        sigmas=np.ones(len(self.xpoints))*betaAndSigma[1]
        


        quadThy,deltaSandMuThy=calcQuadEdgeworthStreaming(self,
                               sigmas=sigmas,
                               lambdas=np.zeros(len(self.xpoints)),
                               kappas=np.zeros(len(self.xpoints)),
                               voidTypes=voidTypes,
                               catType="halo",
                               rMin=rMin,
                               rMax=rMax,
                               modelVelocity=True,
                               calcVelocityForPlot=False,
                               beta=beta,
                               a=1/(1+0.48251),
                               H=130.29171508616594,
                               Integrated=True,
                               nMuBins=100,
                               muDependance=False,
                               degree=0,
                               velFac=1,
                               startIndex=0,
                               numBinFac=numBinFac,
                               endIndex=len(self.xpoints))
        
        self.quadThyConst=quadThy
        self.deltaSandMuThyConst=deltaSandMuThy
        
        return quadThy
    "#########################################################################"
    
    "#########################################################################"
    def optimizeBeta(self,rMin=25,rMax=100,voidTypes=["R","S"],startIndex=4,plotting=False,shiftFac=1):
        
        
        xpoints=self.xpoints
        
        dx=(xpoints[1]-xpoints[0])/2
        Delta=self.calcIntegratedDensityProfile(rMin=rMin,rMax=rMax,voidTypes=["R","S"])
        
        if self.scaling=="Reff":
            weight=-1
        if self.scaling=="R24":
            weight=0
        
        
        trueVel=self.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,numSamples=1000,weight=weight)
        def func(beta,startIndex=startIndex,endIndex=len(self.xpoints)-5):
            
            # startIndex=np.argmax(trueVel[1])-1
            # betaVelInterp=sp.interpolate.CubicSpline((Delta[0])+dx*shiftFac,-1/3*a*H*beta*Delta[1]*(Delta[0]+dx),bc_type = 'natural')
            
            betaVelInterp=interpFunc(Delta,beta,shiftFac=shiftFac)
            chi2=np.sum([((betaVelInterp(xpoints[i]) - trueVel[1][i])/trueVel[2][i])**2 for i in range(startIndex,endIndex)])
            
            return chi2
        
        
        optimObject=scipy.optimize.minimize(func,x0=[0.42], method=None, jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=0.0001, callback=None, options=None)
        
        beta=optimObject.x
        
        
        
        if plotting==True:
            
            betaVelInterp=interpFunc(Delta,beta,shiftFac=shiftFac)
            
            plt.errorbar(trueVel[0],trueVel[1],trueVel[2])
            plt.plot(xpoints,betaVelInterp(xpoints))
            
            plt.legend(["Linear Thy","V_r"])
            plt.title(self.gravType + ", " + str(rMin) + " - " + str(rMax) + ", beta = " + str(int(beta*1000)/1000))
        
            plt.show()
            plt.clf()
            
        return beta,optimObject
    "#########################################################################"
    

"#########################################################################"       

"#########################################################################"
class voidClass:
    
    "#########################################################################"
    def __init__(self,voidInfo,partsRandMu=None,densVelsSigmas=[None,None,None],velZArray=None,scaling="",xpoints=[]):
        
        if isinstance(voidInfo, (list, tuple, np.ndarray)):

            # print(voidInfo)
            self.radius=voidInfo[2]
            self.center=voidInfo[0][3]
            self.realNum=voidInfo[0][0]
            self.voidID=voidInfo[0][1]
            self.nBar=voidInfo[0][2]
            self.voidType=voidInfo[1]
            self.partsRandMu=partsRandMu
            
            self.densProf=densVelsSigmas[0]
            self.velProf=densVelsSigmas[1]
            self.sigmaProf=densVelsSigmas[2]
            self.velZArray=velZArray
            
            # print(self.velProf)
            
            if scaling=="Reff":
                dx=(xpoints[1]-xpoints[0])/2
                shellVols=4./3*np.pi*(self.radius**3)*np.array([(xpoints[i]+dx)**3-(xpoints[i]-dx)**3 for i in range(len(xpoints))]) 
                self.shellVols=shellVols
            
            if scaling=="R24":
                dx=(xpoints[1]-xpoints[0])/2
                shellVols=4./3*np.pi*np.array([(xpoints[i]+dx)**3-(xpoints[i]-dx)**3 for i in range(len(xpoints))]) 
                self.shellVols=shellVols
                
            if isinstance(self.densProf, (list, tuple, np.ndarray)):
                size=len(shellVols)
                self.partsProf=np.array([round((self.densProf[i]+1)*self.shellVols[i]*self.nBar) for i in range(size)])
    "#########################################################################"   
"#########################################################################"            
    

"#########################################################################"     


"Methods which are best to define outside of the catalog class environment"

#######
"#########################################################################"
def coordinateChangeRtoS(cat,r,muR,rMin=25,rMax=100,a=a,H=H,voidTypes=["R","S"]):
    
    
    Q=0
    try:
        if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    
    
    velRInterp=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')

    
    def deltaInterp(x):
        if (x>0 and x<xpoints[-1]):
            return deltaRInterp(x)
        if (x<0):
            return -1.
        return 0.
    
    
    def velInterp(x):
    
        if (x>0 or x<xpoints[-1]):
            return velRInterp(x)
        return 0.
    
    
    
    
    
    
    sPerp=r*np.sqrt((1-muR**2))
    
    sPara=r*muR + muR*velInterp(r)/(a*H)
    
    s=np.sqrt(sPerp**2 + sPara**2)

    muS=sPara/s
    
    print("s = " + str(s))
    print("muS = " + str(muS))
    
    return [s,muS]
"#########################################################################"

"#########################################################################"
def coordinateChangeStoR(cat,s,muS,rMin=25,rMax=100,a=a,H=H,voidTypes=["R","S"]):
    
    
    
    if cat.scaling=="R24":
        weight=0

    if cat.scaling=="Reff":
        weight=-1
    
    xpoints=cat.xpoints

    Q=0
    try:
        if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    
    
    velRInterp=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')

    

    
    def velInterp(x):
    
        if (x>0 or x<xpoints[-1]):
            return velRInterp(x)
        return 0.
    
    
    
    ###


    sPara=s*muS
    sPerp=np.sqrt(s**2-sPara**2)
    
    rPerp=sPerp

        
    def rParaFunc(rPara):
        
        r=np.sqrt(rPerp**2+rPara**2)
        muR=rPara/r
        return rPara+(velInterp(r)*muR)/(a*H)-sPara
    
        
    rPara=fsolve(rParaFunc, x0=[sPara])[0]
    
    r=np.sqrt(rPerp**2+rPara**2)
    muR=rPara/r

    print("r = " + str(r))
    print("muR = " + str(muR))


    return [r,muR]
"#########################################################################"

"#########################################################################"
def numericalDiff(func,x,*args,dx=0.1,order=1):
    "we numerically differentiate x, which is always the first argument of func."

    "first we have to construct the points at which we wish to differentiate"

    "if order is even..."
    if order%2==0:
        numPoints=(order+1)*2+1

    "if order is odd..."
    if order%2==1:
        numPoints=(order+2)*2+1

    # maxToSub=round((numPoints-1)/2)
    # "points is working correctly"
    # points=[x-maxToSub*dx + i*dx for i in range(numPoints)]

    argsToUse=[AA for AA in args]
    
    return diff(func,x0=x,dx=dx,n=order,args=argsToUse,order=numPoints)
"#########################################################################"

"#########################################################################"
def nD(func,x,*args,dx=0.1,order=1):
    "we numerically differentiate x, which is always the first argument of func."

    "first we have to construct the points at which we wish to differentiate"

    "if order is even..."
    if order%2==0:
        numPoints=(order+1)*2+1

    "if order is odd..."
    if order%2==1:
        numPoints=(order+2)*2+1

    # maxToSub=round((numPoints-1)/2)
    # "points is working correctly"
    # points=[x-maxToSub*dx + i*dx for i in range(numPoints)]

    argsToUse=[AA for AA in args]
    
    return diff(func,x0=x,dx=dx,n=order,args=argsToUse,order=numPoints)
"#########################################################################"

"#########################################################################"
def optimizeBetaFixedShift(cat,rMin=35,rMax=100,voidTypes=["R","S"],plotting=True,fitInterior=False,fitExterior=False,fitWholeVoid=True,intFac=1000,shiftFac=1,fitIndex=10):
    
    "shiftFac enters as ...(Delta[0])+(dx*(1+shiftFac))..."
    
    
    shiftFac=shiftFac
    
    xpoints=cat.xpoints
    
    dx=(xpoints[1]-xpoints[0])/2
    # Delta=cat.calcIntegratedDensityProfile(rMin=rMin,rMax=rMax,voidTypes=["R","S"])
    Q=0
    try:
        if cat.rMinDelta==rMin and cat.rMaxDelta==rMax:
            Delta=cat.Delta
            print("loaded Delta")
            Q=1
    except:
        pass
    
    if Q==0:
        print("calculating Delta")
        Delta=cat.calcIntegratedDensityProfile(rMin=rMin,rMax=rMax,voidTypes=["R","S"])
    
    
    if cat.scaling=="Reff":
        weight=-1
    if cat.scaling=="R24":
        weight=0
    
    
    
    
    Q=0
    try:
        if cat.rMinVels==rMin and cat.rMaxVels==rMax:
            trueVel=np.array([cat.xpoints,cat.velProf,cat.velErrors])
            print("loaded Velocity")
            Q=1
    except:
        pass
    
    
    if Q==0:
        print("calculating Velocity")
        trueVel=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,numSamples=1000,weight=weight)
    
    

    
    
    endIndexTemp=len(trueVel[1])-1
    
    if fitInterior==True:
        fitWholeVoid=False


        for i in range(len(trueVel[1])):
            if trueVel[1][i] != 0:
                startIndex=i
                break
    
        
        maxVel=np.max(trueVel)
        endIndex=fitIndex
        if False:
            endIndex=endIndexTemp
            for i in range(startIndex,len(trueVel[1])):
                if (trueVel[1][i] < maxVel/intFac):
                    endIndex=i
                    break
            
        fitTypeText="interior fit only"
        
    
    
    if fitExterior==True:
        fitWholeVoid=False
        
        startIndex=fitIndex
    
        if False:
            maxVel=np.max(trueVel)
            for i in range(startIndex,len(trueVel[1])):
                if (trueVel[1][i] < maxVel/intFac):
                    startIndex=i
                    break    
                
        endIndex=len(trueVel[1])-1
        fitTypeText="exterior fit only"
    

    "################################################"
    if fitWholeVoid==True:
        if False:
            startIndex=np.argmax(trueVel[1])-1
        if True:
            for i in range(len(trueVel[1])):
                if trueVel[1][i] != 0:
                    startIndex=i
                    break
        
        #"return here"#
        endIndex=len(trueVel[1])-1
        fitTypeText="whole void fit"
    "################################################"
    
    
    
    # endIndexTemp=len(trueVel[1])-1
    if cat.scaling=="Reff" and fitExterior==True:
        startIndex=np.min([15,startIndex])
        if rMin >= 40:
            startIndex=15
    
    
    if cat.scaling=="R24" and fitExterior==True:
        startIndex=np.min([16,startIndex])
        if rMin >= 40:
            startIndex=16
    
    if cat.scaling=="Reff" and fitInterior==True:
        endIndex=np.min([15,endIndex])
        if rMin >= 40:
            endIndex=15
    if cat.scaling=="R24" and fitInterior==True:
        endIndex=np.min([16,endIndex])
        if rMin >= 40:
            endIndex = 16
        # startIndex=np.min([16,startIndex])   
    
    
    
    

    print("startIndex = " + str(startIndex) + " (optimBetaOnly)")
    print("endIndex = " + str(endIndex) + " (optimBetaOnly)")
    
    
    
    betaVelInterp=interpFunc(Delta,beta=1,shiftFac=shiftFac)
    
    def func(beta,endIndex=endIndex):
        beta=beta
        # betaVelInterp=sp.interpolate.CubicSpline((Delta[0])+dx*shiftFac,-1/3*a*H*beta*Delta[1]*(Delta[0]+dx),bc_type = 'natural')
        
        # betaVelInterp=interpFunc(Delta,beta,shiftFac=shiftFac)
        chi2=np.sum([((beta*betaVelInterp(xpoints[i]) - trueVel[1][i])/(trueVel[2][i]))**2 for i in range(startIndex,endIndex+1)])
        print(chi2)
        return chi2
    

    optimObject=scipy.optimize.minimize(func,x0=0.35, method=None, jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=0.000001, callback=None, options=None)
    
    beta=optimObject.x[0]
    

    if plotting==True:
        
        beta=beta
        shiftFac=shiftFac
        
        betaAndShift=[beta,shiftFac]
        print("beta = " + str(beta))
        print("shiftFac = " + str(shiftFac))
        print("chi2 = " + str(func(betaAndShift)))     
        print(" ")


        betaVelInterp=interpFunc(Delta,beta,shiftFac=shiftFac)
        
        plt.errorbar(trueVel[0],trueVel[1],trueVel[2])
        plt.plot(xpoints,betaVelInterp(xpoints))
        
        plt.legend(["Linear Thy","V_r"])
        plt.title(cat.gravType + ", " + str(rMin) + " - " + str(rMax) + ", beta = " + str(int(beta*1000)/1000) + ", shiftFac = " + str(int(shiftFac*100)/100) + ", " + fitTypeText)
        
        
        if fitInterior==True:
            plt.plot([xpoints[endIndex],xpoints[endIndex]],[np.min(trueVel[1]),np.max(trueVel[1])],color="black")
        if fitExterior==True:
            plt.plot([xpoints[startIndex],xpoints[startIndex]],[np.min(trueVel[1]),np.max(trueVel[1])],color="black")
        
        
        plt.show()
        plt.clf()
        
    betaVelInterp=interpFunc(Delta,1,shiftFac=shiftFac)
    chi2=np.sum([((beta*betaVelInterp(xpoints[i]) - trueVel[1][i])/(trueVel[2][i]*10))**2 for i in range(startIndex,endIndex+1)])/(len(range(startIndex,endIndex+1))-1)
    
    
    
    print("chi2 per dof with 1 real errors is " + str(chi2))
    print("startIndex is " + str(startIndex))
    print("endIndex is " + str(endIndex))

    
    return beta,optimObject
"#########################################################################"

"chi squared gets divided by (10*np.sqrt(3))**2"

"#########################################################################"
def optimizeSigmaFixedBeta(cat,beta,startIndex,endIndex,guess=None,rMin=0,rMax=100,voidTypes=["R","S"],shiftFac=1,tol=0.05):
    "endIndex is fixed depending on the scaling of the void "

    cat.filterVoids(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
    
    #canonical choices for endIndex thus far
    
    
    if cat.scaling=="Reff":
        sigma=7.        
        # endIndex=28
        scaleByReff=True
        
        
    if cat.scaling=="R24":
        sigma=278.
        # endIndex=45
        scaleByReff=False
        
        
        
    "Construct Quad Data, if not already constructed, with JackKnife errors"    
    Q=0
    try:
        if (cat.rMinQuad==rMin and cat.rMaxQuad==rMax):
            print("Loading Quadrupole")
            quadMoment=np.array([cat.xpoints,cat.quadMoment,cat.quadErrors])
            Q=1
        if Q==0:
            print("quadrupole doesnt match required params, calculating now")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            quadMoment=cat.jackKnifeQuadMoment(rMin=rMin,rMax=rMax,voidTypes=voidTypes,numSamples=1000)
            
    except:
        print("quadrupole not pre-calculated, calculating now")
        quadMoment=cat.jackKnifeQuadMoment(rMin=rMin,rMax=rMax,voidTypes=voidTypes,numSamples=1000)
        
    
    "##################" 
    def sigmaFunc(sigma,beta=beta,quadData=quadMoment,rMin=rMin,rMax=rMax,a=1/(1+0.48251),H=130.29171508616594,voidType=["R","S"],scaleByReff=scaleByReff,startIndex=startIndex,endIndex=endIndex,nMuBins=50):


        sigmas=sigma*np.ones(len(cat.xpoints))
        lambdas=np.zeros(len(cat.xpoints))
        kappas=np.zeros(len(cat.xpoints))
        

        quadThy,deltaSandMu=calcQuadEdgeworthStreaming(cat=cat,sigmas=sigmas,lambdas=lambdas,kappas=kappas,voidTypes=voidTypes,catType="halo",rMin=rMin,rMax=rMax,modelVelocity=True,beta=beta,a=a,H=H,Integrated=True,nMuBins=nMuBins,muDependance=False,degree=0,velFac=1,startIndex=startIndex,endIndex=endIndex,shiftFac=shiftFac)
        plt.plot(quadThy[0],quadThy[1])
        plt.errorbar(quadThy[0],quadData[1],quadData[2])
        plt.title("beta = " + str(beta) + ", sigma = " + str(sigma))
        plt.show()
        plt.clf()

    
        chiSquared=calcChiSquaredForModel(modelData=quadThy[1],realData=quadData[1],dataErrors=quadData[2],startIndex=startIndex,endIndex=endIndex)
        print("chiSquared = " + str(chiSquared))
        print("beta is " + str(beta))
        print("sigma is " + str(sigma))
        
        
        
        return chiSquared
    "##################" 
        
    "##################" 
    def calcChiSquaredForModel(modelData,realData,dataErrors,startIndex=startIndex,endIndex=endIndex):
        "4, 30 or 5, 45"
        if (len(modelData) != len(realData)) and (len(dataErrors) != len(realData)):
            print("input length mismatch")
            return -1
        
        print("start index is " + str(startIndex))
        print("end index is " + str(endIndex))
        chiSquared=0
        for i in range(startIndex,endIndex):
        # for i in range(4,30):
            chiSquared += (modelData[i] - realData[i])**2/(dataErrors[i]**2)
        return chiSquared
    "##################" 
            
    optimObject=scipy.optimize.minimize(sigmaFunc,x0=sigma, method=None, jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=tol, callback=None, options=None)
    
    sigmaOpt=optimObject.x
    
    cat.optimObject=optimObject
    # cat.betaAndSigmaOpt=betaAndSigmaOpt
    print("Successfully terminated optimization!")

    return [beta,sigmaOpt,optimObject.fun]
"#########################################################################"

"#########################################################################" 
def interpFunc(Delta,beta,shiftFac=1):
    z=0.48251
    a=1./(1.+z)

    omega_M=0.3089
    omega_L=1.-omega_M
    H=100.*np.sqrt(omega_L+omega_M*(1./(a**3)))
    
    dx=(Delta[0][1]-Delta[0][0])/2
    try:
        beta=beta[0]
        
    except:
        pass
    
    
    # print(Delta)
    
    # Delta1=np.zeros(len(Delta[1]))
    # Delta1=Delta[1]
    # Delta1[0]=-1.
    
    # for i in range(1,len(Delta[0])):
    #     Delta1[i]=Delta[1][i-1]
    # betaVelInterp=sp.interpolate.CubicSpline((Delta[0]+2*dx),-1/3*a*H*beta*np.array([Delta1[i]*(Delta[0][i]+dx) for i in range(len(Delta[0]))]),bc_type = 'natural')

    print(shiftFac)
    betaVelInterp=sp.interpolate.CubicSpline((Delta[0])+(dx*(1+shiftFac)),-1/3*a*H*beta*Delta[1]*(Delta[0]+dx),bc_type = 'natural')


    return betaVelInterp
"#########################################################################" 

"#########################################################################" 
def calcQuadEdgeworthStreaming(cat,
                               sigmas,
                               lambdas=0,
                               kappas=0,
                               voidTypes=["R","S"],
                               catType="halo",
                               rMin=35,
                               rMax=100,
                               modelVelocity=False,
                               calcVelocityForPlot=False,
                               beta=1,
                               a=1/(1+0.48251)*0+1.,
                               H=130.29171508616594*0+100.,
                               Integrated=True,
                               nMuBins=50,
                               muDependance=False,
                               degree=0,
                               numBinFac=1.,
                               velFac=1,
                               startIndex=0,
                               endIndex=30,
                               shiftFac=0
                               ):
    
    """kappas and lambdas are the coefficients in front of the approapriate He(x) polynomials in the edgeworth expansions (for now)  """
    sqrt2Pi=np.sqrt(2*np.pi)
    numBins=len(cat.xpoints)
    
    
    # "DELETE THIS SOON"
    # endIndex=150
    # numBinFac=3
    
    
    
    
    
    nBins=numBins
    print(numBins)
    outTo=numBins
    xpoints=cat.xpoints
    

    if cat.scaling=="Reff":
        scaleByReff=True
        weight=-1
    if cat.scaling=="R24":
        scaleByReff=False
        weight=0
    

    
    # print("building density profile")
    if scaleByReff==True:
        
        
        
        
        
        Q=0
        try:
            if cat.rMinDens==rMin and cat.rMaxDens==rMax:
                print("loading density Profile")
                # densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                densityProf=np.array([cat.xpoints,cat.densProf,cat.densErrors])
                
                Q=1
            if Q==0:
                print("densityProf doesnt match required params, calculating now")
                densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
        except:
            print("densityProf not pre-calculated, calculating now")
            densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)

                

        # densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
        deltaInterp = sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
        # buildDensityProfilePoisson(cat,voidTypes=voidType,rMin=rMin,rMax=rMax,save=False,outTo=outTo,catType=catType)
        
        def deltaRInterp(x):
            if (x>0 and x<xpoints[-1]):
                return deltaInterp(x)
            # if (x<0):
                # return -1.
            return 0.
        
    
    if scaleByReff==False:
        
        
        Q=0
        try:
            if cat.rMinDens==rMin and cat.rMaxDens==rMax:
                print("loading density Profile")
                # densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                densityProf=np.array([cat.xpoints,cat.densProf,cat.densErrors])
                
                Q=1
            if Q==0:
                print("densityProf doesnt match required params, calculating now")
                densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
        except:
            print("densityProf not pre-calculated, calculating now")
            densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)        

        
        deltaInterp = sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
        
        def deltaRInterp(x):
            if (x>0 and x<xpoints[-1]):
                return deltaInterp(x)
            # if (x<0):
            #     return -1.
            return 0.
    
    
    

    
    if modelVelocity==False:
        # print("building velocity profile")
        if scaleByReff==True:
            
            
            

            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=-1)
            
            
            Q=0
            try:
                if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
                    print("loading Velocity Profile")
                    # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
                    velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
                    
                    
                    Q=1
                if Q==0:
                    print("velocity doesnt match required params, calculating now")
                    velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
                    
            except:
                print("velocity not pre-calculated, calculating now")
                velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
 
        
            
            
            
            
            
            
            velInterp=sp.interpolate.CubicSpline(xpoints,velProf[1]*velFac,bc_type = 'natural') # "This might have to be changed down the line to return 0 for extrapolated values"
            def velProfInterp(x):
                # print(x)
                if (x>0 or x<xpoints[-1]):
                    return velInterp(x)
                return 0.
            
        if scaleByReff==False:
            
            
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=0)
            
            
            Q=0
            try:
                if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
                    print("loading Velocity Profile")
                    # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
                    velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
                    
                    Q=1
                if Q==0:
                    print("velocity doesnt match required params, calculating now")
                    velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
                    
            except:
                print("velocity not pre-calculated, calculating now")
                velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)

        
            
            
            
            
            # velProf=buildVelocityProfileAltMethod(cat,voidTypes=voidType,rMin=rMin,rMax=rMax,weight=0,catType=catType,corFactor=1.,outTo=outTo,children=children,childless=childless,rescale=scaleByReff)
            velInterp=sp.interpolate.CubicSpline(xpoints,velProf[1]*velFac,bc_type = 'natural') # "This might have to be changed down the line to return 0 for extrapolated values"
            def velProfInterp(x):
                if (x>0 and x<xpoints[-1]):
                    return velInterp(x)
                return 0.

    if modelVelocity==True:
        if scaleByReff==True:
            # if supplyDensity == False:
                # Delta=calcIntegratedDensityContrast(cat,voidTypes=voidType,rMin=rMin,rMax=rMax,catType='halo',children=children,childless=childless,scaleByReff=scaleByReff)
                # Delta=cat.calcIntegratedDensityContrast(rMin=rMin,rMax=rMax,voidTypes=voidType)
                
            # if supplyDensity == True:
                # Delta=cat.Delta 
                
            Q=0
            try:
                if cat.rMinDelta==rMin and cat.rMaxDelta==rMax:
                    print("Loading Delta")
                    Delta=cat.Delta
                    Q=1
                if Q==0:
                    print("Delta doesnt match required params, calculating now")
                    Delta=cat.calcIntegratedDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                    
            except:
                print("Delta not pre-calculated, calculating now")
                Delta=cat.calcIntegratedDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)

                    
                

                
                
            dX=(xpoints[1]-xpoints[0])/2
            velInterp=sp.interpolate.CubicSpline(Delta[0]+dX*shiftFac,-1/3*a*H*beta*Delta[1]*(Delta[0] + dX),bc_type = 'natural')
            "TRUE VELOCITY BELOW FOR PLOTTING"
            if calcVelocityForPlot==True:
            # if True:
                
                print("Are you sure you want to calculate velocity for plotting purposes?")
                
                Q=0
                try:
                    if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
                        print("Loading Velocity")
                        # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
                        velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
                        
                        Q=1
                    if Q==0:
                        print("velocity doesnt match required params, calculating now")
                        velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
                        
                except:
                    print("velocity not pre-calculated, calculating now")
                    velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
    
            
            
                
                # velProf=buildVelocityProfileAltMethod(cat,voidTypes=voidType,rMin=rMin,rMax=rMax,weight=-1,catType=catType,corFactor=1.,outTo=outTo,children=children,childless=childless,rescale=scaleByReff)
                
                
                
                
                
                
                plt.plot(xpoints,velProf[1])
                plt.plot(xpoints,velInterp(xpoints))
            
                plt.title("beta = " + str(beta))
                plt.legend(["true velocity","Linear theory"])
                plt.show()
                plt.clf()
            
            
            def velProfInterp(x):
                # print(x)
                if (x>0 or x<xpoints[-1]):
                    return velInterp(x)
                return 0.
        
        
        if scaleByReff==False:
            

            # if supplyDensity == False:
                # Delta=calcIntegratedDensityContrast(cat,voidTypes=voidType,rMin=rMin,rMax=rMax,catType='halo',children=children,childless=childless,scaleByReff=scaleByReff)
            # if supplyDensity == True:
                # Delta=Delta
                
                
            Q=0
            try:
                if cat.rMinDelta==rMin and cat.rMaxDelta==rMax:
                    print("Loading Delta")
                    Delta=cat.Delta
                    Q=1
                if Q==0:
                    print("Delta doesnt match required params, calculating now")
                    Delta=cat.calcIntegratedDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                    
            except:
                print("Delta not pre-calculated, calculating now")
                Delta=cat.calcIntegratedDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)


            velInterp=sp.interpolate.CubicSpline((Delta[0])+shiftFac*0.05*24.,-1/3*a*H*beta*Delta[1]*(Delta[0]+0.05*24),bc_type = 'natural')
            "TRUE VELOCITY BELOW FOR PLOTTING"
            if calcVelocityForPlot==True:
            # if True:
                
                print("Are you sure you want to calculate velocity for plotting purposes?")
                # velProf=buildVelocityProfileAltMethod(cat,voidTypes=voidType,rMin=rMin,rMax=rMax,weight=0,catType=catType,corFactor=1.,outTo=outTo,children=children,childless=childless,rescale=scaleByReff)
                
                Q=0
                try:
                    if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
                        print("Loading Velocity")
                        # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
                        
                        velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
                        Q=1
                    if Q==0:
                        print("velocity doesnt match required params, calculating now")
                        velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
                        
                except:
                    print("velocity not pre-calculated, calculating now")
                    velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)

                
                
                plt.plot(xpoints,velProf[1])
                plt.plot(xpoints,velInterp(xpoints))
            
                plt.title("beta = " + str(beta))
                plt.legend(["true velocity","Linear theory"])
                plt.show()
                plt.clf()
            
            
            def velProfInterp(x):
                if (x>0 and x<xpoints[-1]):
                    return velInterp(x)
                return 0.
    "#####################################################################"



        
    "#####################################################################"
    

    if muDependance==False:
        "muDependance=False means that sigmas must be a 1xlen(xpoints) size array"
        
        #######
        for s in range(len(sigmas)):
            if sigmas[s]==0.2:
                sigmas[s]=0.05
                
        #######
        sigmaInt=sp.interpolate.CubicSpline(xpoints,sigmas,bc_type = 'natural')
        lambdaInt=sp.interpolate.CubicSpline(xpoints,lambdas,bc_type = 'natural')
        kappaInt=sp.interpolate.CubicSpline(xpoints,kappas,bc_type = 'natural')
    
        def sigmaInterp(r,mu):
            if (r>0 and r<xpoints[-1]):
                return np.max([0.05,sigmaInt(r)])
            return sigmas[-1]

    
        def lambdaInterp(r,mu):
            if (r>0 and r<xpoints[-1]):
                return lambdaInt(r)
            return lambdas[-1]
        
        
        def kappaInterp(r,mu):
            if (r>0 and r<xpoints[-1]):
                return kappaInt(r)
            return kappas[-1]
    
    
    
    
    if muDependance==True:
        "Sigmas should be 30x50 or whatever the approapriate size is"
        paramsArraySigmas=np.zeros((len(xpoints),degree+1))
        paramsArrayLambdas=np.zeros((len(xpoints),degree+1))
        paramsArrayKappas=np.zeros((len(xpoints),degree+1))
        
        
        muPointsForFit=(np.linspace(0,1,51)[1:]+np.linspace(0,1,51)[:-1])/2
        for i in range(len(sigmas)):
            
            paramsArraySigmas[i]=np.polyfit(muPointsForFit,sigmas[i],degree)
            paramsArrayLambdas[i]=np.polyfit(muPointsForFit,lambdas[i],degree)
            paramsArrayKappas[i]=np.polyfit(muPointsForFit,kappas[i],degree)
            
        
        interpArraySigmas=[sp.interpolate.CubicSpline(xpoints,[paramsArraySigmas[k][order] for k in range(len(xpoints))],bc_type = 'natural') for order in range(degree+1)]
        interpArrayLambdas=[sp.interpolate.CubicSpline(xpoints,[paramsArrayLambdas[k][order] for k in range(len(xpoints))],bc_type = 'natural') for order in range(degree+1)]
        interpArrayKappas=[sp.interpolate.CubicSpline(xpoints,[paramsArrayKappas[k][order] for k in range(len(xpoints))],bc_type = 'natural') for order in range(degree+1)]
        
        """interpArray will store all coefficients in the polynomial which have been interpolated over as a function of r
        paramsArray[firstIndex] = gives all paramaters for the polynomial at the location indicated by fI
        
        [paramsArray[i][k] for i in range(len(xpoints))] gives all of the coefficients for the kth term 
        in the polynomial across the different r values in x points. k=0 is highest order, k=degree is constand term 
        
        """
    
        def sigmaInterp(r,mu):
            mu=np.abs(mu)
            if (r>0 and r<xpoints[-1]):
                return np.abs(np.sum([interpArraySigmas[i](r)*mu**(degree-i) for i in range(degree+1)]))
                # return mInterp(r)*mu+bInterp(r)
            if (r<0):
                for i in range(100):
                    print("This should never happen!!!!!!!!!")
            return np.abs(np.sum([interpArraySigmas[i](xpoints[-1])*mu**(degree-i) for i in range(degree+1)]))
        
        
        
        def lambdaInterp(r,mu): 
            mu=np.abs(mu)
            if (r>0 and r<xpoints[-1]):
                return np.sum([interpArrayLambdas[i](r)*mu**(degree-i) for i in range(degree+1)])
                # return mInterp(r)*mu+bInterp(r)
            if (r<0):
                for i in range(100):
                    print("This should never happen!!!!!!!!!")
            return np.sum([interpArrayLambdas[i](xpoints[-1])*mu**(degree-i) for i in range(degree+1)])
        
        
        
        def kappaInterp(r,mu): 
            mu=np.abs(mu)
            if (r>0 and r<xpoints[-1]):
                return np.sum([interpArrayKappas[i](r)*mu**(degree-i) for i in range(degree+1)])
                # return mInterp(r)*mu+bInterp(r)
            if (r<0):
                for i in range(100):
                    print("This should never happen!!!!!!!!!")
            return np.sum([interpArrayKappas[i](xpoints[-1])*mu**(degree-i) for i in range(degree+1)])
        
        


            "debugging"
            if False:
                "TEST sigmaInterp"
                # for i in range(len(xpoints)):
                #     plt.plot(muPointsForFit,sigmas[i])
                #     plt.plot(muPointsForFit,[sigmaInterp(r=xpoints[i],mu=mu) for mu in muPointsForFit])
                #     plt.title("r = " + str(xpoints[i]))
                #     plt.ylabel("sigma")
                #     plt.xlabel("mu")
                #     plt.legend(["data","fit of degree " + str(degree)])
                #     plt.show()
                #     plt.clf()
                
                # for i in range(len(muPointsForFit)):
                #     plt.plot(xpoints,[sigmas[j][i] for j in range(len(xpoints))])
                #     plt.plot(xpoints,[sigmaInterp(r=x,mu=muPointsForFit[i]) for x in xpoints])
                #     plt.title("mu = " + str(muPointsForFit[i]))
                #     plt.ylabel("sigma")
                #     plt.xlabel("r/Reff")
                #     plt.legend(["data","fit of degree " + str(degree)])
                #     plt.show()
                #     plt.clf()
                "TEST sigmaInterpV2"
            
        

    ####################################################################
        
    if Integrated==False:
        def getRcoordsFromScoords(s,muS,velProf,vll):
            "function to take the s and muS coordinates, and convert to r and muR accordingly"
            sPara=s*muS
            sPerp=np.sqrt(s**2-sPara**2)
            
            rPerp=sPerp
        
        
            def rParaFunc(rPara):
                
                r=np.sqrt(rPerp**2+rPara**2)
                muR=rPara/r
                return rPara+(velProf(r)*muR+vll)/(a*H)-sPara
            
                
            rPara=fsolve(rParaFunc, x0=[sPara+vll/(a*H)])
            
            r=np.sqrt(rPerp**2+rPara**2)
            muR=rPara/r
    
            return r,muR
    
    
        def integrandOld(s,muS,vel=velProfInterp,deltaR=deltaRInterp,sigmaR=sigmaInterp,a=a,H=H,dr=0.01):
            r,muR=getRcoordsFromScoords(s,muS,velProfInterp,vll=0)
            "Have to fix the line below"
            
            Vr=vel(r)
            dvdr= (vel(r+dr)-vel(r-dr))/(2*dr)
            
            return (1+deltaR(r))/(1+Vr/(r*a*H)+(dvdr-Vr/r)/(a*H)*muR**2)
        
    
    
    # def edgeworth(x,center=0,sigma=0.5,Lambda=0,Kappa=0):
    #     X=x/sigma
    # return Gaussian(x,center=center,sigma=sigma)*(1+Lambda/6*(X**3-3*X)+Kappa/(24)*(X**4-6*X**2+3))
    

           
    
    if Integrated == True:
        def P(vll,sigma,Lambda,Kappa):
            "sigma should be sigma(r) where r is the correct value of r to correspond to vll"
            X=vll/sigma
            return 1/(sqrt2Pi*sigma)*np.e**(-(X**2)/2)*(1+Lambda/6*(X**3-3*X)+Kappa/(24)*(X**4-6*X**2+3))
        
    
    aH=a*H
    
    def integrand(x3,s,muS,vel=velProfInterp,deltaR=deltaRInterp,sigmaR=sigmaInterp,LambdaR=lambdaInterp,KappaR=kappaInterp,a=a,H=H,dr=0.01):
        
        s3=s*muS
        sPerp=np.sqrt(1-muS**2)*s
        xPerp=sPerp
        r=np.sqrt(x3**2+xPerp**2)
        muR=x3/r
        
        # r,muR=getRcoordsFromScoords(s,muS,velProfInterp,vll)
        # dvdr=(vel(r+dr)-vel(r-dr))/(2*dr)
        Vr=vel(r)
        
        
        "Have to fix the line below"
        vll=(s3-x3)*aH-muR*Vr
        # print(vll)
        return P(vll=vll,sigma=sigmaR(r,muR),Lambda=LambdaR(r,muR),Kappa=KappaR(r,muR))*(1.+deltaR(r))*aH#/(1.+Vr/(r*aH)+(dvdr-Vr/r)/(aH)*muR**2)*((aH)+dvdr*muR**2+Vr*(1.-muR**2)/r)
    

    
    
    if Integrated == True:
        def deltaS(s,muS,vel=velProfInterp,deltaR=deltaRInterp,sigmaR=sigmaInterp,a=a,H=H,dr=0.01):
        # return integrand(vll=0,s=s,muS=muS,vel=velProfInterp,deltaR=deltaRInterp,sigmaR=sigmaInterp,a=a,H=H,dr=0.01)
        
            x3Mean=s*muS
            if scaleByReff==True:
                lower=x3Mean-1.5
                upper=x3Mean+1.5
                
            if scaleByReff==False:
                lower=x3Mean-50
                upper=x3Mean+50
            
            return sp.integrate.quad(lambda x3: integrand(x3,s,muS,vel=velProfInterp,deltaR=deltaRInterp,sigmaR=sigmaInterp,a=a,H=H,dr=0.01),lower,upper,limit=200)[0]
            # return sp.integrate.quadrature(lambda vll: integrand(vll,s,muS,vel=velProfInterp,deltaR=deltaRInterp,sigmaR=sigmaInterp,a=a,H=H,dr=0.01),-cutoff,cutoff)[0]
            
            
    if Integrated == False:
        def deltaS(s,muS,vel=velProfInterp,deltaR=deltaRInterp,sigmaR=sigmaInterp,a=a,H=H,dr=0.01):
            return integrandOld(s=s,muS=muS,vel=velProfInterp,deltaR=deltaRInterp,sigmaR=sigmaInterp,a=a,H=H,dr=0.01)
            
        
    muBinEdges=np.linspace(0,1,nMuBins+1)
    
    
    dMu=muBinEdges[1]-muBinEdges[0]
    muPoints=(muBinEdges[0:-1]+muBinEdges[1:])/2
    
    
    numXPoints=int(round(numBins*numBinFac))
    

    
    if scaleByReff==True:
        nBins=numBins
        Max=nBins*0.1
        xpointsNew=(np.linspace(0,Max,numXPoints+1)[1:] + np.linspace(0,Max,numXPoints+1)[:-1])/2
    
    if scaleByReff==False:
        nBins=numBins
        Max=nBins*2.4
        xpointsNew=(np.linspace(0,Max,numXPoints+1)[1:] + np.linspace(0,Max,numXPoints+1)[:-1])/2


    deltaSandMu=np.array([[0. for j in range(len(muPoints))] for i in range(len(xpointsNew))])
    
    
    
    for i in range(startIndex,min(len(deltaSandMu),endIndex)):
        if i%10==0:
            print(i)
        for j in range(len(deltaSandMu[i])):
            #sys.stdout.flush()
            deltaSandMu[i][j]=deltaS(s=xpointsNew[i],muS=muPoints[j])
            
            
            
    
    # deltaSandMu[i][j]=np.array([[deltaS(s=xpoints[i],muS=muPoints[j]) for j in range(len(muPoints))] for i in range(len(xpoints))])

    quad=np.array([np.sum([deltaSandMu[i][j]*(5/2)*(-1+3*muPoints[j]**2)*dMu for j in range(len(muPoints))]) for i in range(len(xpointsNew))])
    
    return np.array([xpointsNew,quad]),deltaSandMu  
"#########################################################################" 

"#########################################################################" 
def heavyTheta(x):
    
    if x>=0:
        return 1
    
    if x<0:
        return 0
    
    else:
        return nan
"#########################################################################" 

"#########################################################################" 
def calcFuncDeriv_delta(cat,r,s,delta,vel,sigma,numMuPoints=100,a=a,H=H,algorithmInt=True):
    "delta, v, and sigma should already be interpolated functions"
    
    def integrand(r,s,muS,delta,vel,sigma):

        r3=np.sqrt(r**2-(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
        
        if r3==0:
            print("FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK")
        if math.isnan(r3)==True:
            print("had a nan")
            r3=0.001
            
        
        
        
        
        
        # print(r3)
        Ap=(r3-muS*s)*a*H+r3/r*vel(r)
        Am=(-r3-muS*s)*a*H-r3/r*vel(r)
        
        Express=5*a*H/4*(3*muS**2-1)*(r/r3)*(1/(sqrt2pi*sigma(r))*(np.exp((-Ap**2)/(2*sigma(r)**2)) + np.exp((-Am**2)/(2*sigma(r)**2))))
        
        return Express
    
    
        
    if algorithmInt==False:
    
        if r>s:
            muPoints=(np.linspace(-1,1,numMuPoints+1)[1:] + np.linspace(-1,1,numMuPoints+1)[:-1])/2
            
            dMu=muPoints[1]-muPoints[0]
            
            funcVals=np.zeros(len(muPoints))
            
            for i in range(len(funcVals)):
                funcVals[i]=integrand(r=r,s=s,muS=muPoints[i],delta=delta,vel=vel,sigma=sigma)
                val=np.sum(funcVals)*dMu
                
                
                
                
            plt.plot(muPoints,funcVals)
                
            return val
                
        
        
        
        
        if r<=s:
                
            numMuPoints1=int(numMuPoints/2)+1
            muPoints1=(np.linspace(-1,-np.sqrt(1-(r/s)**2),numMuPoints1+1)[1:]+ np.linspace(-1,-np.sqrt(1-(r/s)**2),numMuPoints1+1)[:-1])/2
                       
            
            dMu1=muPoints1[1]-muPoints1[0]
            
            funcVals1=np.zeros(len(muPoints1))
            
            for i in range(len(funcVals1)):
                funcVals1[i]=integrand(r=r,s=s,muS=muPoints1[i],delta=delta,vel=vel,sigma=sigma)
                val1=np.sum(funcVals1)*dMu1
            plt.plot(muPoints1,funcVals1)
            
            
            numMuPoints2=int(numMuPoints/2)
            muPoints2=(np.linspace(np.sqrt(1-(r/s)**2),1,numMuPoints2+1)[1:]+ np.linspace(np.sqrt(1-(r/s)**2),1,numMuPoints2+1)[:-1])/2
            
            dMu2=muPoints2[1]-muPoints2[0]
            
            funcVals2=np.zeros(len(muPoints2))
            
            for i in range(len(funcVals2)):
                funcVals2[i]=integrand(r=r,s=s,muS=muPoints2[i],delta=delta,vel=vel,sigma=sigma)
                val2=np.sum(funcVals2)*dMu2
                
                
            plt.plot(muPoints2,funcVals2)
            val=val1+val2
            
            return val

                
                
                
                
    if algorithmInt==True:
        
        if r>s:
            val=sp.integrate.quad(lambda mu: integrand(r=r,s=s,muS=mu ,delta=delta,vel=vel,sigma=sigma),-1,1,limit=200)[0]
            
            return val
                
    
        if r<=s:
            epsilon=1e-8
            
            
            val1=sp.integrate.quad(lambda mu: integrand(r=r,s=s,muS=mu ,delta=delta,vel=vel,sigma=sigma),-1,-np.sqrt(1-(r/s)**2)-epsilon,limit=200)[0]
            
            val2=sp.integrate.quad(lambda mu: integrand(r=r,s=s,muS=mu ,delta=delta,vel=vel,sigma=sigma),np.sqrt(1-(r/s)**2)+epsilon,1,limit=200)[0]
            
            val=val1+val2
            
            return val
"#########################################################################"

"#########################################################################" 
def calcFuncDeriv_vel(cat,r,s,delta,vel,sigma,numMuPoints=100,a=a,H=H,algorithmInt=True):
    "delta, v, and sigma should already be interpolated functions"
    
    def integrand(r,s,muS,delta,vel,sigma):

        r3=np.sqrt(r**2-(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
        
        if r3==0:
            print("FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK")
        if math.isnan(r3)==True:
            print("had a nan")
            r3=0.001
            
        
        
        
        
        
        # print(r3)
        Ap=(r3-muS*s)*a*H+r3/r*vel(r)
        Am=(-r3-muS*s)*a*H-r3/r*vel(r)
        
        Express=5*a*H/4*(3*muS**2-1)*(1)*(-1/(sqrt2pi*sigma(r))*(np.exp((-Ap**2)/(2*sigma(r)**2))*(Ap/sigma(r)**2)*(1) + np.exp((-Am**2)/(2*sigma(r)**2))*((Am/sigma(r)**2)*(-1))))*(1+delta(r))
        
        return Express
    
    
        
    if algorithmInt==False:
    
        if r>s:
            muPoints=(np.linspace(-1,1,numMuPoints+1)[1:] + np.linspace(-1,1,numMuPoints+1)[:-1])/2
            
            dMu=muPoints[1]-muPoints[0]
            
            funcVals=np.zeros(len(muPoints))
            
            for i in range(len(funcVals)):
                funcVals[i]=integrand(r=r,s=s,muS=muPoints[i],delta=delta,vel=vel,sigma=sigma)
                val=np.sum(funcVals)*dMu
                
                
                
                
            plt.plot(muPoints,funcVals)
                
            return val
                
        
        
        
        
        if r<=s:
                
            numMuPoints1=int(numMuPoints/2)+1
            muPoints1=(np.linspace(-1,-np.sqrt(1-(r/s)**2),numMuPoints1+1)[1:]+ np.linspace(-1,-np.sqrt(1-(r/s)**2),numMuPoints1+1)[:-1])/2
                       
            
            dMu1=muPoints1[1]-muPoints1[0]
            
            funcVals1=np.zeros(len(muPoints1))
            
            for i in range(len(funcVals1)):
                funcVals1[i]=integrand(r=r,s=s,muS=muPoints1[i],delta=delta,vel=vel,sigma=sigma)
                val1=np.sum(funcVals1)*dMu1
            plt.plot(muPoints1,funcVals1)
            
            
            numMuPoints2=int(numMuPoints/2)
            muPoints2=(np.linspace(np.sqrt(1-(r/s)**2),1,numMuPoints2+1)[1:]+ np.linspace(np.sqrt(1-(r/s)**2),1,numMuPoints2+1)[:-1])/2
            
            dMu2=muPoints2[1]-muPoints2[0]
            
            funcVals2=np.zeros(len(muPoints2))
            
            for i in range(len(funcVals2)):
                funcVals2[i]=integrand(r=r,s=s,muS=muPoints2[i],delta=delta,vel=vel,sigma=sigma)
                val2=np.sum(funcVals2)*dMu2
                
                
            plt.plot(muPoints2,funcVals2)
            val=val1+val2
            
            return val

                
                
                
                
    if algorithmInt==True:
        
        if r>s:
            val=sp.integrate.quad(lambda mu: integrand(r=r,s=s,muS=mu ,delta=delta,vel=vel,sigma=sigma),-1,1,limit=200)[0]
            
            return val
                
    
        if r<=s:
            epsilon=1e-8
            
            
            val1=sp.integrate.quad(lambda mu: integrand(r=r,s=s,muS=mu ,delta=delta,vel=vel,sigma=sigma),-1,-np.sqrt(1-(r/s)**2)-epsilon,limit=200)[0]
            
            val2=sp.integrate.quad(lambda mu: integrand(r=r,s=s,muS=mu ,delta=delta,vel=vel,sigma=sigma),np.sqrt(1-(r/s)**2)+epsilon,1,limit=200)[0]
            
            val=val1+val2
            
            return val
"#########################################################################"

"#########################################################################" 
def calcFuncDeriv_sigma(cat,r,s,delta,vel,sigma,numMuPoints=100,a=a,H=H,algorithmInt=True):
    "delta, v, and sigma should already be interpolated functions"
    
    def integrand(r,s,muS,delta,vel,sigma):

        r3=np.sqrt(r**2-(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
        
        if r3==0:
            print("FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK FUCK")
        if math.isnan(r3)==True:
            print("had a nan")
            r3=0.001
            
        
        
        
        
        
        # print(r3)
        Ap=(r3-muS*s)*a*H+r3/r*vel(r)
        Am=(-r3-muS*s)*a*H-r3/r*vel(r)
        
        Express=5*a*H/4*(3*muS**2-1)*(r/r3)*(     (-1/(sqrt2pi*sigma(r)**2)*(np.exp((-Ap**2)/(2*sigma(r)**2)) + np.exp((-Am**2)/(2*sigma(r)**2)))) +   (1/(sqrt2pi*sigma(r))*(np.exp((-Ap**2)/(2*sigma(r)**2))*(Ap**2/sigma(r)**3)+ np.exp((-Am**2)/(2*sigma(r)**2))*(Am**2/sigma(r)**3))))*(1+delta(r))
        
        return Express
    
    
        
    if algorithmInt==False:
    
        if r>s:
            muPoints=(np.linspace(-1,1,numMuPoints+1)[1:] + np.linspace(-1,1,numMuPoints+1)[:-1])/2
            
            dMu=muPoints[1]-muPoints[0]
            
            funcVals=np.zeros(len(muPoints))
            
            for i in range(len(funcVals)):
                funcVals[i]=integrand(r=r,s=s,muS=muPoints[i],delta=delta,vel=vel,sigma=sigma)
                val=np.sum(funcVals)*dMu
                
                
                
                
            plt.plot(muPoints,funcVals)
                
            return val
                
        
        
        
        
        if r<=s:
                
            numMuPoints1=int(numMuPoints/2)+1
            muPoints1=(np.linspace(-1,-np.sqrt(1-(r/s)**2),numMuPoints1+1)[1:]+ np.linspace(-1,-np.sqrt(1-(r/s)**2),numMuPoints1+1)[:-1])/2
                       
            
            dMu1=muPoints1[1]-muPoints1[0]
            
            funcVals1=np.zeros(len(muPoints1))
            
            for i in range(len(funcVals1)):
                funcVals1[i]=integrand(r=r,s=s,muS=muPoints1[i],delta=delta,vel=vel,sigma=sigma)
                val1=np.sum(funcVals1)*dMu1
            plt.plot(muPoints1,funcVals1)
            
            
            numMuPoints2=int(numMuPoints/2)
            muPoints2=(np.linspace(np.sqrt(1-(r/s)**2),1,numMuPoints2+1)[1:]+ np.linspace(np.sqrt(1-(r/s)**2),1,numMuPoints2+1)[:-1])/2
            
            dMu2=muPoints2[1]-muPoints2[0]
            
            funcVals2=np.zeros(len(muPoints2))
            
            for i in range(len(funcVals2)):
                funcVals2[i]=integrand(r=r,s=s,muS=muPoints2[i],delta=delta,vel=vel,sigma=sigma)
                val2=np.sum(funcVals2)*dMu2
                
                
            plt.plot(muPoints2,funcVals2)
            val=val1+val2
            
            return val

                
                
                
                
    if algorithmInt==True:
        
        if r>s:
            val=sp.integrate.quad(lambda mu: integrand(r=r,s=s,muS=mu ,delta=delta,vel=vel,sigma=sigma),-1,1,limit=200)[0]
            
            return val
                
    
        if r<=s:
            epsilon=1e-6
            
            
            val1=sp.integrate.quad(lambda mu: integrand(r=r,s=s,muS=mu ,delta=delta,vel=vel,sigma=sigma),-1,-np.sqrt(1-(r/s)**2)-epsilon,limit=200)[0]
            
            val2=sp.integrate.quad(lambda mu: integrand(r=r,s=s,muS=mu ,delta=delta,vel=vel,sigma=sigma),np.sqrt(1-(r/s)**2)+epsilon,1,limit=200)[0]
            
            val=val1+val2
            
            return val
"#########################################################################"

"#########################################################################"
def calcDeltaXi2_delta(cat,cat2,spoints,rMin=25,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=101,voidTypes=["R","S"]):
    
       
    xpoints=cat.xpoints
    

    if cat.scaling=="Reff":
        weight=-1
    if cat.scaling=="R24":
        weight=0
    
    xpoints=cat.xpoints
    
    

    Q=0
    try:
        if cat.rMinDens==rMin and cat.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf=np.array([cat.xpoints,cat.densProf,cat.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    Q=0
    try:
        if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    

    Q=0
    try:
        if cat.rMinSigma==rMin and cat.rMaxSigma==rMax and cat.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf=cat.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
        

    Q=0
    


    Q=0
    try:
        if cat2.rMinDens==rMin and cat2.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf2=np.array([cat2.xpoints,cat2.densProf,cat2.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf2 = cat2.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf2 = cat2.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    





    delta=sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
    vel=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')
    sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf,bc_type = 'natural')
    delta_delta=sp.interpolate.CubicSpline(xpoints,densityProf2[1]-densityProf[1],bc_type = 'natural')
    

    
    
    
    
    
    
    "delta, v, sigma, and delta_var should already be interpolated functions"
    
    def velInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return vel(x)
        return 0.
    
    def deltaInterp(x):
        
        if x<0.2:
            return -1
        
        if (x>0.2 or x<xpoints[-1]):
            return delta(x)
        return 0.
    
    def sigmaInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return np.max([0.2,sigma(x)])
        return sigma(xpoints[-1])
    
    
    def delta_deltaInterp(x):
        if x<0.2:
            return 0
        if (x>0.2 or x<xpoints[-1]):
            return delta_delta(x)
        return 0
    
    
    muPoints=(np.linspace(0,1,numMuPoints+1)[:-1] + np.linspace(0,1,numMuPoints+1)[1:])/2 
    muIntList=np.zeros(len(muPoints))
    dMu=muPoints[1]-muPoints[0]
    
    
    
    
    def integrand(r3,s,muS):

        r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
    

        # print(r3)
        Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)

        
        Express=5*a*H/4*(3*muS**2-1)*(1/(sqrt2pi*sigmaInterp(r))*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))))*delta_deltaInterp(r)
        return Express
    
        
    if cat.scaling=="Reff":
        r3interval=1.5
        
    if cat.scaling=="R24":
        r3interval=50

    
    deltaXi2points=np.zeros(len(spoints))
    for i in range(len(spoints)):
        print(i)
        s=spoints[i]
        for j in range(len(muPoints)):
            muS=muPoints[j]
            
            r3points= np.linspace(s*muS-r3interval,s*muS+r3interval,200)
            # expressPoints=[integrand(r3,s,muS) for r3 in r3points]
            
            # plt.plot(r3points,expressPoints)
            # plt.show()
            # plt.clf()
            
            
            muIntList[j]=2*sp.integrate.quad(lambda r3: integrand(r3,s,muS),s*muS-r3interval,s*muS+r3interval,limit=200)[0] 
            
        deltaXi2points[i]=np.sum(muIntList)*dMu
        
            
        
        
    return spoints,deltaXi2points
"#########################################################################"    

"#########################################################################"
def calcDeltaXi2_vel(cat,cat2,spoints,rMin=25,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=101,voidTypes=["R","S"],otherMoments=False):
    
    
    
    xpoints=cat.xpoints
    

    if cat.scaling=="Reff":
        weight=-1
    if cat.scaling=="R24":
        weight=0
    
    xpoints=cat.xpoints
    
    

    Q=0
    try:
        if cat.rMinDens==rMin and cat.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf=np.array([cat.xpoints,cat.densProf,cat.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    Q=0
    try:
        if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    

    Q=0
    try:
        if cat.rMinSigma==rMin and cat.rMaxSigma==rMax and cat.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf=cat.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
        

    Q=0
    
    try:
        if (cat2.rMinVels==rMin and cat2.rMaxVels==rMax and cat2.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat2.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf2=np.array([cat2.xpoints,cat2.velProf,cat2.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf2=cat2.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf2=cat2.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        

    delta=sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
    vel=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')
    sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf,bc_type = 'natural')
    delta_vel=sp.interpolate.CubicSpline(xpoints,velProf2[1]-velProf[1],bc_type = 'natural')
    

    
    "delta, v, sigma, and delta_var should already be interpolated functions"


    def velInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return vel(x)
        return 0.
    
    def deltaInterp(x):
        
        if x<0.2:
            return -1
        
        if (x>0.2 or x<xpoints[-1]):
            return delta(x)
        return 0.
    
    def sigmaInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return np.max([0.2,sigma(x)])
        return sigma(xpoints[-1])
    
    
    def delta_velInterp(x):
        if x<0.2:
            return 0
        if (x>0.2 or x<xpoints[-1]):
            return delta_vel(x)
        return 0
    
    
    muPoints=(np.linspace(0,1,numMuPoints+1)[:-1] + np.linspace(0,1,numMuPoints+1)[1:])/2 
    muIntList=np.zeros(len(muPoints))
    dMu=muPoints[1]-muPoints[0]
    
    
    
    
    def integrand(r3,s,muS):

        r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
    

        # print(r3)
        Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)

        
        Express=5*a*H/4*(3*muS**2-1)*(r3/r)*(-1/(sqrt2pi*sigmaInterp(r)**3)*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))*Ap ))*(1+deltaInterp(r))*delta_velInterp(r)
        
        return Express
    
        
    if cat.scaling=="Reff":
        r3interval=1.5
        
    if cat.scaling=="R24":
        r3interval=50

    
    deltaXi2points=np.zeros(len(spoints))
    for i in range(len(spoints)):
        print(i)
        s=spoints[i]
        for j in range(len(muPoints)):
            muS=muPoints[j]
            
            r3points= np.linspace(s*muS-r3interval,s*muS+r3interval,200)
            # expressPoints=[integrand(r3,s,muS,delta,vel,sigma,delta_vel) for r3 in r3points]
            
            # plt.plot(r3points,expressPoints)
            # plt.show()
            # plt.clf()
            
            
            muIntList[j]=2*sp.integrate.quad(lambda r3: integrand(r3,s,muS),s*muS-r3interval,s*muS+r3interval,limit=200)[0] 
            
        deltaXi2points[i]=np.sum(muIntList)*dMu
        
            
        
        
    return spoints,deltaXi2points
"#########################################################################"    

"#########################################################################"
def calcDeltaXi2_sigma(cat,cat2,spoints,rMin=25,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=101,voidTypes=["R","S"]):
 
        
    xpoints=cat.xpoints
    

    if cat.scaling=="Reff":
        weight=-1
    if cat.scaling=="R24":
        weight=0
    
    xpoints=cat.xpoints
    
    

    Q=0
    try:
        if cat.rMinDens==rMin and cat.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf=np.array([cat.xpoints,cat.densProf,cat.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    Q=0
    try:
        if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    

    Q=0
    try:
        if cat.rMinSigma==rMin and cat.rMaxSigma==rMax and cat.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf=cat.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
        

    Q=0

    try:
        if cat2.rMinSigma==rMin and cat2.rMaxSigma==rMax and cat2.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf2=cat2.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf2 = cat2.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf2 = cat2.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
            

    delta=sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
    vel=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')
    sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf,bc_type = 'natural')
    delta_sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf2-sigmaProf,bc_type = 'natural')
    
 
    
    "delta, v, sigma, and delta_var should already be interpolated functions"
    xpoints=cat.xpoints
    

    def velInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return vel(x)
        return 0.
    
    def deltaInterp(x):
        
        if x<0.2:
            return -1
        
        if (x>0.2 or x<xpoints[-1]):
            return delta(x)
        return 0.
    
    def sigmaInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return np.max([0.2,sigma(x)])
        return sigma(xpoints[-1])
    
    
    def delta_sigmaInterp(x):
        if x<0.2:
            return 0
        if (x>0.2 or x<xpoints[-1]):
            return delta_sigma(x)
        return delta_sigma(xpoints[-1])
    
    
    muPoints=(np.linspace(0,1,numMuPoints+1)[:-1] + np.linspace(0,1,numMuPoints+1)[1:])/2 
    muIntList=np.zeros(len(muPoints))
    dMu=muPoints[1]-muPoints[0]
    
    
    
    
    def integrand(r3,s,muS):

        r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
    

        # print(r3)
        Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)

        
        Express=5*a*H/4*(3*muS**2-1)*((-1/(sqrt2pi*sigmaInterp(r)**2)*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2)) )) + (1/(sqrt2pi*sigmaInterp(r)))*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))*(Ap**2/(sigmaInterp(r)**3))))*(1+deltaInterp(r))*delta_sigmaInterp(r)
        
        return Express
    
        
    if cat.scaling=="Reff":
        r3interval=0.75
        
    if cat.scaling=="R24":
        r3interval=50

    
    deltaXi2points=np.zeros(len(spoints))
    for i in range(len(spoints)):
        print(i)
        s=spoints[i]
        for j in range(len(muPoints)):
            muS=muPoints[j]
            
            # r3points= np.linspace(s*muS-r3interval,s*muS+r3interval,200)
            # expressPoints=[integrand(r3,s,muS,delta,vel,sigma,delta_sigma) for r3 in r3points]
            
            # plt.plot(r3points,expressPoints)
            # plt.show()
            # plt.clf()
            
            
            muIntList[j]=2*sp.integrate.quad(lambda r3: integrand(r3,s,muS),s*muS-r3interval,s*muS+r3interval,limit=200)[0] 
            
        deltaXi2points[i]=np.sum(muIntList)*dMu
        
            
        
        
    return spoints,deltaXi2points
"#########################################################################"    

"#########################################################################"
def calcDeltaXi2DiffScaling_delta(catR24,catReff,spoints,rMin=35,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=101,voidTypes=["R","S"],returnDensProfs=False,p=3):
    
       
    """take R24 as the base case, convert Reff quantities to R24 ones through the correct multiplication by R_{eff}, and then take
    Reff-R24 as \Delta quantity, and compute functional derivatives"""
    
    
    
    
    xpoints=catR24.xpoints
    
    "Return to this later"
    # if cat.scaling=="Reff":
    #     weight=-1
    # if cat.scaling=="R24":
    #     weight=0
    
    weightR24=0
    weightReff=-1
    
    xpointsR24=catR24.xpoints
    xpointsReff=catReff.xpoints
    
    

    Q=0
    try:
        if catR24.rMinDens==rMin and catR24.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf=np.array([catR24.xpoints,catR24.densProf,catR24.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf = catR24.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf = catR24.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    Q=0
    try:
        if (catR24.rMinVels==rMin and catR24.rMaxVels==rMax and catR24.velWeight==weightR24):
            # print("loading Velocity Profile")
            # velProf=catR24.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([catR24.xpoints,catR24.velProf,catR24.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=catR24.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightR24)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=catR24.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightR24)
        
    

    Q=0
    try:
        if catR24.rMinSigma==rMin and catR24.rMaxSigma==rMax and catR24.sigmaWeight==weightR24:
            # print("loading sigma profile")
            sigmaProf=catR24.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf = catR24.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightR24)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf = catR24.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightR24)   
        
        

    Q=0
    


    Q=0
    try:
        if catReff.rMinDens==rMin and catReff.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf2=np.array([catReff.xpoints,catReff.densProf,catReff.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf2 = catReff.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf2 = catReff.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    "Here, we convert Reff quantities to R24 quantities"
    
    Reff=(np.sum([catReff.voidsToUse[i].radius**p for i in range(len(catReff.voidsToUse))])/len(catReff.voidsToUse))**(1/p)

    "Convert densityProf2, which initially came from Reff voids, to a R24 quantity"
    densityProfReffInterp=sp.interpolate.CubicSpline(xpointsReff*Reff,densityProf2[1],bc_type = 'natural')

    densityProf2=np.array([xpointsR24,densityProfReffInterp(xpointsR24)])

    "Proceede as before"
    
    # "TEST - DELETE THIS WHEN DONE"
    # if True:
    #     densityProf2[1][8:]=densityProf[1][8:]
    
    
    
    
    delta=sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
    vel=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')
    sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf,bc_type = 'natural')
    delta_delta=sp.interpolate.CubicSpline(xpoints,densityProf2[1]-densityProf[1],bc_type = 'natural')
    

    
    plt.plot(densityProf[0],densityProf[1])
    plt.plot(densityProf2[0],densityProf2[1])
    plt.show()
    plt.clf()
    
    plt.plot(densityProf[0],densityProf2[1] - densityProf[1])
    # plt.plot(densityProf2[0],densityProf2[1])
    plt.show()
    plt.clf()
    
    
    
    
    "delta, v, sigma, and delta_var should already be interpolated functions"
    
    def velInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return vel(x)
        return 0.
    
    def deltaInterp(x):
        
        if x<0.2:
            return -1
        
        if (x>0.2 or x<xpoints[-1]):
            return delta(x)
        return 0.
    
    def sigmaInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return np.max([0.2,sigma(x)])
        return sigma(xpoints[-1])
    
    
    def delta_deltaInterp(x):
        if x<0.2:
            return 0
        if (x>0.2 or x<xpoints[-1]):
            return delta_delta(x)
        return 0
    
    
    muPoints=(np.linspace(0,1,numMuPoints+1)[:-1] + np.linspace(0,1,numMuPoints+1)[1:])/2 
    muIntList=np.zeros(len(muPoints))
    dMu=muPoints[1]-muPoints[0]
    
    
    
    
    def integrand(r3,s,muS):

        r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
    

        # print(r3)
        Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)

        
        Express=5*a*H/4*(3*muS**2-1)*(1/(sqrt2pi*sigmaInterp(r))*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))))*delta_deltaInterp(r)
        return Express
    
        
    # if cat.scaling=="Reff":
    #     r3interval=1.5
        
    # if cat.scaling=="R24":
    r3interval=50

    
    deltaXi2points=np.zeros(len(spoints))
    for i in range(len(spoints)):
        print(i)
        s=spoints[i]
        for j in range(len(muPoints)):
            muS=muPoints[j]
            
            r3points= np.linspace(s*muS-r3interval,s*muS+r3interval,200)
            # expressPoints=[integrand(r3,s,muS) for r3 in r3points]
            
            # plt.plot(r3points,expressPoints)
            # plt.show()
            # plt.clf()
            
            
            muIntList[j]=2*sp.integrate.quad(lambda r3: integrand(r3,s,muS),s*muS-r3interval,s*muS+r3interval,limit=200)[0] 
            
        deltaXi2points[i]=np.sum(muIntList)*dMu
        
            
        
    if returnDensProfs==True:
        return spoints,deltaXi2points,densityProf,densityProf2
    
    return spoints,deltaXi2points
"#########################################################################"    

"#########################################################################"
def calcDeltaXi2DiffScaling_vel(catR24,catReff,spoints,rMin=35,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=101,voidTypes=["R","S"],otherMoments=False,p=3):
    
    
    
    xpointsR24=catR24.xpoints
    xpoints=catR24.xpoints
    xpointsReff=catReff.xpoints
    weightR24=0
    weightReff=-1

    # if cat.scaling=="Reff":
    #     weight=-1
    # if cat.scaling=="R24":
    #     weight=0
    

    
    

    Q=0
    try:
        if catR24.rMinDens==rMin and catR24.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf=np.array([catR24.xpoints,catR24.densProf,catR24.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf = catR24.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf = catR24.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    Q=0
    try:
        if (catR24.rMinVels==rMin and catR24.rMaxVels==rMax and catR24.velWeight==weightR24):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([catR24.xpoints,catR24.velProf,catR24.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=catR24.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightR24)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=catR24.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightR24)
        
    

    Q=0
    try:
        if catR24.rMinSigma==rMin and catR24.rMaxSigma==rMax and catR24.sigmaWeight==weightR24:
            # print("loading sigma profile")
            sigmaProf=catR24.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf = catR24.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightR24)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf = catR24.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightR24)   
        
        

    Q=0
    
    try:
        if (catReff.rMinVels==rMin and catReff.rMaxVels==rMax and catReff.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=catReff.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf2=np.array([catReff.xpoints,catReff.velProf,catReff.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf2=catReff.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightReff)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf2=catReff.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightReff)
        




    "Here, we convert Reff quantities to R24 quantities"
    
    Reff=(np.sum([catReff.voidsToUse[i].radius**p for i in range(len(catReff.voidsToUse))])/len(catReff.voidsToUse))**(1/p)

    "Convert velProf2, which initially came from Reff voids, to a R24 quantity"
    velProfReffInterp=sp.interpolate.CubicSpline(xpointsReff*Reff,velProf2[1]*Reff,bc_type = 'natural')

    velProf2=np.array([xpointsR24,velProfReffInterp(xpointsR24)])



    delta=sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
    vel=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')
    sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf,bc_type = 'natural')
    
    
    delta_vel=sp.interpolate.CubicSpline(xpoints,velProf2[1]-velProf[1],bc_type = 'natural')
    
    
    "delta, v, sigma, and delta_var should already be interpolated functions"

    plt.plot(velProf[0],velProf[1])
    plt.plot(velProf2[0],velProf2[1])
    plt.show()
    plt.clf()
    
    plt.plot(velProf[0],velProf[1])
    plt.plot(velProf2[0],velProf2[1])
    plt.show()
    plt.clf()
    
    

    def velInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return vel(x)
        return 0.
    
    def deltaInterp(x):
        
        if x<0.2:
            return -1
        
        if (x>0.2 or x<xpoints[-1]):
            return delta(x)
        return 0.
    
    def sigmaInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return np.max([0.2,sigma(x)])
        return sigma(xpoints[-1])
    
    
    def delta_velInterp(x):
        if x<0.2:
            return 0
        if (x>0.2 or x<xpoints[-1]):
            return delta_vel(x)
        return 0
    
    
    muPoints=(np.linspace(0,1,numMuPoints+1)[:-1] + np.linspace(0,1,numMuPoints+1)[1:])/2 
    muIntList=np.zeros(len(muPoints))
    dMu=muPoints[1]-muPoints[0]
    
    
    
    
    def integrand(r3,s,muS):

        r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
    

        # print(r3)
        Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)

        
        Express=5*a*H/4*(3*muS**2-1)*(r3/r)*(-1/(sqrt2pi*sigmaInterp(r)**3)*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))*Ap ))*(1+deltaInterp(r))*delta_velInterp(r)
        
        return Express
    
        
    # if cat.scaling=="Reff":
    #     r3interval=1.5
        
    # if cat.scaling=="R24":
    r3interval=50

    
    deltaXi2points=np.zeros(len(spoints))
    for i in range(len(spoints)):
        print(i)
        s=spoints[i]
        for j in range(len(muPoints)):
            muS=muPoints[j]
            
            r3points= np.linspace(s*muS-r3interval,s*muS+r3interval,200)
            # expressPoints=[integrand(r3,s,muS,delta,vel,sigma,delta_vel) for r3 in r3points]
            
            # plt.plot(r3points,expressPoints)
            # plt.show()
            # plt.clf()
            
            
            muIntList[j]=2*sp.integrate.quad(lambda r3: integrand(r3,s,muS),s*muS-r3interval,s*muS+r3interval,limit=200)[0] 
            
        deltaXi2points[i]=np.sum(muIntList)*dMu
        
            
        
        
    return spoints,deltaXi2points
"#########################################################################"    

"#########################################################################"
def calcDeltaXi2DiffScaling_sigma(catR24,catReff,spoints,rMin=25,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=101,voidTypes=["R","S"],p=3):
 
        
    xpointsReff=catReff.xpoints
    xpoints=catR24.xpoints
    
    weightReff=-1
    weightR24=0
    
    xpointsR24=catR24.xpoints
    
    

    Q=0
    try:
        if catR24.rMinDens==rMin and catR24.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf=np.array([catR24.xpoints,catR24.densProf,catR24.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf = catR24.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf = catR24.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    Q=0
    try:
        if (catR24.rMinVels==rMin and catR24.rMaxVels==rMax and catR24.velWeight==weightR24):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([catR24.xpoints,catR24.velProf,catR24.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=catR24.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightR24)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=catR24.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightR24)
        
    

    Q=0
    try:
        if catR24.rMinSigma==rMin and catR24.rMaxSigma==rMax and catR24.sigmaWeight==weightR24:
            # print("loading sigma profile")
            sigmaProf=catR24.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf = catR24.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightR24)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf = catR24.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightR24)   
        
        

    Q=0

    try:
        if catReff.rMinSigma==rMin and catReff.rMaxSigma==rMax and catReff.sigmaWeight==weightReff:
            # print("loading sigma profile")
            sigmaProf2=catReff.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf2 = catReff.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightReff)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf2 = catReff.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weightReff)   
        
            
        
        
        
        
    "Here, we convert Reff quantities to R24 quantities"
    
    Reff=(np.sum([catReff.voidsToUse[i].radius**p for i in range(len(catReff.voidsToUse))])/len(catReff.voidsToUse))**(1/p)

    "Convert velProf2, which initially came from Reff voids, to a R24 quantity"
    sigmaProfReffInterp=sp.interpolate.CubicSpline(xpointsReff*Reff,sigmaProf2*Reff,bc_type = 'natural')

    sigmaProf2=sigmaProfReffInterp(xpointsR24)


    
        
        
        
        

    delta=sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
    vel=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')
    sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf,bc_type = 'natural')
    delta_sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf2-sigmaProf,bc_type = 'natural')
    
 
    
 
    
 
    
 
    
 
    "delta, v, sigma, and delta_var should already be interpolated functions"
    # xpoints=cat.xpoints
    

    def velInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return vel(x)
        return 0.
    
    def deltaInterp(x):
        
        if x<0.2:
            return -1
        
        if (x>0.2 or x<xpoints[-1]):
            return delta(x)
        return 0.
    
    def sigmaInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return np.max([0.2,sigma(x)])
        return sigma(xpoints[-1])
    
    
    def delta_sigmaInterp(x):
        if x<0.2:
            return 0
        if (x>0.2 or x<xpoints[-1]):
            return delta_sigma(x)
        return delta_sigma(xpoints[-1])
    
    
    muPoints=(np.linspace(0,1,numMuPoints+1)[:-1] + np.linspace(0,1,numMuPoints+1)[1:])/2 
    muIntList=np.zeros(len(muPoints))
    dMu=muPoints[1]-muPoints[0]
    
    
    
    
    def integrand(r3,s,muS):

        r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
    

        # print(r3)
        Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)

        
        Express=5*a*H/4*(3*muS**2-1)*((-1/(sqrt2pi*sigmaInterp(r)**2)*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2)) )) + (1/(sqrt2pi*sigmaInterp(r)))*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))*(Ap**2/(sigmaInterp(r)**3))))*(1+deltaInterp(r))*delta_sigmaInterp(r)
        
        return Express
    
        

    r3interval=50

    
    deltaXi2points=np.zeros(len(spoints))
    for i in range(len(spoints)):
        print(i)
        s=spoints[i]
        for j in range(len(muPoints)):
            muS=muPoints[j]
            
            # r3points= np.linspace(s*muS-r3interval,s*muS+r3interval,200)
            # expressPoints=[integrand(r3,s,muS,delta,vel,sigma,delta_sigma) for r3 in r3points]
            
            # plt.plot(r3points,expressPoints)
            # plt.show()
            # plt.clf()
            
            
            muIntList[j]=2*sp.integrate.quad(lambda r3: integrand(r3,s,muS),s*muS-r3interval,s*muS+r3interval,limit=200)[0] 
            
        deltaXi2points[i]=np.sum(muIntList)*dMu
        
            
        
        
    return spoints,deltaXi2points
"#########################################################################"    

"#########################################################################"
def calcDeltaMoments_delta(cat,cat2,spoints,rMin=35,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=101,voidTypes=["R","S"]):
    
    
    
    xpoints=cat.xpoints
    

    if cat.scaling=="Reff":
        weight=-1
    if cat.scaling=="R24":
        weight=0
    
    xpoints=cat.xpoints
    
    

    Q=0
    try:
        if cat.rMinDens==rMin and cat.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf=np.array([cat.xpoints,cat.densProf,cat.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    Q=0
    try:
        if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    

    Q=0
    try:
        if cat.rMinSigma==rMin and cat.rMaxSigma==rMax and cat.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf=cat.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
        

    Q=0
    


    Q=0
    try:
        if cat2.rMinDens==rMin and cat2.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf2=np.array([cat2.xpoints,cat2.densProf,cat2.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf2 = cat2.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf2 = cat2.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    



    delta=sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
    vel=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')
    sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf,bc_type = 'natural')
    delta_delta=sp.interpolate.CubicSpline(xpoints,densityProf2[1]-densityProf[1],bc_type = 'natural')
    

    
    "delta, v, sigma, and delta_var should already be interpolated functions"


    def velInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return vel(x)
        return 0.
    
    def deltaInterp(x):
        
        if x<0.2:
            return -1
        
        if (x>0.2 or x<xpoints[-1]):
            return delta(x)
        return 0.
    
    def sigmaInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return np.max([0.2,sigma(x)])
        return sigma(xpoints[-1])
    
    
    
    def delta_deltaInterp(x):
        if x<0.2:
            return 0
        if (x>0.2 or x<xpoints[-1]):
            return delta_delta(x)
        return 0
    
    

        
    # def integrand(r3,s,muS):

    #     r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
    

    #     # print(r3)
    #     Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)

        
    #     Express=5*a*H/4*(3*muS**2-1)*(1/(sqrt2pi*sigmaInterp(r))*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))))*delta_deltaInterp(r)
    #     return Express
    
    

    muPoints=(np.linspace(0,1,numMuPoints+1)[:-1] + np.linspace(0,1,numMuPoints+1)[1:])/2 
    muIntList=np.zeros(len(muPoints))
    dMu=muPoints[1]-muPoints[0]
    
    
    def integrand(r3,s,muS):

        r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
    

        # print(r3)
        Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)

        
        #Express=5*a*H/4*(3*muS**2-1)*(r3/r)*(-1/(sqrt2pi*sigmaInterp(r)**3)*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))*Ap ))*(1+deltaInterp(r))*delta_velInterp(r)
        Express=a*H*(1/(sqrt2pi*sigmaInterp(r))*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))))*delta_deltaInterp(r)
        return Express
    
        
    if cat.scaling=="Reff":
        r3interval=1.5
        
    if cat.scaling=="R24":
        r3interval=50

    deltaXi0points=np.zeros(len(spoints))
    deltaXi2points=np.zeros(len(spoints))
    deltaXi4points=np.zeros(len(spoints))
    deltaSandMu=np.zeros((len(spoints),len(muPoints)))
                
    for i in range(len(spoints)):
        print(i)
        s=spoints[i]
        for j in range(len(muPoints)):
            muS=muPoints[j]
            
            r3points= np.linspace(s*muS-r3interval,s*muS+r3interval,200)
            # expressPoints=[integrand(r3,s,muS,delta,vel,sigma,delta_vel) for r3 in r3points]
            
            # plt.plot(r3points,expressPoints)
            # plt.show()
            # plt.clf()
            
            
            muIntList[j]=2*sp.integrate.quad(lambda r3: integrand(r3,s,muS),s*muS-r3interval,s*muS+r3interval,limit=200)[0] 
            
        deltaSandMu[i]=muIntList
            
        deltaXi0points[i]=np.sum([muIntList[i]*(2*(0) + 1)/2 * (1) for i in range(len(muIntList))])*dMu
        deltaXi2points[i]=np.sum([muIntList[i]*(2*(2) + 1)/2 * (1/2)*(3*muPoints[i]**2 - 1) for i in range(len(muIntList))])*dMu
        deltaXi4points[i]=np.sum([muIntList[i]*(2*(4) + 1)/2 * (1/8)*(35*muPoints[i]**4 - 30*muPoints[i]**2 + 3) for i in range(len(muIntList))])*dMu
            
        
        
    return spoints,deltaXi0points,deltaXi2points,deltaXi4points,deltaSandMu
"#########################################################################"    

"#########################################################################"
def calcDeltaMoments_vel(cat,cat2,spoints,rMin=35,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=101,voidTypes=["R","S"]):
    
    
    
    xpoints=cat.xpoints
    

    if cat.scaling=="Reff":
        weight=-1
    if cat.scaling=="R24":
        weight=0
    
    xpoints=cat.xpoints
    
    

    Q=0
    try:
        if cat.rMinDens==rMin and cat.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf=np.array([cat.xpoints,cat.densProf,cat.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    Q=0
    try:
        if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    

    Q=0
    try:
        if cat.rMinSigma==rMin and cat.rMaxSigma==rMax and cat.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf=cat.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
        

    Q=0
    
    try:
        if (cat2.rMinVels==rMin and cat2.rMaxVels==rMax and cat2.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat2.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf2=np.array([cat2.xpoints,cat2.velProf,cat2.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf2=cat2.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf2=cat2.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        

    delta=sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
    vel=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')
    sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf,bc_type = 'natural')
    delta_vel=sp.interpolate.CubicSpline(xpoints,velProf2[1]-velProf[1],bc_type = 'natural')
    

    
    "delta, v, sigma, and delta_var should already be interpolated functions"


    def velInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return vel(x)
        return 0.
    
    def deltaInterp(x):
        
        if x<0.2:
            return -1
        
        if (x>0.2 or x<xpoints[-1]):
            return delta(x)
        return 0.
    
    def sigmaInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return np.max([0.2,sigma(x)])
        return sigma(xpoints[-1])
    
    
    def delta_velInterp(x):
        if x<0.2:
            return 0
        if (x>0.2 or x<xpoints[-1]):
            return delta_vel(x)
        return 0
    
    
    muPoints=(np.linspace(0,1,numMuPoints+1)[:-1] + np.linspace(0,1,numMuPoints+1)[1:])/2 
    muIntList=np.zeros(len(muPoints))
    dMu=muPoints[1]-muPoints[0]
    
    
    
    
    def integrand(r3,s,muS):

        r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
    

        # print(r3)
        Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)

        
        #Express=5*a*H/4*(3*muS**2-1)*(r3/r)*(-1/(sqrt2pi*sigmaInterp(r)**3)*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))*Ap ))*(1+deltaInterp(r))*delta_velInterp(r)
        Express=a*H*(r3/r)*(-1/(sqrt2pi*sigmaInterp(r)**3)*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))*Ap ))*(1+deltaInterp(r))*delta_velInterp(r)
        return Express
    
        
    if cat.scaling=="Reff":
        r3interval=1.5
        
    if cat.scaling=="R24":
        r3interval=50

    deltaXi0points=np.zeros(len(spoints))
    deltaXi2points=np.zeros(len(spoints))
    deltaXi4points=np.zeros(len(spoints))
    deltaSandMu=np.zeros((len(spoints),len(muPoints)))
                
    for i in range(len(spoints)):
        print(i)
        s=spoints[i]
        for j in range(len(muPoints)):
            muS=muPoints[j]
            
            r3points= np.linspace(s*muS-r3interval,s*muS+r3interval,200)
            # expressPoints=[integrand(r3,s,muS,delta,vel,sigma,delta_vel) for r3 in r3points]
            
            # plt.plot(r3points,expressPoints)
            # plt.show()
            # plt.clf()
            
            
            muIntList[j]=2*sp.integrate.quad(lambda r3: integrand(r3,s,muS),s*muS-r3interval,s*muS+r3interval,limit=200)[0] 
            
        deltaSandMu[i]=muIntList
            
        deltaXi0points[i]=np.sum([muIntList[i]*(2*(0) + 1)/2 * (1) for i in range(len(muIntList))])*dMu
        deltaXi2points[i]=np.sum([muIntList[i]*(2*(2) + 1)/2 * (1/2)*(3*muPoints[i]**2 - 1) for i in range(len(muIntList))])*dMu
        deltaXi4points[i]=np.sum([muIntList[i]*(2*(4) + 1)/2 * (1/8)*(35*muPoints[i]**4 - 30*muPoints[i]**2 + 3) for i in range(len(muIntList))])*dMu
            
        
        
    return spoints,deltaXi0points,deltaXi2points,deltaXi4points,deltaSandMu
"#########################################################################"    

"#########################################################################"
def calcDeltaMoments_sigma(cat,cat2,spoints,rMin=35,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=101,voidTypes=["R","S"]):
 
        
    xpoints=cat.xpoints
    

    if cat.scaling=="Reff":
        weight=-1
    if cat.scaling=="R24":
        weight=0
    
    xpoints=cat.xpoints
    
    

    Q=0
    try:
        if cat.rMinDens==rMin and cat.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf=np.array([cat.xpoints,cat.densProf,cat.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    Q=0
    try:
        if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    

    Q=0
    try:
        if cat.rMinSigma==rMin and cat.rMaxSigma==rMax and cat.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf=cat.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
        

    Q=0

    try:
        if cat2.rMinSigma==rMin and cat2.rMaxSigma==rMax and cat2.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf2=cat2.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf2 = cat2.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf2 = cat2.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
            

    delta=sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
    vel=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')
    sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf,bc_type = 'natural')
    delta_sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf2-sigmaProf,bc_type = 'natural')
    
 
    
    "delta, v, sigma, and delta_var should already be interpolated functions"
    xpoints=cat.xpoints
    

    def velInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return vel(x)
        return 0.
    
    def deltaInterp(x):
        
        if x<0.2:
            return -1
        
        if (x>0.2 or x<xpoints[-1]):
            return delta(x)
        return 0.
    
    def sigmaInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return np.max([0.2,sigma(x)])
        return sigma(xpoints[-1])
    
    
    def delta_sigmaInterp(x):
        if x<0.2:
            return 0
        if (x>0.2 or x<xpoints[-1]):
            return delta_sigma(x)
        return delta_sigma(xpoints[-1])
    
    
    muPoints=(np.linspace(0,1,numMuPoints+1)[:-1] + np.linspace(0,1,numMuPoints+1)[1:])/2 
    muIntList=np.zeros(len(muPoints))
    dMu=muPoints[1]-muPoints[0]
    
    
    
    
    def integrand(r3,s,muS):

        r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
    

        # print(r3)
        Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)

        
        #Express=5*a*H/4*(3*muS**2-1)*((-1/(sqrt2pi*sigmaInterp(r)**2)*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2)) )) + (1/(sqrt2pi*sigmaInterp(r)))*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))*(Ap**2/(sigmaInterp(r)**3))))*(1+deltaInterp(r))*delta_sigmaInterp(r)
        Express=a*H*((-1/(sqrt2pi*sigmaInterp(r)**2)*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2)) )) + (1/(sqrt2pi*sigmaInterp(r)))*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))*(Ap**2/(sigmaInterp(r)**3))))*(1+deltaInterp(r))*delta_sigmaInterp(r)
        
        
        return Express
    
        
    if cat.scaling=="Reff":
        r3interval=0.75
        
    if cat.scaling=="R24":
        r3interval=50

    
    # deltaXi2points=np.zeros(len(spoints))
    # for i in range(len(spoints)):
    #     print(i)
    #     s=spoints[i]
    #     for j in range(len(muPoints)):
    #         muS=muPoints[j]
            
    #         # r3points= np.linspace(s*muS-r3interval,s*muS+r3interval,200)
    #         # expressPoints=[integrand(r3,s,muS,delta,vel,sigma,delta_sigma) for r3 in r3points]
            
    #         # plt.plot(r3points,expressPoints)
    #         # plt.show()
    #         # plt.clf()
            
            
    #         muIntList[j]=2*sp.integrate.quad(lambda r3: integrand(r3,s,muS),s*muS-r3interval,s*muS+r3interval,limit=200)[0] 
            
    #     deltaXi2points[i]=np.sum(muIntList)*dMu
        
            
        
        
    # return spoints,deltaXi2points
    
    deltaXi0points=np.zeros(len(spoints))
    deltaXi2points=np.zeros(len(spoints))
    deltaXi4points=np.zeros(len(spoints))
    
    deltaSandMu=np.zeros((len(spoints),len(muPoints)))

    
    
    
    
    for i in range(len(spoints)):
        print(i)
        s=spoints[i]
        for j in range(len(muPoints)):
            muS=muPoints[j]
            
            r3points= np.linspace(s*muS-r3interval,s*muS+r3interval,200)
            # expressPoints=[integrand(r3,s,muS,delta,vel,sigma,delta_vel) for r3 in r3points]
            
            # plt.plot(r3points,expressPoints)
            # plt.show()
            # plt.clf()
            
            muIntList[j]=2*sp.integrate.quad(lambda r3: integrand(r3,s,muS),s*muS-r3interval,s*muS+r3interval,limit=200)[0] 
            
        deltaSandMu[i]=muIntList
        deltaXi0points[i]=np.sum([muIntList[i]*(2*(0) + 1)/2 * (1) for i in range(len(muIntList))])*dMu
        deltaXi2points[i]=np.sum([muIntList[i]*(2*(2) + 1)/2 * (1/2)*(3*muPoints[i]**2 - 1) for i in range(len(muIntList))])*dMu
        deltaXi4points[i]=np.sum([muIntList[i]*(2*(4) + 1)/2 * (1/8)*(35*muPoints[i]**4 - 30*muPoints[i]**2 + 3) for i in range(len(muIntList))])*dMu
            
        
        
    return spoints,deltaXi0points,deltaXi2points,deltaXi4points,deltaSandMu
"#########################################################################"    

"###############################################################################################################"

"#########################################################################"
def calcDeltaXi2_omegaHoldSigma(cat,cat2,spoints,rMin=25,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=101,voidTypes=["R","S"]):
 
        
    xpoints=cat.xpoints
    

    if cat.scaling=="Reff":
        weight=-1
    if cat.scaling=="R24":
        weight=0
    
    xpoints=cat.xpoints
    
    

    Q=0
    try:
        if cat.rMinDens==rMin and cat.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf=np.array([cat.xpoints,cat.densProf,cat.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    Q=0
    try:
        if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    

    Q=0
    try:
        if cat.rMinSigma==rMin and cat.rMaxSigma==rMax and cat.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf=cat.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
        

    Q=0

    try:
        if cat2.rMinSigma==rMin and cat2.rMaxSigma==rMax and cat2.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf2=cat2.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf2 = cat2.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf2 = cat2.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
            
        
    Q=0
    try:
        if (cat2.rMinVels==rMin and cat2.rMaxVels==rMax and cat2.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf2=np.array([cat2.xpoints,cat2.velProf,cat2.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf2=cat2.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf2=cat2.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    
    
        
        
        
        
        
        

    delta=sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
    vel=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')
    sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf,bc_type = 'natural')
    sigma2=sp.interpolate.CubicSpline(xpoints,sigmaProf2,bc_type = 'natural')
    vel2=sp.interpolate.CubicSpline(xpoints,velProf2[1],bc_type = 'natural')
    
 
    
    "delta, v, sigma, and delta_var should already be interpolated functions"
    xpoints=cat.xpoints
    

    def velInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return vel(x)
        return 0.
    
    def deltaInterp(x):
        
        if x<0.2:
            return -1
        
        if (x>0.2 or x<xpoints[-1]):
            return delta(x)
        return 0.
    
    def sigmaInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return np.max([0.2,sigma(x)])
        return sigma(xpoints[-1])
    
    
    def omegaInterp(x):
        return velInterp(x)/sigmaInterp(x)
    
    
    
    def sigmaInterp2(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return np.max([0.2,sigma2(x)])
        return sigma2(xpoints[-1])
    
    
    def velInterp2(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return vel2(x)
        return 0.
    

    def delta_omegaInterp(x):
        if x<0.2:
            return 0
        if (x>0.2 or x<xpoints[-1]):
            return velInterp2(x)/sigmaInterp2(x) - velInterp(x)/sigmaInterp(x)
        return velInterp2(xpoints[-1])/sigmaInterp2(xpoints[-1]) - velInterp(xpoints[-1])/sigmaInterp(xpoints[-1])
    
    
    
    
    
    
    
    muPoints=(np.linspace(0,1,numMuPoints+1)[:-1] + np.linspace(0,1,numMuPoints+1)[1:])/2 
    muIntList=np.zeros(len(muPoints))
    dMu=muPoints[1]-muPoints[0]
    
    
    
    
    def integrand(r3,s,muS):

        r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
    

        Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)
        
        # Express=5*a*H/4*(3*muS**2-1)*((-1/(sqrt2pi*sigmaInterp(r)**2)*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2)) )) + (1/(sqrt2pi*sigmaInterp(r)))*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))*(Ap**2/(sigmaInterp(r)**3))))*(1+deltaInterp(r))*delta_sigmaInterp(r)
        
        Express=5*a*H/4*(3*muS**2-1)*(-((Ap*r3/r)*np.exp(-Ap**2/(2.*sigmaInterp(r)**2))/(sqrt2pi*sigmaInterp(r)**2)))*(1 + deltaInterp(r))*delta_omegaInterp(r)
        
        return Express
    
        
    if cat.scaling=="Reff":
        r3interval=0.75
        
    if cat.scaling=="R24":
        r3interval=50

    
    deltaXi2points=np.zeros(len(spoints))
    for i in range(len(spoints)):
        print(i)
        s=spoints[i]
        for j in range(len(muPoints)):
            muS=muPoints[j]
            
            # r3points= np.linspace(s*muS-r3interval,s*muS+r3interval,200)
            # expressPoints=[integrand(r3,s,muS,delta,vel,sigma,delta_sigma) for r3 in r3points]
            
            # plt.plot(r3points,expressPoints)
            # plt.show()
            # plt.clf()
            
            
            muIntList[j]=2*sp.integrate.quad(lambda r3: integrand(r3,s,muS),s*muS-r3interval,s*muS+r3interval,limit=200)[0] 
            
        deltaXi2points[i]=np.sum(muIntList)*dMu
        
            
        
        
    return spoints,deltaXi2points
"#########################################################################"    

"#########################################################################"
def calcDeltaXi2_sigmaHoldOmega(cat,cat2,spoints,rMin=25,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=101,voidTypes=["R","S"]):
 
        
    xpoints=cat.xpoints
    

    if cat.scaling=="Reff":
        weight=-1
    if cat.scaling=="R24":
        weight=0
    
    xpoints=cat.xpoints
    
    

    Q=0
    try:
        if cat.rMinDens==rMin and cat.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf=np.array([cat.xpoints,cat.densProf,cat.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    Q=0
    try:
        if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    

    Q=0
    try:
        if cat.rMinSigma==rMin and cat.rMaxSigma==rMax and cat.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf=cat.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
        

    Q=0

    try:
        if cat2.rMinSigma==rMin and cat2.rMaxSigma==rMax and cat2.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf2=cat2.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf2 = cat2.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf2 = cat2.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
            
        
    Q=0
    try:
        if (cat2.rMinVels==rMin and cat2.rMaxVels==rMax and cat2.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf2=np.array([cat2.xpoints,cat2.velProf,cat2.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf2=cat2.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf2=cat2.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    
    
        
        
        
        
        
        

    delta=sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
    vel=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')
    sigma=sp.interpolate.CubicSpline(xpoints,sigmaProf,bc_type = 'natural')
    sigma2=sp.interpolate.CubicSpline(xpoints,sigmaProf2,bc_type = 'natural')
    vel2=sp.interpolate.CubicSpline(xpoints,velProf2[1],bc_type = 'natural')
    
 
    
    "delta, v, sigma, and delta_var should already be interpolated functions"
    xpoints=cat.xpoints
    

    def velInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return vel(x)
        return 0.
    
    def deltaInterp(x):
        
        if x<0.2:
            return -1
        
        if (x>0.2 or x<xpoints[-1]):
            return delta(x)
        return 0.
    
    def sigmaInterp(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return np.max([0.2,sigma(x)])
        return sigma(xpoints[-1])
    
    
    def omegaInterp(x):
        return velInterp(x)/sigmaInterp(x)
    
    
    
    def sigmaInterp2(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return np.max([0.2,sigma2(x)])
        return sigma2(xpoints[-1])
    
    
    def velInterp2(x):
        # print(x)
        if (x>0 or x<xpoints[-1]):
            return vel2(x)
        return 0.
    

    def delta_omegaInterp(x):
        if x<0.2:
            return 0
        if (x>0.2 or x<xpoints[-1]):
            return velInterp2(x)/sigmaInterp2(x) - velInterp(x)/sigmaInterp(x)
        return velInterp2(xpoints[-1])/sigmaInterp2(xpoints[-1]) - velInterp(xpoints[-1])/sigmaInterp(xpoints[-1])
    
    
    
    
    
    
    
    muPoints=(np.linspace(0,1,numMuPoints+1)[:-1] + np.linspace(0,1,numMuPoints+1)[1:])/2 
    muIntList=np.zeros(len(muPoints))
    dMu=muPoints[1]-muPoints[0]
    
    
    
    
    def integrand(r3,s,muS):

        r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
    

        Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)
        
        
        Express=5*a*H/4*(3*muS**2-1)*(np.exp(-Ap**2/(2.*sigmaInterp(r)**2)))*((-(1/(sqrt2pi*sigmaInterp(r)**2)) + (Ap**2/sigmaInterp(r)**3 - (Ap*r3*omegaInterp(r))/(r*sigmaInterp(r)**2))/(sqrt2pi*sigmaInterp(r))))*(1 + deltaInterp(r))*(sigmaInterp2(r) - sigmaInterp(r))
        return Express
    
        
    if cat.scaling=="Reff":
        r3interval=0.75
        
    if cat.scaling=="R24":
        r3interval=50

    
    deltaXi2points=np.zeros(len(spoints))
    for i in range(len(spoints)):
        print(i)
        s=spoints[i]
        for j in range(len(muPoints)):
            muS=muPoints[j]
            
            # r3points= np.linspace(s*muS-r3interval,s*muS+r3interval,200)
            # expressPoints=[integrand(r3,s,muS,delta,vel,sigma,delta_sigma) for r3 in r3points]
            
            # plt.plot(r3points,expressPoints)
            # plt.show()
            # plt.clf()
            
            
            muIntList[j]=2*sp.integrate.quad(lambda r3: integrand(r3,s,muS),s*muS-r3interval,s*muS+r3interval,limit=200)[0] 
            
        deltaXi2points[i]=np.sum(muIntList)*dMu
        
            
        
        
    return spoints,deltaXi2points
"#########################################################################" 

"#####"


"###############################################################################################################"

"#########################################################################"
def plotArgDeltaXi2_sigma(cat,cat2,s,muS,rMin=25,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"]):

    xpoints=cat.xpoints
    

    if cat.scaling=="Reff":
        weight=-1
    if cat.scaling=="R24":
        weight=0
    
    xpoints=cat.xpoints
    
    

    Q=0
    try:
        if cat.rMinDens==rMin and cat.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf=np.array([cat.xpoints,cat.densProf,cat.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    Q=0
    try:
        if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    

    Q=0
    try:
        if cat.rMinSigma==rMin and cat.rMaxSigma==rMax and cat.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf=cat.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
        
    Q=0
    try:
        if cat2.rMinSigma==rMin and cat2.rMaxSigma==rMax and cat2.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf2=cat2.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf2 = cat2.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf2 = cat2.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
            
        
        
        
        
        
    # densityProf=cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax)
    # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,weight=weight)
    # sigmaProf=cat.calcSigmaV_r(rMin=rMin,rMax=rMax,weight=weight)
    
    # for i in range(len(sigmaProf)):
    #     if sigmaProf[i]==0:
    #         sigmaProf[i]=0.2
    
    deltaRInterp=sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
    velRInterp=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')
    sigmaRInterp=sp.interpolate.CubicSpline(xpoints,sigmaProf,bc_type = 'natural')
    deltaSigmaRInterp=sp.interpolate.CubicSpline(xpoints,sigmaProf2-sigmaProf,bc_type = 'natural')
    
    def deltaInterp(x):
        if (x>0 and x<xpoints[-1]):
            return deltaRInterp(x)
        if (x<0):
            return -1.
        return 0.
    
    
    def velInterp(x):
    
        if (x>0 or x<xpoints[-1]):
            return velRInterp(x)
        return 0.
    
    
    
    def sigmaInterp(r):
        if (r>0 and r<xpoints[-1]):
            return sigmaRInterp(r)
        if r<0:
            return 0
        return sigmaProf[-1]
    
    
    
    def deltaSigmaInterp(r):
        if (r>0 and r<xpoints[-1]):
            return deltaSigmaRInterp(r)
        if r<0:
            return 0
        return sigmaProf[-1]-sigmaProf2[-1]
    
    
    
    
    # muPoints=(np.linspace(0,1,numMuPoints+1)[:-1] + np.linspace(0,1,numMuPoints+1)[1:])/2 
    # muIntList=np.zeros(len(muPoints))
    # dMu=muPoints[1]-muPoints[0]
    
    
    
    
    def integrand(r3,s,muS,delta,vel,sigma,delta_sigmaInterp):

        r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) #might have to put some "try" statements here if we have any numerical BS
    

        # print(r3)
        Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)

        
        term1=a*H*((-1/(sqrt2pi*sigmaInterp(r)**2)*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2)) )))#*(1+delta(r))
        term2=a*H*(1/(sqrt2pi*sigmaInterp(r)))*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2))*(Ap**2/(sigmaInterp(r)**3)))#*(1+delta(r))
        
        

        return term1,term2,delta_sigmaInterp(r)
    
    
    coords0=coordinateChangeStoR(cat,s=s,muS=muS,rMin=rMin,rMax=rMax)
    
    r0=coords0[0]
    muR0=coords0[1]


    def integrandConst(r3,s,muS,delta,vel,sigma,delta_sigmaInterp):
        
        Ap=(r3-muS*s)*a*H+muR0*velInterp(r0)
    
        term1=a*H*((-1/(sqrt2pi*sigmaInterp(r0)**2)*(np.exp((-Ap**2)/(2*sigmaInterp(r0)**2)) )))#*(1+delta(r0))
        term2=a*H*(1/(sqrt2pi*sigmaInterp(r0)))*(np.exp((-Ap**2)/(2*sigmaInterp(r0)**2))*(Ap**2/(sigmaInterp(r0)**3)))#*(1+delta(r0))
    
    
        return term1 + term2
    
    
        
    if cat.scaling=="Reff":
        r3interval=0.75
        
    if cat.scaling=="R24":
        r3interval=50

    
    sPerp=np.sqrt((s**2)*(1.-muS**2))




    r3points= np.linspace(s*muS-r3interval,s*muS+r3interval,400)
    expressPoints=[integrand(r3,s,muS,delta=deltaInterp,vel=velInterp,sigma=sigmaInterp,delta_sigmaInterp=deltaSigmaInterp) for r3 in r3points]
    term1Points=np.array([expressPoints[i][0] for i in range(len(expressPoints))])
    term2Points=np.array([expressPoints[i][1] for i in range(len(expressPoints))])
    deltaSigmaPoints=np.array([expressPoints[i][2] for i in range(len(expressPoints))])
    
    densProfPoints=np.array([deltaInterp(np.sqrt(r3**2 + sPerp**2))+1 for r3 in r3points])
    
    
    integrandPoints=term1Points + term2Points
    integrandConstPoints=np.array([integrandConst(r3,s,muS,delta=deltaInterp,vel=velInterp,sigma=sigmaInterp,delta_sigmaInterp=deltaSigmaInterp) for r3 in r3points])
    
    
    dr3=r3points[1]-r3points[0]
    val=(integrandPoints@densProfPoints)*dr3
    
    
    fullVal=(integrandPoints*densProfPoints)@deltaSigmaPoints*dr3
    
    print("Sigma Integral at s = " + str(s) + ", muS = " + str(muS) + " is " + str(val))
    print("Sigma fullVal  at s = " + str(s) + ", muS = " + str(muS) + " is " + str(fullVal))

        
    return [r3points,integrandPoints,densProfPoints,deltaSigmaPoints,integrandConstPoints],val,fullVal
"#########################################################################"    

"#########################################################################"
def plotArgDeltaXi2_vel(cat,cat2,s,muS,rMin=25,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"]):

    xpoints=cat.xpoints
    

    if cat.scaling=="Reff":
        weight=-1
    if cat.scaling=="R24":
        weight=0
    
    xpoints=cat.xpoints
    
    

    Q=0
    try:
        if cat.rMinDens==rMin and cat.rMaxDens==rMax:
            # print("loading density Profile")
            densityProf=np.array([cat.xpoints,cat.densProf,cat.densErrors])
            
            Q=1
        if Q==0:
            print("densityProf doesnt match required params, calculating now")
            densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
            
    except:
        print("densityProf not pre-calculated, calculating now")
        densityProf = cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes)
                
    

    Q=0
    try:
        if (cat.rMinVels==rMin and cat.rMaxVels==rMax and cat.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf=np.array([cat.xpoints,cat.velProf,cat.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        
    

    Q=0
    try:
        if cat.rMinSigma==rMin and cat.rMaxSigma==rMax and cat.sigmaWeight==weight:
            # print("loading sigma profile")
            sigmaProf=cat.sigmaV_r
            
            Q=1
        if Q==0:
            print("Sigma profile doesnt match required params, calculating now")
            sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
            
    except:
        print("sigma Prof not pre-calculated, calculating now")
        sigmaProf = cat.calcSigmaV_r(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)   
        
        

    Q=0
    
    try:
        if (cat2.rMinVels==rMin and cat2.rMaxVels==rMax and cat2.velWeight==weight):
            # print("loading Velocity Profile")
            # velProf=cat2.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            velProf2=np.array([cat2.xpoints,cat2.velProf,cat2.velErrors])
            
            
            Q=1
        if Q==0:
            print("velocity doesnt match required params, calculating now")
            velProf2=cat2.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
            
    except:
            print("velocity not pre-calculated, calculating now")
            velProf2=cat2.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,voidTypes=voidTypes,weight=weight)
        

        
        
        
    # densityProf=cat.jackKnifeDensityProfile(rMin=rMin,rMax=rMax)
    # velProf=cat.jackKnifeVelocityProfile(rMin=rMin,rMax=rMax,weight=weight)
    # sigmaProf=cat.calcSigmaV_r(rMin=rMin,rMax=rMax,weight=weight)
    
    # for i in range(len(sigmaProf)):
    #     if sigmaProf[i]==0:
    #         sigmaProf[i]=0.2
    
    deltaRInterp=sp.interpolate.CubicSpline(xpoints,densityProf[1],bc_type = 'natural')
    velRInterp=sp.interpolate.CubicSpline(xpoints,velProf[1],bc_type = 'natural')
    sigmaRInterp=sp.interpolate.CubicSpline(xpoints,sigmaProf,bc_type = 'natural')
    deltaVelRInterp=sp.interpolate.CubicSpline(xpoints,velProf2[1]-velProf[1],bc_type = 'natural')
    
    def deltaInterp(x):
        if (x>0 and x<xpoints[-1]):
            return deltaRInterp(x)
        if (x<0):
            return -1.
        return 0.
    
    
    def velInterp(x):
    
        if (x>0 or x<xpoints[-1]):
            return velRInterp(x)
        return 0.
    
    
    
    def sigmaInterp(r):
        if (r>0 and r<xpoints[-1]):
            return sigmaRInterp(r)
        if r<0:
            return 0
        return sigmaProf[-1]
    
    
    
    def deltaVelInterp(r):
        if (r>0 and r<xpoints[-1]):
            return deltaVelRInterp(r)
        if r<0:
            return 0
        return velProf2[1][-1]-velProf[1][-1]
    
    
    
    
    # muPoints=(np.linspace(0,1,numMuPoints+1)[:-1] + np.linspace(0,1,numMuPoints+1)[1:])/2 
    # muIntList=np.zeros(len(muPoints))
    # dMu=muPoints[1]-muPoints[0]
    
    
    
    
    def integrand(r3,s,muS,delta,vel,sigma,delta_velInterp):

        r=np.sqrt(r3**2+(s**2)*(1.-muS**2)) 

        Ap=(r3-muS*s)*a*H+r3/r*velInterp(r)
        
        term1=a*H*(r3/r*(-1/(sqrt2pi*sigmaInterp(r)**3)*(np.exp((-Ap**2)/(2*sigmaInterp(r)**2)) )*Ap))#*(1+delta(r))

        return term1,delta_velInterp(r)
    
    
    coords0=coordinateChangeStoR(cat,s=s,muS=muS,rMin=rMin,rMax=rMax)
    
    r0=coords0[0]
    muR0=coords0[1]



    def integrandConst(r3,s,muS,delta,vel,sigma,delta_velInterp):
        "First, solve for the r3 that makes the expression in the exponential equal to 0"

        Ap=(r3-muS*s)*a*H+muR0*velInterp(r0)
        
        "Might want to double check on the muR0 factor out in front here"
        term1=a*H*(muR0*(-1/(sqrt2pi*sigmaInterp(r0)**3)*(np.exp((-Ap**2)/(2*sigmaInterp(r0)**2)) )*Ap))#*(1+delta(r))

        return term1
    

  
    if cat.scaling=="Reff":
        r3interval=0.75
        
    if cat.scaling=="R24":
        r3interval=50

    sPerp=np.sqrt((s**2)*(1.-muS**2))

    r3points = np.linspace(s*muS-r3interval,s*muS+r3interval,400)
    expressPoints=[integrand(r3,s,muS,delta=deltaInterp,vel=velInterp,sigma=sigmaInterp,delta_velInterp=deltaVelInterp) for r3 in r3points]
    term1Points=np.array([expressPoints[i][0] for i in range(len(expressPoints))])
    deltaVelPoints=np.array([expressPoints[i][1] for i in range(len(expressPoints))])
    
    
    densProfPoints=np.array([deltaInterp(np.sqrt(r3**2 + sPerp**2))+1 for r3 in r3points])
    
    
    integrandPoints=term1Points
    integrandConstPoints=np.array([integrandConst(r3,s,muS,delta=deltaInterp,vel=velInterp,sigma=sigmaInterp,delta_velInterp=deltaVelInterp) for r3 in r3points])
    
    
    
    dr3=r3points[1]-r3points[0]
    val=(integrandPoints@densProfPoints)*dr3
    
    fullVal=(integrandPoints*densProfPoints)@deltaVelPoints*dr3
    print("Vel Integral at s = " + str(s) + ", muS = " + str(muS) + " is " + str(val))
    print("Vel FullVal  at s = " + str(s) + ", muS = " + str(muS) + " is " + str(fullVal))
    

        
    return [r3points,integrandPoints,densProfPoints,deltaVelPoints,integrandConstPoints],val,fullVal
"#########################################################################"    






"Load Halo Voids for LOS=Z or others"
numsToLoad=range(1,101)
if True:
    
    "############# Void Info Files #############"
    
    'voidInfo_CC_R24_z0p5'
    GR_voidInfo_CC_R24_z0p5=[np.load(voidDir + "GR" + str(num) + "_voidInfo_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    F5_voidInfo_CC_R24_z0p5=[np.load(voidDir + "F5" + str(num) + "_voidInfo_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    N1_voidInfo_CC_R24_z0p5=[np.load(voidDir + "N1" + str(num) + "_voidInfo_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]


    "############# Void Density and Velocity Profiles #############"

    'densVelsSigmas_CC_R24_z0p5'
    GR_densVelsSigmas_CC_R24_z0p5=[np.load(voidDir + "GR" + str(num) + "_densVelsSigmas_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    F5_densVelsSigmas_CC_R24_z0p5=[np.load(voidDir + "F5" + str(num) + "_densVelsSigmas_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    N1_densVelsSigmas_CC_R24_z0p5=[np.load(voidDir + "N1" + str(num) + "_densVelsSigmas_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]


    "############# Assemble Z catalogs here #############"    
    catGR_CC_R24_z0p5=catalogClass(voidInfos=GR_voidInfo_CC_R24_z0p5,densVelsSigmas=GR_densVelsSigmas_CC_R24_z0p5,scaling="R24",centerType="CC",gravType="GR")
    catF5_CC_R24_z0p5=catalogClass(voidInfos=F5_voidInfo_CC_R24_z0p5,densVelsSigmas=F5_densVelsSigmas_CC_R24_z0p5,scaling="R24",centerType="CC",gravType="F5")
    catN1_CC_R24_z0p5=catalogClass(voidInfos=N1_voidInfo_CC_R24_z0p5,densVelsSigmas=N1_densVelsSigmas_CC_R24_z0p5,scaling="R24",centerType="CC",gravType="N1")




"Load 'HOD' Voids for LOS=Z"
numsToLoad=range(1,101)
if True: 
    voidDirHOD="/Users/christopherwilson/Desktop/Research/GLAM/125/catalogs/"
    "############# Void Info Files #############"
    
    'voidInfo_CC_R24_z0p5'
    GRHOD_voidInfo_CC_R24_z0p5=[np.load(voidDirHOD + "GRHOD" + str(num) + "_voidInfo_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    F5HOD_voidInfo_CC_R24_z0p5=[np.load(voidDirHOD + "F5HOD" + str(num) + "_voidInfo_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    N1HOD_voidInfo_CC_R24_z0p5=[np.load(voidDirHOD + "N1HOD" + str(num) + "_voidInfo_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    
 
    "############# Void Density and Velocity Profiles #############"


    'densVelsSigmas_CC_R24_z0p5'
    GRHOD_densVelsSigmas_CC_R24_z0p5=[np.load(voidDirHOD + "GRHOD" + str(num) + "_densVelsSigmas_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    F5HOD_densVelsSigmas_CC_R24_z0p5=[np.load(voidDirHOD + "F5HOD" + str(num) + "_densVelsSigmas_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    N1HOD_densVelsSigmas_CC_R24_z0p5=[np.load(voidDirHOD + "N1HOD" + str(num) + "_densVelsSigmas_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]




    "############# Void muZProfZ files #############"
    'muZProfZ_CC_R24_z0p5'
    GRHOD_muZProfZ_CC_R24_z0p5=[np.load(voidDirHOD + "GRHOD" + str(num) + "_muZProfZ_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    F5HOD_muZProfZ_CC_R24_z0p5=[np.load(voidDirHOD + "F5HOD" + str(num) + "_muZProfZ_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    N1HOD_muZProfZ_CC_R24_z0p5=[np.load(voidDirHOD + "N1HOD" + str(num) + "_muZProfZ_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]


    "############# Assemble catalogs here #############"
    # start=time.time()
    catGRHOD_CC_R24_z0p5=catalogClass(voidInfos=GRHOD_voidInfo_CC_R24_z0p5,partsRandMu=GRHOD_muZProfZ_CC_R24_z0p5,densVelsSigmas=GRHOD_densVelsSigmas_CC_R24_z0p5,scaling="R24",centerType="CC",gravType="GRHOD")
    catF5HOD_CC_R24_z0p5=catalogClass(voidInfos=F5HOD_voidInfo_CC_R24_z0p5,partsRandMu=F5HOD_muZProfZ_CC_R24_z0p5,densVelsSigmas=F5HOD_densVelsSigmas_CC_R24_z0p5,scaling="R24",centerType="CC",gravType="F5HOD")
    catN1HOD_CC_R24_z0p5=catalogClass(voidInfos=N1HOD_voidInfo_CC_R24_z0p5,partsRandMu=N1HOD_muZProfZ_CC_R24_z0p5,densVelsSigmas=N1HOD_densVelsSigmas_CC_R24_z0p5,scaling="R24",centerType="CC",gravType="N1HOD")
    
    

"Load 'HOD' Voids, combining lines of sight, LOS=XYZ"
numsToLoad=range(1,101)
if True: 
    voidDirHOD="/Users/christopherwilson/Desktop/Research/GLAM/125/catalogs/"
    "############# Void Info Files #############"
    

    'voidInfo_CC_R24_z0p5'
    GRHOD_voidInfo_CC_R24_z0p5=[np.load(voidDirHOD + "GRHOD" + str(num) + "_voidInfo_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    F5HOD_voidInfo_CC_R24_z0p5=[np.load(voidDirHOD + "F5HOD" + str(num) + "_voidInfo_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    N1HOD_voidInfo_CC_R24_z0p5=[np.load(voidDirHOD + "N1HOD" + str(num) + "_voidInfo_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]


    "############# Void Density and Velocity Profiles #############"
    'densVelsSigmas_CC_R24_z0p5'
    GRHOD_densVelsSigmas_CC_R24_z0p5=[np.load(voidDirHOD + "GRHOD" + str(num) + "_densVelsSigmas_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    F5HOD_densVelsSigmas_CC_R24_z0p5=[np.load(voidDirHOD + "F5HOD" + str(num) + "_densVelsSigmas_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    N1HOD_densVelsSigmas_CC_R24_z0p5=[np.load(voidDirHOD + "N1HOD" + str(num) + "_densVelsSigmas_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]

    "############# Void muXProfZ files #############"
    'muXProfZ_CC_R24_z0p5'
    GRHOD_muXProfZ_CC_R24_z0p5=[np.load(voidDirHOD + "GRHOD" + str(num) + "_muXProfZ_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    F5HOD_muXProfZ_CC_R24_z0p5=[np.load(voidDirHOD + "F5HOD" + str(num) + "_muXProfZ_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    N1HOD_muXProfZ_CC_R24_z0p5=[np.load(voidDirHOD + "N1HOD" + str(num) + "_muXProfZ_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]

    "############# Void muYProfZ files #############"
    'muYProfZ_CC_R24_z0p5'
    GRHOD_muYProfZ_CC_R24_z0p5=[np.load(voidDirHOD + "GRHOD" + str(num) + "_muYProfZ_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    F5HOD_muYProfZ_CC_R24_z0p5=[np.load(voidDirHOD + "F5HOD" + str(num) + "_muYProfZ_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    N1HOD_muYProfZ_CC_R24_z0p5=[np.load(voidDirHOD + "N1HOD" + str(num) + "_muYProfZ_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]

    "############# Void muZProfZ files #############"
    'muZProfZ_CC_R24_z0p5'
    GRHOD_muZProfZ_CC_R24_z0p5=[np.load(voidDirHOD + "GRHOD" + str(num) + "_muZProfZ_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    F5HOD_muZProfZ_CC_R24_z0p5=[np.load(voidDirHOD + "F5HOD" + str(num) + "_muZProfZ_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]
    N1HOD_muZProfZ_CC_R24_z0p5=[np.load(voidDirHOD + "N1HOD" + str(num) + "_muZProfZ_CC_R24_z0p5.npy",allow_pickle=True,encoding="latin1") for num in numsToLoad]



    "############# Void muXYZProfZ creation ############# - Option 1 - each void gets the average of its X + Y + Z"
    
    # 'muZProfZ_CC_R24_z0p5'
    # GRHOD_muXYZProfZ_CC_R24_z0p5=[[np.zeros((50,50)) for j in range(len(GRHOD_voidInfo_CC_R24_z0p5[i]))] for i in range(len(numsToLoad))]
    # for i in range(len(GRHOD_muXYZProfZ_CC_R24_z0p5)):
    #     for j in range(len(GRHOD_muXYZProfZ_CC_R24_z0p5[i])):
    #         GRHOD_muXYZProfZ_CC_R24_z0p5[i][j] = np.array(GRHOD_muXProfZ_CC_R24_z0p5[i][j]) + np.array(GRHOD_muYProfZ_CC_R24_z0p5[i][j]) + np.array(GRHOD_muZProfZ_CC_R24_z0p5[i][j])

    # F5HOD_muXYZProfZ_CC_R24_z0p5=[[np.zeros((50,50)) for j in range(len(F5HOD_voidInfo_CC_R24_z0p5[i]))] for i in range(len(numsToLoad))]
    # for i in range(len(F5HOD_muXYZProfZ_CC_R24_z0p5)):
    #     for j in range(len(F5HOD_muXYZProfZ_CC_R24_z0p5[i])):
    #         F5HOD_muXYZProfZ_CC_R24_z0p5[i][j] = np.array(F5HOD_muXProfZ_CC_R24_z0p5[i][j]) + np.array(F5HOD_muYProfZ_CC_R24_z0p5[i][j]) + np.array(F5HOD_muZProfZ_CC_R24_z0p5[i][j])

    # N1HOD_muXYZProfZ_CC_R24_z0p5=[[np.zeros((50,50)) for j in range(len(N1HOD_voidInfo_CC_R24_z0p5[i]))] for i in range(len(numsToLoad))]
    # for i in range(len(N1HOD_muXYZProfZ_CC_R24_z0p5)):
    #     for j in range(len(N1HOD_muXYZProfZ_CC_R24_z0p5[i])):
    #         N1HOD_muXYZProfZ_CC_R24_z0p5[i][j] = np.array(N1HOD_muXProfZ_CC_R24_z0p5[i][j]) + np.array(N1HOD_muYProfZ_CC_R24_z0p5[i][j]) + np.array(N1HOD_muZProfZ_CC_R24_z0p5[i][j])



    
    "############# Void muXYZProfZ creation ############# - Option 2 - X, Y, and Z treated like seperate voids"    
    'muZProfZ_CC_R24_z0p5'
    
    
    GRHOD_muXYZProfZ_CC_R24_z0p5=[[] for i in range(300)]
    GRHOD_muXYZProfZ_CC_R24_z0p5[:100]=GRHOD_muXProfZ_CC_R24_z0p5
    GRHOD_muXYZProfZ_CC_R24_z0p5[100:200]=GRHOD_muYProfZ_CC_R24_z0p5
    GRHOD_muXYZProfZ_CC_R24_z0p5[200:]=GRHOD_muZProfZ_CC_R24_z0p5
    
    GRHOD_voidInfoXYZ_CC_R24_z0p5=[[] for i in range(300)]
    GRHOD_voidInfoXYZ_CC_R24_z0p5[:100]=GRHOD_voidInfo_CC_R24_z0p5
    GRHOD_voidInfoXYZ_CC_R24_z0p5[100:200]=GRHOD_voidInfo_CC_R24_z0p5
    GRHOD_voidInfoXYZ_CC_R24_z0p5[200:]=GRHOD_voidInfo_CC_R24_z0p5

    GRHOD_densVelsSigmasXYZ_CC_R24_z0p5=[[] for i in range(300)]
    GRHOD_densVelsSigmasXYZ_CC_R24_z0p5[:100]=GRHOD_densVelsSigmas_CC_R24_z0p5
    GRHOD_densVelsSigmasXYZ_CC_R24_z0p5[100:200]=GRHOD_densVelsSigmas_CC_R24_z0p5
    GRHOD_densVelsSigmasXYZ_CC_R24_z0p5[200:]=GRHOD_densVelsSigmas_CC_R24_z0p5
    
        
    F5HOD_muXYZProfZ_CC_R24_z0p5=[[] for i in range(300)]
    F5HOD_muXYZProfZ_CC_R24_z0p5[:100]=F5HOD_muXProfZ_CC_R24_z0p5
    F5HOD_muXYZProfZ_CC_R24_z0p5[100:200]=F5HOD_muYProfZ_CC_R24_z0p5
    F5HOD_muXYZProfZ_CC_R24_z0p5[200:]=F5HOD_muZProfZ_CC_R24_z0p5
    
    F5HOD_voidInfoXYZ_CC_R24_z0p5=[[] for i in range(300)]
    F5HOD_voidInfoXYZ_CC_R24_z0p5[:100]=F5HOD_voidInfo_CC_R24_z0p5
    F5HOD_voidInfoXYZ_CC_R24_z0p5[100:200]=F5HOD_voidInfo_CC_R24_z0p5
    F5HOD_voidInfoXYZ_CC_R24_z0p5[200:]=F5HOD_voidInfo_CC_R24_z0p5

    F5HOD_densVelsSigmasXYZ_CC_R24_z0p5=[[] for i in range(300)]
    F5HOD_densVelsSigmasXYZ_CC_R24_z0p5[:100]=F5HOD_densVelsSigmas_CC_R24_z0p5
    F5HOD_densVelsSigmasXYZ_CC_R24_z0p5[100:200]=F5HOD_densVelsSigmas_CC_R24_z0p5
    F5HOD_densVelsSigmasXYZ_CC_R24_z0p5[200:]=F5HOD_densVelsSigmas_CC_R24_z0p5
    
        
    N1HOD_muXYZProfZ_CC_R24_z0p5=[[] for i in range(300)]
    N1HOD_muXYZProfZ_CC_R24_z0p5[:100]=N1HOD_muXProfZ_CC_R24_z0p5
    N1HOD_muXYZProfZ_CC_R24_z0p5[100:200]=N1HOD_muYProfZ_CC_R24_z0p5
    N1HOD_muXYZProfZ_CC_R24_z0p5[200:]=N1HOD_muZProfZ_CC_R24_z0p5
    
    N1HOD_voidInfoXYZ_CC_R24_z0p5=[[] for i in range(300)]
    N1HOD_voidInfoXYZ_CC_R24_z0p5[:100]=N1HOD_voidInfo_CC_R24_z0p5
    N1HOD_voidInfoXYZ_CC_R24_z0p5[100:200]=N1HOD_voidInfo_CC_R24_z0p5
    N1HOD_voidInfoXYZ_CC_R24_z0p5[200:]=N1HOD_voidInfo_CC_R24_z0p5

    N1HOD_densVelsSigmasXYZ_CC_R24_z0p5=[[] for i in range(300)]
    N1HOD_densVelsSigmasXYZ_CC_R24_z0p5[:100]=N1HOD_densVelsSigmas_CC_R24_z0p5
    N1HOD_densVelsSigmasXYZ_CC_R24_z0p5[100:200]=N1HOD_densVelsSigmas_CC_R24_z0p5
    N1HOD_densVelsSigmasXYZ_CC_R24_z0p5[200:]=N1HOD_densVelsSigmas_CC_R24_z0p5
    
    
    


    

    
    

    print("beginning assembly")

    "############# Assemble catalogs here ############# - option 1"
    # start=time.time()
    # catGRHODXYZ_CC_R24_z0p5=catalogClass(voidInfos=GRHOD_voidInfo_CC_R24_z0p5,partsRandMu=GRHOD_muXYZProfZ_CC_R24_z0p5,densVelsSigmas=GRHOD_densVelsSigmas_CC_R24_z0p5,scaling="R24",centerType="CC",gravType="GRHOD")
    # catF5HODXYZ_CC_R24_z0p5=catalogClass(voidInfos=F5HOD_voidInfo_CC_R24_z0p5,partsRandMu=F5HOD_muXYZProfZ_CC_R24_z0p5,densVelsSigmas=F5HOD_densVelsSigmas_CC_R24_z0p5,scaling="R24",centerType="CC",gravType="F5HOD")
    # catN1HODXYZ_CC_R24_z0p5=catalogClass(voidInfos=N1HOD_voidInfo_CC_R24_z0p5,partsRandMu=N1HOD_muXYZProfZ_CC_R24_z0p5,densVelsSigmas=N1HOD_densVelsSigmas_CC_R24_z0p5,scaling="R24",centerType="CC",gravType="N1HOD")
    
    
    "############# Assemble catalogs here ############# - option 2 - used in paper"
    # start=time.time()
    catGRHODXYZ_CC_R24_z0p5=catalogClass(voidInfos=GRHOD_voidInfoXYZ_CC_R24_z0p5,partsRandMu=GRHOD_muXYZProfZ_CC_R24_z0p5,densVelsSigmas=GRHOD_densVelsSigmasXYZ_CC_R24_z0p5,scaling="R24",centerType="CC",gravType="GRHOD")
    catF5HODXYZ_CC_R24_z0p5=catalogClass(voidInfos=F5HOD_voidInfoXYZ_CC_R24_z0p5,partsRandMu=F5HOD_muXYZProfZ_CC_R24_z0p5,densVelsSigmas=F5HOD_densVelsSigmasXYZ_CC_R24_z0p5,scaling="R24",centerType="CC",gravType="F5HOD")
    catN1HODXYZ_CC_R24_z0p5=catalogClass(voidInfos=N1HOD_voidInfoXYZ_CC_R24_z0p5,partsRandMu=N1HOD_muXYZProfZ_CC_R24_z0p5,densVelsSigmas=N1HOD_densVelsSigmasXYZ_CC_R24_z0p5,scaling="R24",centerType="CC",gravType="N1HOD")
    
    Q=0















