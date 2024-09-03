#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 12:00:54 2023

@author: christopherwilson

These blocks of code will create all of the figures in the main section of arXiv:2212.02569. 


"""


from voids_MG_GSM_methods import *

#%%


colors_dic = {"black": "#000000", "orange": "#ff7f00", "skyblue": "#56B4E9", "green": "#4daf4a", "yellow": "#F0E442", "blue": "#377eb8", "red": "#D55E00", "purple":"#CC79A7"}




cmap_wBLACK = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
cmap =  [ "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]



#%%

"Fig 1 - void size function for Halos and HOD's, along with relative error wrt GR"
if True:
    
    sns.set_theme(style="whitegrid")
    numCols=2
    numRows=1
    scale=1.6
    
    fig, axs = plt.subplots(nrows=numRows, ncols=numCols, sharex=False, sharey=False, figsize=(4*numCols*scale,0.7*3*numRows*scale))
    
    
    halos_nBins=26
    halos_rMax=65
    
    fac=10**7
    
    
    "Halos"
    GR_VSF_vol=catGR_CC_R24_z0p5.jackKnifeVSF_VolNorm(rMax=halos_rMax,nBins=halos_nBins)
    F5_VSF_vol=catF5_CC_R24_z0p5.jackKnifeVSF_VolNorm(rMax=halos_rMax,nBins=halos_nBins)
    N1_VSF_vol=catN1_CC_R24_z0p5.jackKnifeVSF_VolNorm(rMax=halos_rMax,nBins=halos_nBins)

    axs[0].errorbar(GR_VSF_vol[0],GR_VSF_vol[1]*fac,GR_VSF_vol[2]*10*fac,capsize=2,linewidth=1,color=colors_dic['blue'])
    axs[0].errorbar(F5_VSF_vol[0],F5_VSF_vol[1]*fac,F5_VSF_vol[2]*10*fac,capsize=2,linewidth=1,color=colors_dic['orange'])
    axs[0].errorbar(N1_VSF_vol[0],N1_VSF_vol[1]*fac,F5_VSF_vol[2]*10*fac,capsize=2,linewidth=1,color=colors_dic['green'])
    

    axs[0].legend(["GR - Halos","F5 - Halos","N1 - Halos"],fontsize=12)
    

    "HOD"
    GRHOD_VSF_vol=catGRHOD_CC_R24_z0p5.jackKnifeVSF_VolNorm(nBins=round(80/halos_rMax*halos_nBins))
    F5HOD_VSF_vol=catF5HOD_CC_R24_z0p5.jackKnifeVSF_VolNorm(nBins=round(80/halos_rMax*halos_nBins))
    N1HOD_VSF_vol=catN1HOD_CC_R24_z0p5.jackKnifeVSF_VolNorm(nBins=round(80/halos_rMax*halos_nBins))

    axs[1].errorbar(GRHOD_VSF_vol[0],GRHOD_VSF_vol[1]*fac,GRHOD_VSF_vol[2]*10*fac,capsize=2,linewidth=1,linestyle="-",color=colors_dic['blue'])
    axs[1].errorbar(F5HOD_VSF_vol[0],F5HOD_VSF_vol[1]*fac,F5HOD_VSF_vol[2]*10*fac,capsize=2,linewidth=1,linestyle="-",color=colors_dic['orange'])
    axs[1].errorbar(N1HOD_VSF_vol[0],N1HOD_VSF_vol[1]*fac,F5HOD_VSF_vol[2]*10*fac,capsize=2,linewidth=1,linestyle="-",color=colors_dic['green'])
    axs[1].legend(["GR - HOD","F5 - HOD","N1 - HOD"],fontsize=12)
    
    

    axs[0].set_xlabel(r"R $[Mpc/h]$",fontsize=12)
    axs[1].set_xlabel(r"R $[Mpc/h]$",fontsize=12)
    

    "xlims"
    
    axs[0].set_xlim([0,halos_rMax])
    axs[1].set_xlim([0,80])
    
    plt.text(-107.0,-0.004,r"$ \frac{d\ n(R_{eff} < R)}{d\ r}$",fontsize=16,rotation = "vertical")
    plt.text(-106.0,0.55,r"$[10^{-7}(Mpc/h)^{-4}]$",fontsize = 12,rotation = "vertical")
    

    # saveName="Fig_1_Final_cb"
    # plt.savefig(fname= plotDir + saveName + ".pdf",bbox_inches="tight")
    
    plt.show()
    plt.clf()

#%%
"########## completed ##########"
"Fig 3 - densities, velocities, and sigmas"
if True:
    
    
    sns.set_theme(style="whitegrid")
    numCols=3
    numRows=2
    scale=1.6
    
    fig, axs = plt.subplots(nrows=numRows, ncols=numCols, sharex=False, sharey=False, figsize=(4*numCols*scale,0.7*3*numRows*scale),gridspec_kw={'height_ratios':[1,1]})
    
    

    "#### densities ####"
    densGR_S=catGRHOD_CC_R24_z0p5.jackKnifeDensityProfile(rMin=0,rMax=35,numSamples=1000)
    densGR_L=catGRHOD_CC_R24_z0p5.jackKnifeDensityProfile(rMin=35,rMax=100,numSamples=1000)
    
    densF5_S=catF5HOD_CC_R24_z0p5.jackKnifeDensityProfile(rMin=0,rMax=35,numSamples=1000)
    densF5_L=catF5HOD_CC_R24_z0p5.jackKnifeDensityProfile(rMin=35,rMax=100,numSamples=1000)
    
    densN1_S=catN1HOD_CC_R24_z0p5.jackKnifeDensityProfile(rMin=0,rMax=35,numSamples=1000)
    densN1_L=catN1HOD_CC_R24_z0p5.jackKnifeDensityProfile(rMin=35,rMax=100,numSamples=1000)
    
    
    "#### radial velocities ####"
    velGR_S=catGRHOD_CC_R24_z0p5.jackKnifeVelocityProfile(rMin=0,rMax=35,numSamples=1000,weight=0)
    velGR_L=catGRHOD_CC_R24_z0p5.jackKnifeVelocityProfile(rMin=35,rMax=100,numSamples=1000,weight=0)
    
    velF5_S=catF5HOD_CC_R24_z0p5.jackKnifeVelocityProfile(rMin=0,rMax=35,numSamples=1000,weight=0)
    velF5_L=catF5HOD_CC_R24_z0p5.jackKnifeVelocityProfile(rMin=35,rMax=100,numSamples=1000,weight=0)
    
    velN1_S=catN1HOD_CC_R24_z0p5.jackKnifeVelocityProfile(rMin=0,rMax=35,numSamples=1000,weight=0)
    velN1_L=catN1HOD_CC_R24_z0p5.jackKnifeVelocityProfile(rMin=35,rMax=100,numSamples=1000,weight=0)
    
    
    
    "#### sigmas ####"
    sigmasGR_S=catGRHOD_CC_R24_z0p5.calcSigmaV_r(rMin=0,rMax=35,weight=0)
    sigmasGR_L=catGRHOD_CC_R24_z0p5.calcSigmaV_r(rMin=35,rMax=100,weight=0)
    
    sigmasF5_S=catF5HOD_CC_R24_z0p5.calcSigmaV_r(rMin=0,rMax=35,weight=0)
    sigmasF5_L=catF5HOD_CC_R24_z0p5.calcSigmaV_r(rMin=35,rMax=100,weight=0)
    
    sigmasN1_S=catN1HOD_CC_R24_z0p5.calcSigmaV_r(rMin=0,rMax=35,weight=0)
    sigmasN1_L=catN1HOD_CC_R24_z0p5.calcSigmaV_r(rMin=35,rMax=100,weight=0)
    


    xpoints=catGRHOD_CC_R24_z0p5.xpoints
    
   
    "GR"
    axs[0][0].errorbar(densGR_L[0],1+densGR_L[1],densGR_L[2]*0,color=colors_dic['blue'])
    axs[0][0].errorbar(densGR_S[0],1+densGR_S[1],densGR_S[2]*0,color=colors_dic['blue'],linestyle="--")

    axs[0][1].errorbar(velGR_L[0][6:],velGR_L[1][6:],velGR_L[2][6:]*0,color=colors_dic['blue'])
    axs[0][1].errorbar(velGR_S[0][4:],velGR_S[1][4:],velGR_S[2][4:]*0,color=colors_dic['blue'],linestyle="--")
    
    axs[0][2].plot(xpoints[6:],sigmasGR_L[6:],color=colors_dic['blue'])
    axs[0][2].plot(xpoints[4:],sigmasGR_S[4:],color=colors_dic['blue'],linestyle="--")
    
    
    "F5"
    axs[1][0].errorbar(densF5_L[0],densF5_L[1] - densGR_L[1],densF5_L[2]*0,color=colors_dic['orange'])
    axs[1][0].errorbar(densF5_S[0],densF5_S[1] - densGR_S[1],densF5_S[2]*0,color=colors_dic['orange'],linestyle="--")
    
    axs[1][1].errorbar(velF5_L[0][6:],velF5_L[1][6:] - velGR_L[1][6:],velF5_L[2][6:]*0,color=colors_dic['orange'])
    axs[1][1].errorbar(velF5_S[0][4:],velF5_S[1][4:] - velGR_S[1][4:],velF5_S[2][4:]*0,color=colors_dic['orange'],linestyle="--")
    
    axs[1][2].plot(xpoints[6:],sigmasF5_L[6:] - sigmasGR_L[6:],color=colors_dic['orange'])
    axs[1][2].plot(xpoints[4:],sigmasF5_S[4:] - sigmasGR_S[4:],color=colors_dic['orange'],linestyle="--")
    
    
    "N1"
    axs[1][0].errorbar(densN1_L[0],densN1_L[1] - densGR_L[1],densGR_L[2]*0,color=colors_dic['green'])
    axs[1][0].errorbar(densN1_S[0],densN1_S[1] - densGR_S[1],densN1_S[2]*0,color=colors_dic['green'],linestyle="--")
    
    axs[1][1].errorbar(velN1_L[0][6:],velN1_L[1][6:] - velGR_L[1][6:],velN1_L[2][6:]*0,color=colors_dic['green'])
    axs[1][1].errorbar(velN1_S[0][4:],velN1_S[1][4:] - velGR_S[1][4:],velN1_S[2][4:]*0,color=colors_dic['green'],linestyle="--")
    
    axs[1][2].plot(xpoints[6:],sigmasN1_L[6:] - sigmasGR_L[6:],color=colors_dic['green'])
    axs[1][2].plot(xpoints[4:],sigmasN1_S[4:] - sigmasGR_S[4:],color=colors_dic['green'],linestyle="--")
     
    
    
    "grey lines to signal no particles interior to a point"
    axs[0][1].plot([xpoints[4],xpoints[4]],[-1000,2000],color="grey",linestyle="--",label="_nolegend_")
    axs[0][2].plot([xpoints[4],xpoints[4]],[-1000,2000],color="grey",linestyle="--",label="_nolegend_")
    
    axs[1][1].plot([xpoints[4],xpoints[4]],[-1000,2000],color="grey",linestyle="--",label="_nolegend_")
    axs[1][2].plot([xpoints[4],xpoints[4]],[-1000,2000],color="grey",linestyle="--",label="_nolegend_")
    
    
    axs[0][1].plot([xpoints[6],xpoints[6]],[-1000,2000],color="grey",linestyle="-",label="_nolegend_")
    axs[0][2].plot([xpoints[6],xpoints[6]],[-1000,2000],color="grey",linestyle="-",label="_nolegend_")
    
    axs[1][1].plot([xpoints[6],xpoints[6]],[-1000,2000],color="grey",linestyle="-",label="_nolegend_")
    axs[1][2].plot([xpoints[6],xpoints[6]],[-1000,2000],color="grey",linestyle="-",label="_nolegend_")
    
    
    axs[1][0].set_xlabel(r"$r$ [Mpc/h]",fontsize=14)
    axs[1][1].set_xlabel(r"$r$ [Mpc/h]",fontsize=14)
    axs[1][2].set_xlabel(r"$r$ [Mpc/h]",fontsize=14)
    
    
    
    axs[0][0].set_ylabel(r"$1 + \delta^{r}$",fontsize=14)
    axs[0][1].set_ylabel(r"$v_{r}$ [km/s]",fontsize=14)
    axs[0][2].set_ylabel(r"$\sigma_{v_{\parallel}}$ [km/s]",fontsize=14)
    
    
    axs[1][0].set_ylabel(r"$\delta^{r}_{MG} - \delta^{r}_{GR}$",fontsize=14)
    axs[1][1].set_ylabel(r"$v_{r,MG} - v_{r,GR}$ [km/s]",fontsize=14)
    axs[1][2].set_ylabel(r"$\sigma_{v_{\parallel},MG} - \sigma_{v_{\parallel},GR}$ [km/s]",fontsize=14)
    
    

    axs[0][0].set_xlim([0,80])
    axs[0][1].set_xlim([0,80])
    axs[0][2].set_xlim([0,80])
    axs[1][0].set_xlim([0,80])
    axs[1][1].set_xlim([0,80])
    axs[1][2].set_xlim([0,80])
    
    

    axs[0][0].set_ylim([0,1.4])
    axs[0][1].set_ylim([-100,225])
    axs[0][2].set_ylim([200,350])
    
    axs[1][0].set_ylim([-0.0175,0.01])
    axs[1][1].set_ylim([-10,30])
    axs[1][2].set_ylim([-0,40])
    
    

    
    axs[0][0].legend([r"GR, $R_{eff} \geq$ 35 Mpc/h",r"GR, $R_{eff} <$ 35 Mpc/h"],fontsize=10,loc="best")



    axs[1][0].legend([r"F5 - GR, $R_{eff} \geq$ 35 Mpc/h",r"F5 - GR, $R_{eff} <$ 35 Mpc/h",
                      r"N1 - GR, $R_{eff} \geq$ 35 Mpc/h",r"N1 - GR, $R_{eff} <$ 35 Mpc/h"],fontsize=10,loc="best")
    



    # saveName="Fig_3_Final_cb"

    # plt.savefig(fname= plotDir + saveName + ".pdf",bbox_inches="tight")
    

    plt.show()
    plt.clf()
    
    


#%%
"########## completed ##########" 
"Fig 2 - Evolution across scale of quantities, and effects of the HOD"
if True:
    sns.set_theme(style="whitegrid")
    
    numCols=2
    numRows=1
    scale=1.6
    
    fig, axs = plt.subplots(nrows=numRows, ncols=numCols, sharex=False, sharey=False, figsize=(4*numCols*scale,0.7*3*numRows*scale))
    

    "Halos"
    rMinList=np.array([10,20,30,40,50,60])
    rMaxList=np.array([20,30,40,50,60,70])
    
    
    rMedianArray=np.array([15,25,35,45,55,65])
    

    
    
    velProfArrayGR_R24=[catGR_CC_R24_z0p5.jackKnifeVelocityProfile(rMin=rMinList[i],rMax=rMaxList[i],weight=0) for i in range(len(rMinList))]
    velProfArrayF5_R24=[catF5_CC_R24_z0p5.jackKnifeVelocityProfile(rMin=rMinList[i],rMax=rMaxList[i],weight=0) for i in range(len(rMinList))]
    velProfArrayN1_R24=[catN1_CC_R24_z0p5.jackKnifeVelocityProfile(rMin=rMinList[i],rMax=rMaxList[i],weight=0) for i in range(len(rMinList))]
    
    

    
    velMaxArrayGR_R24=[np.max(velProfArrayGR_R24[i][1]) for i in range(len(rMinList))]
    velMaxArrayF5_R24=[np.max(velProfArrayF5_R24[i][1]) for i in range(len(rMinList))]
    velMaxArrayN1_R24=[np.max(velProfArrayN1_R24[i][1]) for i in range(len(rMinList))]
    
    

    
    velRatioArrayF5_R24=[velMaxArrayF5_R24[i]/velMaxArrayGR_R24[i]for i in range(len(rMinList))]
    velRatioArrayN1_R24=[velMaxArrayN1_R24[i]/velMaxArrayGR_R24[i]for i in range(len(rMinList))]
    
    
    axs[0].plot(rMedianArray[:-1],velRatioArrayF5_R24[:-1],color=colors_dic['orange'],linestyle="-",marker="*")
    axs[0].plot(rMedianArray[:-1],velRatioArrayN1_R24[:-1],color=colors_dic['green'],linestyle="-",marker="*")
    

    axs[0].set_ylabel(r"$v_{r,peak,MG} \ /\  v_{r,peak,GR}$",fontsize = 12)
    axs[0].set_xlabel(r"$R_{eff}$ [Mpc/h]",fontsize=14)
    axs[0].set_yticks([1,1.05,1.1,1.15,1.2])




    "HOD"
    rMinList=np.array([10,20,30,40,50,60]) 
    rMaxList=np.array([20,30,40,50,60,70])
    
    
    rMedianArray=np.array([15,25,35,45,55,65])
    
    velProfArrayGRHOD_R24=[catGRHOD_CC_R24_z0p5.jackKnifeVelocityProfile(rMin=rMinList[i],rMax=rMaxList[i],weight=0) for i in range(len(rMinList))]
    velProfArrayF5HOD_R24=[catF5HOD_CC_R24_z0p5.jackKnifeVelocityProfile(rMin=rMinList[i],rMax=rMaxList[i],weight=0) for i in range(len(rMinList))]
    velProfArrayN1HOD_R24=[catN1HOD_CC_R24_z0p5.jackKnifeVelocityProfile(rMin=rMinList[i],rMax=rMaxList[i],weight=0) for i in range(len(rMinList))]
    
    

    velMaxArrayGRHOD_R24=[np.max(velProfArrayGRHOD_R24[i][1]) for i in range(len(rMinList))]
    velMaxArrayF5HOD_R24=[np.max(velProfArrayF5HOD_R24[i][1]) for i in range(len(rMinList))]
    velMaxArrayN1HOD_R24=[np.max(velProfArrayN1HOD_R24[i][1]) for i in range(len(rMinList))]
    
    
    
    velRatioArrayF5HOD_R24=[velMaxArrayF5HOD_R24[i]/velMaxArrayGRHOD_R24[i]for i in range(len(rMinList))]
    velRatioArrayN1HOD_R24=[velMaxArrayN1HOD_R24[i]/velMaxArrayGRHOD_R24[i]for i in range(len(rMinList))]
    

    
    axs[0].plot(rMedianArray[:-1],velRatioArrayF5HOD_R24[:-1],color=colors_dic['orange'],linestyle="--",marker="*")
    axs[0].plot(rMedianArray[:-1],velRatioArrayN1HOD_R24[:-1],color=colors_dic['green'],linestyle="--",marker="*")
    


    axs[0].set_yticks([1,1.05,1.1,1.15,1.2])



    "Sigmas from Halos"
    
    rMinList=[10,20,30,40,50,60]
    rMaxList=[20,30,40,50,60,70]
    
    
    rMedianArray=np.array([15,25,35,45,55,65])
    
    sigmaProfArrayGR_R24=[catGR_CC_R24_z0p5.calcSigmaV_r(rMin=rMinList[i],rMax=rMaxList[i],weight=0) for i in range(len(rMinList))]
    sigmaProfArrayF5_R24=[catF5_CC_R24_z0p5.calcSigmaV_r(rMin=rMinList[i],rMax=rMaxList[i],weight=0) for i in range(len(rMinList))]
    sigmaProfArrayN1_R24=[catN1_CC_R24_z0p5.calcSigmaV_r(rMin=rMinList[i],rMax=rMaxList[i],weight=0) for i in range(len(rMinList))]
        
    sigmaAvgArrayGR_R24=[np.mean(sigmaProfArrayGR_R24[i][25:]) for i in range(len(rMinList))]
    sigmaAvgArrayF5_R24=[np.mean(sigmaProfArrayF5_R24[i][25:]) for i in range(len(rMinList))]
    sigmaAvgArrayN1_R24=[np.mean(sigmaProfArrayN1_R24[i][25:]) for i in range(len(rMinList))]
    
    
    sigmaRatioArrayF5_R24=[sigmaAvgArrayF5_R24[i]/sigmaAvgArrayGR_R24[i]for i in range(len(rMinList))]
    sigmaRatioArrayN1_R24=[sigmaAvgArrayN1_R24[i]/sigmaAvgArrayGR_R24[i]for i in range(len(rMinList))]
    
    sigmaRatioArrayF5_R24=[sigmaAvgArrayF5_R24[i]/sigmaAvgArrayGR_R24[i]for i in range(len(rMinList))]
    sigmaRatioArrayN1_R24=[sigmaAvgArrayN1_R24[i]/sigmaAvgArrayGR_R24[i]for i in range(len(rMinList))]
    
    
    axs[1].plot(rMedianArray[:-1],sigmaRatioArrayF5_R24[:-1],color=colors_dic['orange'],linestyle="-",marker="*")
    axs[1].plot(rMedianArray[:-1],sigmaRatioArrayN1_R24[:-1],color=colors_dic['green'],linestyle="-",marker="*")
    


    "Sigmas from HOD"
    
    rMinList=[10,20,30,40,50,60]
    rMaxList=[20,30,40,50,60,70]
    
    
    rMedianArray=np.array([15,25,35,45,55,65])
    
    sigmaProfArrayGRHOD_R24=[catGRHOD_CC_R24_z0p5.calcSigmaV_r(rMin=rMinList[i],rMax=rMaxList[i],weight=0) for i in range(len(rMinList))]
    sigmaProfArrayF5HOD_R24=[catF5HOD_CC_R24_z0p5.calcSigmaV_r(rMin=rMinList[i],rMax=rMaxList[i],weight=0) for i in range(len(rMinList))]
    sigmaProfArrayN1HOD_R24=[catN1HOD_CC_R24_z0p5.calcSigmaV_r(rMin=rMinList[i],rMax=rMaxList[i],weight=0) for i in range(len(rMinList))]
    

    sigmaAvgArrayGRHOD_R24=[np.mean(sigmaProfArrayGRHOD_R24[i][25:])+1 for i in range(len(rMinList))]
    sigmaAvgArrayF5HOD_R24=[np.mean(sigmaProfArrayF5HOD_R24[i][25:])+1 for i in range(len(rMinList))]
    sigmaAvgArrayN1HOD_R24=[np.mean(sigmaProfArrayN1HOD_R24[i][25:])+1 for i in range(len(rMinList))]
    
    
    sigmaRatioArrayF5HOD_R24=[sigmaAvgArrayF5HOD_R24[i]/sigmaAvgArrayGRHOD_R24[i]for i in range(len(rMinList))]
    sigmaRatioArrayN1HOD_R24=[sigmaAvgArrayN1HOD_R24[i]/sigmaAvgArrayGRHOD_R24[i]for i in range(len(rMinList))]
    

    
    axs[1].plot(rMedianArray[:-1],sigmaRatioArrayF5HOD_R24[:-1],color=colors_dic['orange'],linestyle="--",marker="*")
    axs[1].plot(rMedianArray[:-1],sigmaRatioArrayN1HOD_R24[:-1],color=colors_dic['green'],linestyle="--",marker="*")
    

    axs[1].legend([r"F5 / GR, - Halos",r"N1 / GR, - Halos",r"F5 / GR, - HOD",r"N1 / GR, - HOD"],fontsize=11.5,ncol=2,loc=[0.02,0.05])
    axs[1].set_ylabel(r"$\sigma_{v_{\parallel},MG} \ /\  \sigma_{v_{\parallel},GR}$",fontsize = 12)
    axs[1].set_xlabel(r"$R_{eff}$ [Mpc/h]",fontsize=14)
    axs[1].set_ylim([1.0,1.11])
    axs[0].set_ylim([1.0,1.2])
    
    
    # saveName="Fig_2_Final_cb"
    # plt.savefig(fname= plotDir + saveName + ".pdf",bbox_inches="tight")
    
    plt.show()
    plt.clf()
    


#%%
"########## completed ##########"
"Fig 4 - quadData, without errors"
if True:
    
    
    "here"
    sns.set_theme(style="whitegrid")
    saveName="Fig_3"    

    numCols=2
    numRows=2
    scale=1.6
    
    fig, axs = plt.subplots(nrows=numRows, ncols=numCols, sharex=False, sharey=False, figsize=(4*numCols*scale,0.7*3*numRows*scale),gridspec_kw={'height_ratios':[3,1]})
    



    quadData_GRHOD=catGRHODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=35,rMax=100)
    quadData_F5HOD=catF5HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=35,rMax=100)
    quadData_N1HOD=catN1HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=35,rMax=100)
    
    
    quadData_GRHOD_S=catGRHODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=0,rMax=35)
    quadData_F5HOD_S=catF5HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=0,rMax=35)
    quadData_N1HOD_S=catN1HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=0,rMax=35)
    
        
    errorGR=quadData_GRHOD[2]
    errorF5=quadData_F5HOD[2]
    errorN1=quadData_N1HOD[2]
    
    
    endIndex=30
    errorEveryLim=2
    
    
    "errorbars"

    axs[0][0].errorbar(quadData_GRHOD[0],quadData_GRHOD[1],color=colors_dic['blue'],capsize=2,errorevery=np.append(np.arange(3,endIndex,errorEveryLim),[5]),marker="o",markersize=0,zorder=10,linewidth=1)
    axs[0][0].errorbar(quadData_F5HOD[0],quadData_F5HOD[1],color=colors_dic['orange'],capsize=2,errorevery=np.append(np.arange(4,endIndex,errorEveryLim),[6]),marker="^",markersize=0,zorder=11,linewidth=1)



    axs[0][1].errorbar(quadData_GRHOD[0],quadData_GRHOD[1],color=colors_dic['blue'],capsize=2,errorevery=np.append(np.arange(3,endIndex,errorEveryLim),[1]),marker="o",markersize=0,zorder=10,linewidth=1)
    axs[0][1].errorbar(quadData_N1HOD[0],quadData_N1HOD[1],color=colors_dic['green'],capsize=2,errorevery=np.append(np.arange(4,endIndex,errorEveryLim),[1]),marker="^",markersize=0,zorder=11,linewidth=1)


    axs[1][0].plot(quadData_GRHOD[0],(quadData_F5HOD[1] - quadData_GRHOD[1]),color=colors_dic['orange'],linewidth=1)
    axs[1][1].plot(quadData_GRHOD[0],(quadData_N1HOD[1] - quadData_GRHOD[1]),color=colors_dic['green'],linewidth=1)


    axs[0][0].errorbar(quadData_GRHOD_S[0],quadData_GRHOD_S[1],color=colors_dic['blue'],capsize=2,errorevery=np.append(np.arange(1,endIndex,errorEveryLim),[1,3,7]),linestyle="-.",marker="o",markersize=0,zorder=5,linewidth=1)
    axs[0][0].errorbar(quadData_F5HOD_S[0],quadData_F5HOD_S[1],color=colors_dic['orange'],capsize=2,errorevery=np.append(np.arange(2,endIndex,errorEveryLim),[2,4,8]),linestyle="-.",marker="^",markersize=0,zorder=6,linewidth=1)


    axs[0][1].errorbar(quadData_GRHOD_S[0],quadData_GRHOD_S[1],color=colors_dic['blue'],capsize=2,errorevery=np.append(np.arange(1,endIndex,errorEveryLim),[1,3,7]),linestyle="-.",marker="o",markersize=0,zorder=5,linewidth=1)
    axs[0][1].errorbar(quadData_N1HOD_S[0],quadData_N1HOD_S[1],color=colors_dic['green'],capsize=2,errorevery=np.append(np.arange(2,endIndex,errorEveryLim),[2,4,8]),linestyle="-.",marker="^",markersize=0,zorder=6,linewidth=1)

    axs[1][0].plot(quadData_GRHOD_S[0],(quadData_F5HOD_S[1] - quadData_GRHOD_S[1]),color=colors_dic['orange'],linestyle="--",linewidth=1)
    axs[1][1].plot(quadData_GRHOD_S[0],(quadData_N1HOD_S[1] - quadData_GRHOD_S[1]),color=colors_dic['green'],linestyle="--",linewidth=1)



    axs[0][0].set_ylabel(r"$\xi_{2}$",fontsize=14)
    axs[1][0].set_xlabel(r"$s$ [Mpc/h]",fontsize=14)
    axs[1][1].set_xlabel(r"$s$ [Mpc/h]",fontsize=14)
    
    axs[1][0].set_ylabel(r"$\Delta \xi_{2}$",fontsize=14)
    
   

    axs[0][0].legend([r"GR, $R_{eff} \geq 35$ Mpc/h",r"F5, $R_{eff} \geq 35$ Mpc/h",
                      r"GR, $R_{eff} < 35$ Mpc/h",r"F5, $R_{eff} < 35$ Mpc/h"],fontsize=12)
    axs[0][1].legend([r"GR, $R_{eff} \geq 35$ Mpc/h",r"N1, $R_{eff} \geq 35$ Mpc/h",
                      r"GR, $R_{eff} < 35$ Mpc/h",r"N1, $R_{eff} < 35$ Mpc/h"],fontsize=12)
                  
    
    axs[0][0].set_ylim([-0.175,0.09])
    axs[0][1].set_ylim([-0.175,0.09])
    
    axs[1][0].set_ylim([-0.015,0.015])
    axs[1][1].set_ylim([-0.015,0.015])
    
    
    endLoc=80
    axs[0][0].set_xlim([0,endLoc])
    axs[0][1].set_xlim([0,endLoc])
    
    axs[1][0].set_xlim([0,endLoc])
    axs[1][1].set_xlim([0,endLoc])
    
    
    

    
    # saveName="Fig_4_Final_cb"
    # plt.savefig(fname= plotDir + saveName + ".pdf",bbox_inches="tight")
    
    plt.show()
    plt.clf()
    


#%%
"########## completed ##########"
"Fig 5 - Unscaled Data vs. GSM, Large and Small Voids, HOD, F5 and N1, with relative errors shown below"

if True:
    "here"
    sns.set_theme(style="whitegrid")
    
    rMin=35
    rMax=100
    voidTypes=["R","S"]
    
    small=True


    numCols=2
    numRows=2
    scale=1.6
    
    fig, axs = plt.subplots(nrows=numRows, ncols=numCols, sharex=False, sharey=False, figsize=(4*numCols*scale,0.7*3*numRows*scale),gridspec_kw={'height_ratios':[3,1]})
    
    

    quadThy_F5HOD_R24_L=catF5HODXYZ_CC_R24_z0p5.calcFullGSM(rMin=35,rMax=100,numBinFac=10)
    quadThy_N1HOD_R24_L=catN1HODXYZ_CC_R24_z0p5.calcFullGSM(rMin=35,rMax=100,numBinFac=10)
    
    
    quadThy_F5HOD_R24_S=catF5HODXYZ_CC_R24_z0p5.calcFullGSM(rMin=0,rMax=35,numBinFac=10)
    quadThy_N1HOD_R24_S=catN1HODXYZ_CC_R24_z0p5.calcFullGSM(rMin=0,rMax=35,numBinFac=10)
    

    quadData_F5HOD_L=catF5HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=35,rMax=100,numSamples=2000)
    quadData_N1HOD_L=catN1HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=35,rMax=100,numSamples=2000)
    
    

    quadData_F5HOD_S=catF5HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=0,rMax=35,numSamples=2000)
    quadData_N1HOD_S=catN1HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=0,rMax=35,numSamples=2000)


    
    
    xpoints=catGRHODXYZ_CC_R24_z0p5.xpoints
    
    
    s=2


        

    quadThy_F5HOD_50_L=catF5HODXYZ_CC_R24_z0p5.calcFullGSM(rMin=35,rMax=100)
    quadThy_N1HOD_50_L=catN1HODXYZ_CC_R24_z0p5.calcFullGSM(rMin=35,rMax=100)
    

    quadThy_F5HOD_50_S=catF5HODXYZ_CC_R24_z0p5.calcFullGSM(rMin=0,rMax=35)
    quadThy_N1HOD_50_S=catN1HODXYZ_CC_R24_z0p5.calcFullGSM(rMin=0,rMax=35)
    
    Q=0
    

    
    axs[0][0].errorbar(quadData_F5HOD_L[0],quadData_F5HOD_L[1],quadData_F5HOD_L[2],color="black",capsize=4,linestyle="",marker="o",markersize=4)
    axs[0][0].errorbar(quadThy_F5HOD_R24_L[0],quadThy_F5HOD_R24_L[1],color=colors_dic['orange'],linestyle="-")
    
    axs[0][0].errorbar(quadData_F5HOD_S[0],quadData_F5HOD_S[1],quadData_F5HOD_S[2],color="blue",capsize=4,linestyle="",marker="o",markersize=4)
    axs[0][0].errorbar(quadThy_F5HOD_R24_S[0],quadThy_F5HOD_R24_S[1],color=colors_dic['orange'],linestyle="--")
    
    



    axs[0][1].errorbar(quadData_N1HOD_L[0],quadData_N1HOD_L[1],quadData_N1HOD_L[2],color="black",capsize=4,linestyle="",marker="o",markersize=4)
    axs[0][1].errorbar(quadThy_N1HOD_R24_L[0],quadThy_N1HOD_R24_L[1],color=colors_dic['green'],linestyle="-")

    axs[0][1].errorbar(quadData_N1HOD_S[0],quadData_N1HOD_S[1],quadData_N1HOD_S[2],color="blue",capsize=4,linestyle="",marker="o",markersize=4)
    axs[0][1].errorbar(quadThy_N1HOD_R24_S[0],quadThy_N1HOD_R24_S[1],color=colors_dic['green'],linestyle="--")



    axs[1][0].plot(quadData_F5HOD_L[0],(quadThy_F5HOD_50_L[1] - quadData_F5HOD_L[1])/(quadData_F5HOD_L[2]),color=colors_dic['orange'])
    axs[1][0].plot(quadData_F5HOD_S[0],(quadThy_F5HOD_50_S[1] - quadData_F5HOD_S[1])/(quadData_F5HOD_S[2]),color=colors_dic['orange'],linestyle="--")
    axs[1][0].plot([0,80],[0,0],color="black",linestyle="-",label="_nolegend_")
    axs[1][0].plot([0,80],[1,1],color="black",linestyle="--",label="_nolegend_")
    axs[1][0].plot([0,80],[-1,-1],color="black",linestyle="--",label="_nolegend_")



    axs[1][1].plot(quadData_N1HOD_L[0],(quadThy_N1HOD_50_L[1] - quadData_N1HOD_L[1])/(quadData_N1HOD_L[2]),color=colors_dic['green'])
    axs[1][1].plot(quadData_N1HOD_S[0],(quadThy_N1HOD_50_S[1] - quadData_N1HOD_S[1])/(quadData_N1HOD_S[2]),color=colors_dic['green'],linestyle="--")
    
    axs[1][1].plot([0,80],[0,0],color="black",linestyle="-",label="_nolegend_")
    axs[1][1].plot([0,80],[1,1],color="black",linestyle="--",label="_nolegend_")
    axs[1][1].plot([0,80],[-1,-1],color="black",linestyle="--",label="_nolegend_")


    axs[0][0].set_ylabel(r"$\xi_{2}$",fontsize=14)
    axs[1][0].set_xlabel(r"$s$ [Mpc/h]",fontsize=14)
    axs[1][1].set_xlabel(r"$s$ [Mpc/h]",fontsize=14)
    
    axs[1][0].set_ylabel(r"$\frac{\Delta \xi_{2}}{ \sigma_{\xi_{2}}}$",fontsize=16)
    
    axs[0][0].legend([r"F5 - Data - $R_{eff}\ \geq$ 35 Mpc/h", r"F5 - GSM - $R_{eff}\ \geq$ 35 Mpc/h",
                      r"F5 - Data - $R_{eff}\ <$ 35 Mpc/h", r"F5 - GSM - $R_{eff}\ <$ 35 Mpc/h"],fontsize=12)
    
    
    axs[0][1].legend([r"N1 - Data - $R_{eff}\ \geq$ 35 Mpc/h", r"N1 - GSM - $R_{eff}\ \geq$ 35 Mpc/h",
                      r"N1 - Data - $R_{eff}\ <$ 35 Mpc/h ", r"N1 - GSM - $R_{eff}\ <$ 35 Mpc/h"],fontsize=12)
    
    
    
    axs[1][0].legend(["$R_{eff}\ \geq$ 35 Mpc/h","$R_{eff}\ <$ 35 Mpc/h"],fontsize=10,ncol=2,frameon=False,loc=(0.265,0.75))
    axs[1][1].legend(["$R_{eff}\ \geq$ 35 Mpc/h","$R_{eff}\ <$ 35 Mpc/h"],fontsize=10,ncol=2,frameon=False,loc=(0.265,0.75))
    

    
    axs[0][0].set_ylim([-0.175,0.07])
    axs[0][1].set_ylim([-0.175,0.07])
    
    axs[1][0].set_ylim([-4,4])
    axs[1][1].set_ylim([-4,4])
    axs[1][0].set_yticks([-4,-2,0,2,4])
    axs[1][1].set_yticks([-4,-2,0,2,4])
    
    dumb=[0,80]
    
    axs[0][0].set_xlim(dumb)
    axs[0][1].set_xlim(dumb)
    
    axs[1][0].set_xlim(dumb)
    axs[1][1].set_xlim(dumb)
    
    
    

    saveName="Fig_5_Final_cb"

    plt.show()
    plt.clf()




#%%
"########## completed ##########"
"Fig 6 -- Functional Derivatives"

if True:
    
    sns.set_theme(style="whitegrid")


    
    sns.set_theme(style="whitegrid")
    saveName="Fig_6_v2_sns" 

    numCols=2
    numRows=2
    scale=1.6
    
    fig, axs = plt.subplots(nrows=numRows, ncols=numCols, sharex=False, sharey=False, figsize=(4*numCols*scale,0.7*3*numRows*scale),gridspec_kw={'height_ratios':[3,1]})
    
    ax = fig.add_subplot(111,facecolor="none") 
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    ax.set_yticks([])
    ax.set_xticks([])
    


    Spoints=(np.linspace(0,120,301)[0:-1] + np.linspace(0,120,301)[1:])/2
    quadDelta_delta_F5HOD_R24_S=calcDeltaXi2_delta(catGRHOD_CC_R24_z0p5,catF5HOD_CC_R24_z0p5,spoints=Spoints,rMin=0,rMax=35,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"])
    quadDelta_vel_F5HOD_R24_S=calcDeltaXi2_vel(catGRHOD_CC_R24_z0p5,catF5HOD_CC_R24_z0p5,spoints=Spoints,rMin=0,rMax=35,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"])
    quadDelta_sigma_F5HOD_R24_S=calcDeltaXi2_sigma(catGRHOD_CC_R24_z0p5,catF5HOD_CC_R24_z0p5,spoints=Spoints,rMin=0,rMax=35,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"])
        
    quadDelta_delta_N1HOD_R24_L=calcDeltaXi2_delta(catGRHOD_CC_R24_z0p5,catN1HOD_CC_R24_z0p5,spoints=Spoints,rMin=35,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"])
    quadDelta_vel_N1HOD_R24_L=calcDeltaXi2_vel(catGRHOD_CC_R24_z0p5,catN1HOD_CC_R24_z0p5,spoints=Spoints,rMin=35,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"])
    quadDelta_sigma_N1HOD_R24_L=calcDeltaXi2_sigma(catGRHOD_CC_R24_z0p5,catN1HOD_CC_R24_z0p5,spoints=Spoints,rMin=35,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"])
    
        
        
    endIndex=200
    saveName="Fig_6b_sns"
    axs[0][0].plot(quadDelta_delta_F5HOD_R24_S[0][:endIndex],quadDelta_delta_F5HOD_R24_S[1][:endIndex],color=colors_dic['blue'])
    axs[0][0].plot(quadDelta_vel_F5HOD_R24_S[0][:endIndex],quadDelta_vel_F5HOD_R24_S[1][:endIndex],color=colors_dic['orange'])
    axs[0][0].plot(quadDelta_sigma_F5HOD_R24_S[0][:endIndex],quadDelta_sigma_F5HOD_R24_S[1][:endIndex],color=colors_dic['green'])


    axs[0][1].plot(quadDelta_delta_N1HOD_R24_L[0][:endIndex],quadDelta_delta_N1HOD_R24_L[1][:endIndex],color=colors_dic['blue'])
    axs[0][1].plot(quadDelta_vel_N1HOD_R24_L[0][:endIndex],quadDelta_vel_N1HOD_R24_L[1][:endIndex],color=colors_dic['orange'])
    axs[0][1].plot(quadDelta_sigma_N1HOD_R24_L[0][:endIndex],quadDelta_sigma_N1HOD_R24_L[1][:endIndex],color=colors_dic['green'])
    optionB=True
    optionA=False
    
    
    axs[0][0].set_title(r"F5, $R_{eff}\ <$ 35 Mpc/h",fontsize=15,pad=10)
    axs[0][1].set_title(r"N1, $R_{eff}\ \geq$ 35 Mpc/h",fontsize=15,pad=10)
       


    axs[1][0].set_xlabel(r"$s$ [Mpc/h]",fontsize=14)
    axs[1][1].set_xlabel(r"$s$ [Mpc/h]",fontsize=14)
    
    axs[0][0].set_ylabel(r"$\Delta_{x}\xi_{2}$",fontsize=14)
    
    axs[0][0].set_ylim([-0.01,0.01])
    axs[0][1].set_ylim([-0.01,0.01])
    
    axs[0][0].set_xlim([0,80])
    axs[0][1].set_xlim([0,80])
    
    axs[1][0].set_xlim([0,80])
    axs[1][1].set_xlim([0,80])
    




    "#########################################################################################"

    quadData_GRHOD_L=catGRHODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=35,rMax=100)
    quadData_GRHOD_S=catGRHODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=0,rMax=35)
    
    quadData_F5HOD_S=catF5HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=0,rMax=35)
    quadData_N1HOD_L=catN1HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=35,rMax=100)
    
    xpoints=catGRHODXYZ_CC_R24_z0p5.xpoints
    

    startIndex=2
    SP=quadDelta_delta_N1HOD_R24_L[0][startIndex:]



    axs[1][0].set_ylabel(r"$\Delta\xi_{2}$",fontsize=14)
    
    totalF5_S=quadDelta_delta_F5HOD_R24_S[1][startIndex:] + quadDelta_vel_F5HOD_R24_S[1][startIndex:] + quadDelta_sigma_F5HOD_R24_S[1][startIndex:]
    axs[1][0].plot(SP,totalF5_S,color="black")
    axs[1][0].plot(xpoints,(quadData_F5HOD_S[1] - quadData_GRHOD_S[1]),color="violet")
    
    
    axs[1][0].legend([r"$\Delta \xi_{2}$ - GSM" ,r"$\Delta \xi_{2}$ - Data"],ncol=2,loc=4,frameon=False,fontsize=11.35)
    axs[1][0].set_ylim([-0.015,0.015])


    
    
    totalN1=quadDelta_delta_N1HOD_R24_L[1][startIndex:] + quadDelta_vel_N1HOD_R24_L[1][startIndex:] + quadDelta_sigma_N1HOD_R24_L[1][startIndex:]
    axs[1][1].plot(SP,totalN1,color="black")
    axs[1][1].plot(xpoints,(quadData_N1HOD_L[1] - quadData_GRHOD_L[1]),color="violet")
    
    axs[1][1].legend([r"$\Delta \xi_{2}$ - GSM" ,r"$\Delta \xi_{2}$ - Data"],ncol=2,loc=4,frameon=False,fontsize=11.35)
    axs[1][1].set_ylim([-0.015,0.015])

    

    "#########################################################################################"

    axs[0][0].legend([r"$\Delta_{\delta^{r}}\xi_{2,F5}$",
                r"$\Delta_{v_{r}}\xi_{2,F5}$",
                r"$\Delta_{\sigma_{v_{\parallel}}}\xi_{2,F5}$"],ncol=1,fontsize=15,loc="best",frameon=True)
    



    axs[0][1].legend([r"$\Delta_{\delta^{r}}\xi_{2,N1}$",
                r"$\Delta_{v_{r}}\xi_{2,N1}$",
                r"$\Delta_{\sigma_{v_{\parallel}}}\xi_{2,N1}$"],ncol=1,fontsize=15,loc="best",frameon=True)



    saveName="Fig_6_Final_cb"
    # plt.savefig(fname= plotDir + saveName + ".pdf",bbox_inches="tight")

#%%
"Fig_7"
"block to locally optimize and calculate all relevant quantities for confidence ellipses -- HOD - On one plot"
if True:
    
    "This plot takes a while to run..."
    error="new"
    numCols=1
    numRows=1
    scale=1.6
    fig, axs = plt.subplots(nrows=numRows, ncols=numCols, sharex=False, sharey=False, figsize=(4*numCols*scale,0.7*3*numRows*scale))
    
    sns.set_theme(style="whitegrid")
    


    
    "First, the R24 voids"
    
    
    rMin=35
    rMax=100
    startIndex=3
    
    shiftFac=1
    tol=.1
    XYZ=False
    
    save=False
    load=False
    loadHessian=False
    errors = error
    
    
    

    if False:
        catGRHOD_R24=catGRHODXYZ_CC_R24_z0p5
        catF5HOD_R24=catF5HODXYZ_CC_R24_z0p5
        catN1HOD_R24=catN1HODXYZ_CC_R24_z0p5
        
        
    if True:
        catGRHOD_R24=catGRHOD_CC_R24_z0p5
        catF5HOD_R24=catF5HOD_CC_R24_z0p5
        catN1HOD_R24=catN1HOD_CC_R24_z0p5
        
        
        



    endIndex=45
    ylimits=[235,340]
    
    "GRHOD"
    dAListGRHOD_R24,optParamsGRHOD_R24=catGRHOD_R24.calcConfidenceEllipses(startIndex=startIndex,endIndex=endIndex,rMin=rMin,rMax=rMax,errors=errors,load=load,XYZ=XYZ,shiftFac=shiftFac,save=save,tol=tol,loadHessian=loadHessian)
    betaGRHOD_R24=optimizeBetaFixedShift(catGRHOD_R24,rMin=rMin,rMax=rMax,plotting=False,shiftFac=shiftFac)
    dAListGRHOD_R24_95=catGRHOD_R24.dAList_95
                                                                                                                                  
  
    "F5HOD"
    dAListF5HOD_R24,optParamsF5HOD_R24=catF5HOD_R24.calcConfidenceEllipses(startIndex=startIndex,endIndex=endIndex,rMin=rMin,rMax=rMax,errors=errors,load=load,XYZ=XYZ,shiftFac=shiftFac,save=save,tol=tol,loadHessian=loadHessian)
    betaF5HOD_R24=optimizeBetaFixedShift(catF5HOD_R24,rMin=rMin,rMax=rMax,plotting=False,shiftFac=shiftFac)
    dAListF5HOD_R24_95=catF5HOD_R24.dAList_95


    "N1HOD"
    dAListN1HOD_R24,optParamsN1HOD_R24=catN1HOD_R24.calcConfidenceEllipses(startIndex=startIndex,endIndex=endIndex,rMin=rMin,rMax=rMax,errors=errors,load=load,XYZ=XYZ,shiftFac=shiftFac,save=save,tol=tol,loadHessian=loadHessian)
    betaN1HOD_R24=optimizeBetaFixedShift(catN1HOD_R24,rMin=rMin,rMax=rMax,plotting=False,shiftFac=shiftFac)
    dAListN1HOD_R24_95=catN1HOD_R24.dAList_95

  
    numCols=1
    numRows=1
    scale=1.6
    

    
    "Plotting" 
    axs.plot(dAListGRHOD_R24[0],dAListGRHOD_R24[1],color=colors_dic['blue'])
    axs.plot([optParamsGRHOD_R24[0]] ,[optParamsGRHOD_R24[1]],marker="*",color=colors_dic['blue'])
    xVals=dAListGRHOD_R24[0]
    yVals=dAListGRHOD_R24[1]
    axs.plot([betaGRHOD_R24[0],betaGRHOD_R24[0]],ylimits,color=colors_dic['blue'],linestyle="--")
        
            
    axs.plot(dAListF5HOD_R24[0],dAListF5HOD_R24[1],color=colors_dic['orange'])
    axs.plot([optParamsF5HOD_R24[0]] ,[optParamsF5HOD_R24[1]],marker="*",color=colors_dic['orange'])
    xVals=dAListF5HOD_R24[0]
    yVals=dAListF5HOD_R24[1]
    axs.plot([betaF5HOD_R24[0],betaF5HOD_R24[0]],ylimits,color=colors_dic['orange'],linestyle="--")
        
    
            
    axs.plot(dAListN1HOD_R24[0],dAListN1HOD_R24[1],color=colors_dic['green'])
    axs.plot([optParamsN1HOD_R24[0]] ,[optParamsN1HOD_R24[1]],marker="*",color=colors_dic['green'])
    xVals=dAListN1HOD_R24[0]
    yVals=dAListN1HOD_R24[1]
    axs.plot([betaN1HOD_R24[0],betaN1HOD_R24[0]],ylimits,color=colors_dic['green'],linestyle="--")
        
        
    axs.set_xlabel(r"$\beta$",fontsize=12)
    axs.set_ylabel(r"$\sigma_{0}$ $[km/s]$",fontsize=12)

    
    axs.set_title(r"$R_{eff} \geq 35$ Mpc/h",fontsize=12)
    axs.legend([r"GR - $(\beta,\sigma_{0})$",r"F5 - $(\beta,\sigma_{0})$",r"N1 - $(\beta,\sigma_{0})$"],fontsize=11.5)
    
    axs.set_ylim([240,340])
        
    
    lines = axs.get_lines()
    legend1 = axs.legend([lines[i] for i in [1,4,7]], [r"GR - $(\beta,\sigma_{0})$",r"F5 - $(\beta,\sigma_{0})$",r"N1 - $(\beta,\sigma_{0})$"], loc=4,fontsize=10)
    legend2 = axs.legend([lines[i] for i in [2,5,8]], [r"GR - $\beta$ from $v_{r}$",r"F5 - $\beta$ from $v_{r}$",r"N1 - $\beta$ from $v_{r}$"], loc=2,fontsize=10)
    axs.add_artist(legend1)
    axs.add_artist(legend2)
            

    saveName="Fig_7_Final_cb"
    
    # plt.savefig(fname= plotDir + saveName + ".pdf",bbox_inches="tight")

    # plt.show()
    # plt.clf()
   

   








