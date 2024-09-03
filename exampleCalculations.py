#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 17:38:24 2024

@author: christopherwilson
"""

from voids_MG_GSM_methods import *

#%%
"Should now have loaded the catalogs objects"
"""
catGR_CC_R24_z0p5 - voids identified from halos in GR, includes real space void density, velocity, and velocity dispersion profiles but no redshift space quantities
catF5_CC_R24_z0p5 - same but for F5
catN1_CC_R24_z0p5 - same but for N1

catGRHOD_CC_R24_z0p5 - voids identified from HOD mock galaxies in GR. Indludes real space void density, velocity, and velocity dispersion profiles along with redshift space tracer locations to compute Xi_{0} and Xi_{2}. Uses only LOS=Z
catF5HOD_CC_R24_z0p5 - same but for F5
catN1HOD_CC_R24_z0p5 - same but for F6

catGRHODXYZ_CC_R24_z0p5 - same as above but also with LOS=X,Y, and Z independantly so that a mean Xi_{0,2} may be calculated more accurately. For the most accurate means, use these. For most accurate errors, use above 
catF5HODXYZ_CC_R24_z0p5 - same as above but for F5
catN1HODXYZ_CC_R24_z0p5 - same as above but for N1
"""



"Compute halo Void size function"
halos_nBins=26
halos_rMax=65
halos_rMin=0

GR_VSF_vol=catGR_CC_R24_z0p5.jackKnifeVSF_VolNorm(rMin=halos_rMin, rMax=halos_rMax, nBins=halos_nBins)
F5_VSF_vol=catF5_CC_R24_z0p5.jackKnifeVSF_VolNorm(rMin=halos_rMin, rMax=halos_rMax, nBins=halos_nBins)
N1_VSF_vol=catN1_CC_R24_z0p5.jackKnifeVSF_VolNorm(rMin=halos_rMin, rMax=halos_rMax, nBins=halos_nBins)


"Compute HOD Void size function"
hod_nBins=32
hod_rMax=80
hod_rMin=0

GRHOD_VSF_vol=catGRHOD_CC_R24_z0p5.jackKnifeVSF_VolNorm(rMin=hod_rMin, rMax=hod_rMax, nBins=hod_nBins)
F5HOD_VSF_vol=catF5HOD_CC_R24_z0p5.jackKnifeVSF_VolNorm(rMin=hod_rMin, rMax=hod_rMax, nBins=hod_nBins)
N1HOD_VSF_vol=catN1HOD_CC_R24_z0p5.jackKnifeVSF_VolNorm(rMin=hod_rMin, rMax=hod_rMax, nBins=hod_nBins)




"Compute densities, velocities, and velocity dispersions around void centers for small void and large void populations"
"Code below is for HOD catalogs only, but can be modified for halos by changing catalogs as above"

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




"calculate redshift space monopoles and quadrupoles - using only LOS=Z"
"this will give a less accurate mean quadrupole moment, but the error is representative of the standard error on the mean for a total observational volume of 100 (Mpc/h)^3"
"can scale this error accordingly to represent different observational volumes"
quadData_GRHOD=catGRHODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=35,rMax=100)
quadData_F5HOD=catF5HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=35,rMax=100)
quadData_N1HOD=catN1HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=35,rMax=100)


quadData_GRHOD_S=catGRHODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=0,rMax=35)
quadData_F5HOD_S=catF5HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=0,rMax=35)
quadData_N1HOD_S=catN1HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=0,rMax=35)




"calculate redshift space monopoles and quadrupoles - combining all 3 lines of sight"
"since we have combined all 3 lines of sight, errors may be more correlated than usual and should not be trusted. Mean will be more accurate however"
quadData_GRHOD=catGRHODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=35,rMax=100)
quadData_F5HOD=catF5HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=35,rMax=100)
quadData_N1HOD=catN1HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=35,rMax=100)


quadData_GRHOD_S=catGRHODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=0,rMax=35)
quadData_F5HOD_S=catF5HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=0,rMax=35)
quadData_N1HOD_S=catN1HODXYZ_CC_R24_z0p5.jackKnifeQuadMoment(rMin=0,rMax=35)



"Calculate a predicted redshift space from the gaussian streaming model - can use either HOD or HODXYZ for this, and the answer shouldn't change"
quadThy_GRHOD_R24_L=catGRHOD_CC_R24_z0p5.calcFullGSM(rMin=35,rMax=100,numBinFac=10)
quadThy_F5HOD_R24_L=catF5HOD_CC_R24_z0p5.calcFullGSM(rMin=35,rMax=100,numBinFac=10)
quadThy_N1HOD_R24_L=catN1HOD_CC_R24_z0p5.calcFullGSM(rMin=35,rMax=100,numBinFac=10)

quadThy_GRHOD_R24_S=catGRHOD_CC_R24_z0p5.calcFullGSM(rMin=0,rMax=35,numBinFac=10)
quadThy_F5HOD_R24_S=catF5HOD_CC_R24_z0p5.calcFullGSM(rMin=0,rMax=35,numBinFac=10)
quadThy_N1HOD_R24_S=catN1HOD_CC_R24_z0p5.calcFullGSM(rMin=0,rMax=35,numBinFac=10)



"Calculate changes to Xi_{2} through functional derivatives of the GSM, with two different catalogs so that changes in \delta, v_{r}, or \sigma may be computed"
Spoints=(np.linspace(0,120,301)[0:-1] + np.linspace(0,120,301)[1:])/2
quadDelta_delta_F5HOD_R24_S=calcDeltaXi2_delta(catGRHOD_CC_R24_z0p5,catF5HOD_CC_R24_z0p5,spoints=Spoints,rMin=0,rMax=35,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"])
quadDelta_vel_F5HOD_R24_S=calcDeltaXi2_vel(catGRHOD_CC_R24_z0p5,catF5HOD_CC_R24_z0p5,spoints=Spoints,rMin=0,rMax=35,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"])
quadDelta_sigma_F5HOD_R24_S=calcDeltaXi2_sigma(catGRHOD_CC_R24_z0p5,catF5HOD_CC_R24_z0p5,spoints=Spoints,rMin=0,rMax=35,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"])
    
quadDelta_delta_N1HOD_R24_L=calcDeltaXi2_delta(catGRHOD_CC_R24_z0p5,catN1HOD_CC_R24_z0p5,spoints=Spoints,rMin=35,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"])
quadDelta_vel_N1HOD_R24_L=calcDeltaXi2_vel(catGRHOD_CC_R24_z0p5,catN1HOD_CC_R24_z0p5,spoints=Spoints,rMin=35,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"])
quadDelta_sigma_N1HOD_R24_L=calcDeltaXi2_sigma(catGRHOD_CC_R24_z0p5,catN1HOD_CC_R24_z0p5,spoints=Spoints,rMin=35,rMax=100,a=a,H=H,algorithmInt=True,numMuPoints=51,voidTypes=["R","S"])















