
#!/usr/bin/env cdat
"""
This script computes the four Scalar Measures of the Fidelity of Model
Cloud Simulations described in Section 4 of Klein et al. (2013) for an example CMIP5 model

Returns:
E_TCA:    total cloud amount error 
E_CP:     Cloud properties error 
E_LWCP:   LW-relevant Cloud properties error 
E_SWCP:   SW-relevant Cloud properties error 

Reference: Klein, S.A., Y. Zhang, M.D. Zelinka, R.N. Pincus, J.Boyle, and P.J. Gleckler, 2013: 
Are climate model simulations of clouds improving? An evaluation using the ISCCP simulator.  
J. Geophys. Res. 118, 1329-1342. doi: 10.1002/jgrd.50141
"""

#IMPORT STUFF:
#=====================
import cdms2 as cdms
import cdutil
import MV2 as MV
import numpy as np
import pylab as pl

###########################################################################
# HELPFUL FUNCTIONS FOLLOW 
###########################################################################

###########################################################################
def add_cyclic(data):
    # Add Cyclic point around 360 degrees longitude:

    # This function assumes that your longitudes range from 0 to 360, not -180 to 180
    lons=data.getLongitude()[:]
    dx=np.gradient(lons)[-1]
    data2 = data(longitude=(0, dx+np.max(lons)), squeeze=True)    
    return data2

###########################################################################
def reshape_generic(orig_data,data_to_match):

    # this function resizes and tiles orig_data the same shape as data_to_match
    # orig_data must have fewer dimensions than data_to_match

    A=orig_data.shape
    B=data_to_match.shape
    ndim_new = data_to_match.ndim

    # find index where B disagrees with A
    #disagree=np.setdiff1d(B,A)
    agree=np.in1d(B,A)
    j=[]
    for i in range(len(B)):
        ndim_orig = orig_data.ndim
        if agree[i]==False:
            j=np.append(j,i)
            new = np.expand_dims(orig_data,axis=ndim_orig)
            NEW = np.tile(new,[B[i]])
            try:
                new_mask = np.expand_dims(orig_data.mask,axis=ndim_orig)
                MASK = np.tile(new_mask,[B[i]])
                orig_data = np.ma.array(NEW,mask=MASK)
            except:
                orig_data = np.ma.array(NEW)

    # need to move axes around
    for i in range(len(B)):
        C=orig_data.shape
        if C[i]!=B[i]:
            orig_data = np.moveaxis(orig_data, i, B.index(C[i])) # (a, source, destination)

    return orig_data

###########################################################################
def nanarray(vector):

    # this generates a masked array with the size given by vector
    # example: vector = (90,144,28)

    # similar to this=NaN*ones(x,y,z) in matlab

    this=MV.zeros(vector)
    this=MV.masked_where(this==0,this)

    return this

###########################################################################
def map_SWkern_to_lon(Ksw,albcsmap):

    from scipy.interpolate import interp1d
    ## Map each location's clear-sky surface albedo to the correct albedo bin
    # Ksw is size 12,7,7,lats,3
    # albcsmap is size A,lats,lons
    albcs=np.arange(0.0,1.5,0.5) 
    A=albcsmap.shape[0]
    TT=Ksw.shape[1]
    PP=Ksw.shape[2]
    lenlat=Ksw.shape[3]
    lenlon=albcsmap.shape[2]
    SWkernel_map=nanarray((A,TT,PP,lenlat,lenlon))
    for M in range(A):
        MM=M
        while MM>11:
            MM=MM-12
        for LA in range(lenlat):
            alon=albcsmap[M,LA,:] 
            # interp1d can't handle mask but it can deal with NaN (?)
            try:
                alon2=MV.where(alon.mask,np.nan,alon)   
            except:
                alon2=alon
            if np.ma.count(alon2)>1: # at least 1 unmasked value
                if len(pl.find(Ksw[MM,:,:,LA,:]>0))==0:
                    SWkernel_map[M,:,:,LA,:] = 0
                else:
                    f = interp1d(albcs,Ksw[MM,:,:,LA,:],axis=2)
                    ynew = f(alon2.data)
                    ynew=MV.masked_where(alon2.mask,ynew)
                    SWkernel_map[M,:,:,LA,:] = ynew
            else:
                continue

    return SWkernel_map

###########################################################################
# MAIN ROUTINE FOLLOWS
###########################################################################

datadir = '/work/zelinka1/git/klein2013-cloud-error-metrics/data/'

# Load in the Zelinka et al 2012 kernels:
f=cdms.open(datadir+'cloud_kernels2.nc')
LWkernel0=f('LWkernel')
SWkernel0=f('SWkernel')
f.close()

# Take only the portion of the kernel histogram where there are obs (ignore first tau bin)
SWkernel = SWkernel0[:,1:,:]
LWkernel = LWkernel0[:,1:,:]
del(LWkernel0,SWkernel0)

albcs=np.arange(0.0,1.5,0.5) # the clear-sky albedos over which the kernel is computed

######################################################
############# Load in ISCCP observations #############
######################################################
f=cdms.open(datadir+'AC_clisccp.nc','r')
obs_clisccp=f('avgclisccp',squeeze=1)   
f.close()

grid = obs_clisccp.getGrid()

# Flip the CTP dimension to go SFC to TOA, set to units of %, and ignore the 1st TAU bin:
obs_clisccp_grd = 100*obs_clisccp[:,1:,-1::-1,:]

######################################################
############# Load in MODIS observations #############
######################################################
f=cdms.open(datadir+'AC_clmodis.nc') 
obs_clmodis=f('avgclisccp',squeeze=1)  # not a typo
f.close()

# Flip the CTP dimension to go SFC to TOA, set to units of %, and ignore the 1st TAU bin:
obs_clmodis_grd = 100*obs_clmodis[:,1:,-1::-1,:]

agg_mod_clisccp_bias = nanarray((12,2,3,60,144)) # (month, tau_bins, CTP_bins, lat, lon)
agg_mod_SW_bias = nanarray((12,2,3,60,144)) # (month, tau_bins, CTP_bins, lat, lon)
agg_mod_LW_bias = nanarray((12,2,3,60,144)) # (month, tau_bins, CTP_bins, lat, lon)

agg_obs_clisccp_bias=nanarray((12,2,3,60,144))
agg_obs_SW_bias=nanarray((12,2,3,60,144))
agg_obs_LW_bias=nanarray((12,2,3,60,144))

######################################################
############# Load in CLISCCP from model #############
######################################################
# Grab a random AMIP simulation:
f=cdms.open(datadir+'clisccp_cfMon_MPI-ESM-LR_amip_r1i1p1_197901-200812.nc','r')
clisccp0=f('clisccp')
f.close()

# Compute Climatological Annual Cycle:
clisccp = cdutil.ANNUALCYCLE.climatology(clisccp0) #(12,...)
del(clisccp0)

# Remove the thinnest optical depth bin so as to compare properly with obs:
clisccp=clisccp[:,1:,:,:]

# Make sure cloud fractions are in percent  
sumclisccp=MV.sum(MV.sum(clisccp,2),1)
if np.max(sumclisccp) <= 1.:
    clisccp = clisccp*100.          

######################################################
########## Compute clear-sky surface albedo ##########
######################################################
f=cdms.open(datadir+'rsdscs_Amon_MPI-ESM-LR_amip_r1i1p1_197901-200812.nc','r')
rsdscs0 = f('rsdscs',squeeze=1) # Clearsky downwelling solar flux at surface
f=cdms.open(datadir+'rsuscs_Amon_MPI-ESM-LR_amip_r1i1p1_197901-200812.nc','r')
rsuscs0 = f('rsuscs',squeeze=1) # Clearsky upwelling solar flux at surface
f.close()

# Compute Climatological Annual Cycle:
rsdscs = cdutil.ANNUALCYCLE.climatology(rsdscs0) #(12,...)
rsuscs = cdutil.ANNUALCYCLE.climatology(rsuscs0) #(12,...)

albcs = rsuscs/rsdscs
albcs=MV.where(albcs>1.,1,albcs) # where(condition, x, y) is x where condition is true, y otherwise
albcs=MV.where(albcs<0.,0,albcs)

 
# Regrid everything to the kernel grid:
albcs = add_cyclic(albcs)
albcs_grd = albcs.regrid(grid,regridTool="esmf",regridMethod = "linear")
clisccp = add_cyclic(clisccp)
clisccp_grd = clisccp.regrid(grid,regridTool="esmf",regridMethod = "linear")

## Use average control albcs to map SW kernel to appropriate longitudes
SWkernel_map = map_SWkern_to_lon(SWkernel,albcs_grd)

# LW kernel does not depend on albcs, just repeat the final dimension over longitudes:
A=SWkernel_map.shape[0]
LWkernel_map0=np.tile(np.tile(LWkernel[:,:,:,:,0],(1,1,1,1,1)),(144,1,1,1,1))(order=[1,2,3,4,0])
LWkernel_map=nanarray(SWkernel_map.shape)
for a in range(A):
    aa=a
    while aa>11:
        aa=aa-12
    LWkernel_map[a,:] = LWkernel_map0[aa,:]

## Compute Cloud Fraction Histogram Anomalies w.r.t. observations
clisccp_bias = clisccp_grd - obs_clisccp_grd

## Multiply Anomalies by Kernels
SW0 = SWkernel_map*clisccp_bias
LW = LWkernel_map*clisccp_bias

## Set the SW cloud feedbacks to zero in the polar night
# The sun is down if every bin of the SW kernel is zero:
sundown=MV.sum(MV.sum(SWkernel_map,axis=2),axis=1)  #MO,90,144
repsundown=np.tile(np.tile(sundown,(1,1,1,1,1)),(7,6,1,1,1))(order=[2,1,0,3,4])
SW1 = MV.where(repsundown==0, 0, SW0) # where(condition, x, y) is x where condition is true, y otherwise
SW = MV.where(repsundown.mask, 0, SW1) # where(condition, x, y) is x where condition is true, y otherwise

# SW and LW contain the SW and LW radiation anomalies contributed from cloud anomalies in each bin of the histogram
LW.setAxisList(clisccp_bias.getAxisList())
SW.setAxisList(clisccp_bias.getAxisList())
LWkernel_map.setAxisList(clisccp_bias.getAxisList())
SWkernel_map.setAxisList(clisccp_bias.getAxisList())

########################################################
######### Compute Klein et al (2013) metrics ########### 
########################################################
eq60 = cdutil.region.domain(latitude=(-60.,60.)) # equatorward of 60

## E_TCA (TOTAL CLOUD AMOUNT METRIC)
# take only clouds with tau>1.3 between 60S-60N
obs_clisccp_eq60 = eq60.select(obs_clisccp_grd[:,1:,:]) 
obs_clmodis_eq60 = eq60.select(obs_clmodis_grd[:,1:,:]) 
mod_clisccp_eq60 = eq60.select(clisccp_grd[:,1:,:])  

# sum over CTP and TAU:
mod_cltisccp_eq60 = MV.sum(MV.sum(mod_clisccp_eq60,axis=1),axis=1) # (time, lat, lon)
obs_cltisccp_eq60 = MV.sum(MV.sum(obs_clisccp_eq60,axis=1),axis=1) # (time, lat, lon)
obs_cltisccp_eq60 = MV.sum(MV.sum(obs_clisccp_eq60,axis=1),axis=1) # (time, lat, lon)
obs_cltmodis_eq60 = MV.sum(MV.sum(obs_clmodis_eq60,axis=1),axis=1) # (time, lat, lon)

# Create map of area weights
lat=grid.getLatitude()[:]
coslats=np.cos(lat*np.pi/180)
w_k0=coslats/np.sum(np.cos(lat*np.pi/180))
w_k = w_k0[15:-15] # equatorward of 60

########################################################
# E_TCA for Model minus ISCCP:
########################################################
# 1) Denominator (Eq. 3 in Klein et al. (2013))
avg = cdutil.averager(MV.average(obs_cltisccp_eq60,axis=0), axis='xy', weights='weighted') # (scalar)
rep_avg = reshape_generic(avg,obs_cltisccp_eq60) # (time, lat, lon)
anom = obs_cltisccp_eq60 - rep_avg # anomaly of obs from its spatio-temporal mean
area_wts2 = reshape_generic(w_k,anom) # (time, lat, lon)
E_TCA_denom = np.ma.sqrt(MV.sum(area_wts2*anom**2)) # (scalar)
#E_TCA_denom = np.ma.sqrt(cdutil.averager(MV.average(anom**2,axis=0), axis='xy', weights='weighted')) # (scalar)

# 2) Numerator
anom = mod_cltisccp_eq60 - obs_cltisccp_eq60  # (time, lat, lon)
E_TCA_numer = np.ma.sqrt(MV.sum(area_wts2*anom**2)) # (scalar)
#E_TCA_numer = np.ma.sqrt(cdutil.averager(MV.average(anom**2,axis=0), axis='xy', weights='weighted')) # (scalar)
E_TCA_mod = E_TCA_numer/E_TCA_denom
     
# E_TCA for MODIS minus ISCCP (where they overlap):
# 2) Numerator
anom = obs_cltmodis_eq60 - obs_cltisccp_eq60  # (time, lat, lon)
E_TCA_numer = np.ma.sqrt(MV.sum(area_wts2*anom**2)) # (scalar)
#E_TCA_numer = np.ma.sqrt(cdutil.averager(MV.average(anom**2,axis=0), axis='xy', weights='weighted')) # (scalar)
E_TCA_obs = E_TCA_numer/E_TCA_denom

########################################################
# CLOUD PROPERTY METRICS
########################################################
# take only clouds with tau>3.6 between 60S-60N
clisccp_bias_eq60 = eq60.select(clisccp_bias[:,2:,:]) 
obs_clisccp_eq60 = eq60.select(obs_clisccp_grd[:,2:,:]) 
mod_clisccp_eq60 = eq60.select(clisccp_grd[:,2:,:])  
LWkernel_eq60 = eq60.select(LWkernel_map[:,2:,:])  
SWkernel_eq60 = eq60.select(SWkernel_map[:,2:,:])  

# Compute anomaly of obs histogram from its spatio-temporal mean
avg_obs_clisccp_eq60 = cdutil.averager(obs_clisccp_eq60, axis='xy', weights='weighted') # (time,TAU,CTP)
rep_avg_obs_clisccp_eq60 = reshape_generic(avg_obs_clisccp_eq60,obs_clisccp_eq60) # (time, TAU, CTP, lat, lon)
anom_obs_clisccp_eq60 = obs_clisccp_eq60 - rep_avg_obs_clisccp_eq60 # anomaly of obs from its spatio-temporal mean

## Compute radiative impacts of cloud fraction anomalies
mod_SW_bias = eq60.select(SW[:,2:,:])  
obs_SW_bias = anom_obs_clisccp_eq60*SWkernel_eq60
mod_LW_bias = eq60.select(LW[:,2:,:])  
obs_LW_bias = anom_obs_clisccp_eq60*LWkernel_eq60

## Aggregate high, mid, and low clouds over medium and thick ISCCP ranges
Psec_name = ['LO','MID','HI']
Psections=[slice(0,2),slice(2,4),slice(4,7)]
Psec_dic=dict(zip(Psec_name,Psections))
Tsec_name = ['MED','THICK']
Tsections=[slice(0,2),slice(2,4)]
Tsec_dic=dict(zip(Tsec_name,Tsections))

tt=-1
for Tsec in Tsec_name:
    tt+=1
    TT=Tsec_dic[Tsec]
    pp=-1
    for Psec in Psec_name:
        pp+=1
        PP=Psec_dic[Psec]
        agg_obs_SW_bias[:,tt,pp,:] = MV.sum(MV.sum(obs_SW_bias[:,TT,PP,:],axis=1),axis=1)
        agg_mod_SW_bias[:,tt,pp,:] = MV.sum(MV.sum(mod_SW_bias[:,TT,PP,:],axis=1),axis=1)
        agg_obs_LW_bias[:,tt,pp,:] = MV.sum(MV.sum(obs_LW_bias[:,TT,PP,:],axis=1),axis=1)
        agg_mod_LW_bias[:,tt,pp,:] = MV.sum(MV.sum(mod_LW_bias[:,TT,PP,:],axis=1),axis=1)
        agg_obs_clisccp_bias[:,tt,pp,:] = MV.sum(MV.sum(anom_obs_clisccp_eq60[:,TT,PP,:],axis=1),axis=1)
        agg_mod_clisccp_bias[:,tt,pp,:] = MV.sum(MV.sum(clisccp_bias_eq60[:,TT,PP,:],axis=1),axis=1)

## Compute E_ctp-tau -- Cloud properties error 
ctot = MV.sum(MV.sum(agg_mod_clisccp_bias**2,axis=1),axis=1)/6;
ctot.setAxisList(clisccp_bias_eq60[:,0,0,:].getAxisList())
E_ctpt_numer = np.ma.sqrt(MV.sum(area_wts2*ctot)) # (scalar)
#E_ctpt_numer = np.ma.sqrt(cdutil.averager(ctot, axis='xy', weights='weighted')) # (time)

ctot = MV.sum(MV.sum(agg_obs_clisccp_bias**2,axis=1),axis=1)/6;
ctot.setAxisList(clisccp_bias_eq60[:,0,0,:].getAxisList())
E_ctpt_denom = np.ma.sqrt(MV.sum(area_wts2*ctot)) # (scalar)
#E_ctpt_denom = np.ma.sqrt(cdutil.averager(ctot, axis='xy', weights='weighted')) # (time)

E_ctpt_mod = E_ctpt_numer/E_ctpt_denom

## Compute E_LW -- LW-relevant cloud properties error 
ctot = MV.sum(MV.sum(agg_mod_LW_bias**2,axis=1),axis=1)/6;
ctot.setAxisList(clisccp_bias_eq60[:,0,0,:].getAxisList())
E_LW_numer = np.ma.sqrt(MV.sum(area_wts2*ctot)) # (scalar)
#E_LW_numer = np.ma.sqrt(cdutil.averager(ctot, axis='xy', weights='weighted')) # (time)

ctot = MV.sum(MV.sum(agg_obs_LW_bias**2,axis=1),axis=1)/6;
ctot.setAxisList(clisccp_bias_eq60[:,0,0,:].getAxisList())
E_LW_denom = np.ma.sqrt(MV.sum(area_wts2*ctot)) # (scalar)
#E_LW_denom = np.ma.sqrt(cdutil.averager(ctot, axis='xy', weights='weighted')) # (time)

E_LW_mod = E_LW_numer/E_LW_denom

## Compute E_SW -- SW-relevant cloud properties error 
ctot = MV.sum(MV.sum(agg_mod_SW_bias**2,axis=1),axis=1)/6;
ctot.setAxisList(clisccp_bias_eq60[:,0,0,:].getAxisList())
E_SW_numer = np.ma.sqrt(MV.sum(area_wts2*ctot)) # (scalar)
#E_SW_numer = np.ma.sqrt(cdutil.averager(ctot, axis='xy', weights='weighted')) # (time)

ctot = MV.sum(MV.sum(agg_obs_SW_bias**2,axis=1),axis=1)/6;
ctot.setAxisList(clisccp_bias_eq60[:,0,0,:].getAxisList())
E_SW_denom = np.ma.sqrt(MV.sum(area_wts2*ctot)) # (scalar)
#E_SW_denom = np.ma.sqrt(cdutil.averager(ctot, axis='xy', weights='weighted')) # (time)

E_SW_mod = E_SW_numer/E_SW_denom

######## SANITY CHECK ######## 
print('E_TCA (model minus ISCCP): '+str(E_TCA_mod))
print('E_TCA (MODIS minus ISCCP): '+str(E_TCA_obs))
print('E_CTP_TAU: '+str(E_ctpt_mod))
print('E_LW: '+str(E_LW_mod))
print('E_SW: '+str(E_SW_mod))

#print('E_TCA (model minus ISCCP): '+str(np.ma.average(E_TCA_mod)))
#print('E_TCA (MODIS minus ISCCP): '+str(np.ma.average(E_TCA_obs)))
#print('E_CTP_TAU: '+str(np.ma.average(E_ctpt_mod)))
#print('E_LW: '+str(np.ma.average(E_LW_mod)))
#print('E_SW: '+str(np.ma.average(E_SW_mod)))

