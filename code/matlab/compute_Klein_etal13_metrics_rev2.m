% _rev2: Revised by Mark Zelinka on 23 June 2014:
% regridding using the interp function has been modified

% _rev: Revised by Mark Zelinka on 20 June 2014:
% It was discovered that the calculation of Eq. 3 of Klein et al. (2013)
% was being performed erroneously in the previous version. 
% Although the effect is negligible, this has been fixed here

% Original script written by Mark Zelinka (zelinka1@llnl.gov) on 7 October 2013

% This script computes the four Scalar Measures of the Fidelity of Model
% Cloud Simulations described in Section 4 of Klein et al. (2013)

% Reference: Klein, S.A., Y. Zhang, M.D. Zelinka, R.N. Pincus, J.Boyle, and P.J. Gleckler, 2013: 
% Are climate model simulations of clouds improving? An evaluation using the ISCCP simulator.  
% J. Geophys. Res. 118, 1329-1342. doi: 10.1002/jgrd.50141

% This script is written to perform the calculations for a singe model 
% (MPI-ESM-LR) and for MODIS observations

% Data that are used in this script:
% 1. observed ISCCP pctau histograms
% 2. observed MODIS pctau histograms (not essential)
% 3. model clisccp fields
% 4. model rsuscs fields
% 5. model rsdscs fields
% 6. cloud radiative kernels

% Other scripts that are called by this script:
% 1. get_netcdf_data.m
% 2. count_wt_mean.m

% Lines that include the comment "USER: MODIFY THIS LINE" are those in which 
% the user must update the path to the relevant netcdf files
   
%% Load in the cloud radiative kernels of Zelinka et al. (2012)
% Kernels are in units of W/m2/%
% Size: (mo, tau_midpt, p_midpt, lat, albcs) 
% The following are size (mo,tau,CTP,lat,albcs)
% where mo=1:12, tau, CTP, and albedo bins are given above

variable='/g/g19/zelinka1/CMIP3_kernels/cloud_kernels2.nc'; % USER: MODIFY THIS LINE 
kern_lat = get_netcdf_data(variable,'lat');  
coslat=cos(pi*kern_lat./180);
pm60=find(kern_lat>=-60 & kern_lat<60);
kern_lon=1.25:2.5:360;
LWkernel = get_netcdf_data(variable,'LWkernel'); 
SWkernel = get_netcdf_data(variable,'SWkernel'); 
albcs_midpt = get_netcdf_data(variable,'albcs'); % 0, 0.5, 1.0

%% Load model clisccp, rsuscs, rsutcs
variable='/p/lscratche/zelinka1/cmip5/clisccp/clisccp_cfMon_MPI-ESM-LR_amip_r1i1p1_197901-200812.nc'; % USER: MODIFY THIS LINE 
clisccp = get_netcdf_data(variable,'clisccp'); % size (time,tau,CTP,lat,lon); units: percent
clisccp(clisccp>500)=NaN; 

variable='/p/lscratche/zelinka1/cmip5/rsdscs/rsdscs_Amon_MPI-ESM-LR_amip_r1i1p1_197901-200812.nc'; % USER: MODIFY THIS LINE 
rsdscs = get_netcdf_data(variable,'rsdscs'); % size (time,lat,lon)
rsdscs(rsdscs>500)=NaN; 

variable='/p/lscratche/zelinka1/cmip5/rsuscs/rsuscs_Amon_MPI-ESM-LR_amip_r1i1p1_197901-200812.nc'; % USER: MODIFY THIS LINE 
lat = get_netcdf_data(variable,'lat'); % size (time,lat,lon)
lon = get_netcdf_data(variable,'lon'); % size (time,lat,lon)
rsuscs = get_netcdf_data(variable,'rsuscs'); % size (time,lat,lon)
rsuscs(rsuscs>500)=NaN; 

% Compute clear-sky surface albedo
albcs=rsuscs./rsdscs; 
clear rsuscs rsdscs

%% Compute anual cycles
[a,b,c,d,e]=size(clisccp);
avg_albcs=nan*ones(12,d,e);
avg_clisccp=nan*ones(12,7,7,d,e);
% before blindly doing this for every model, make sure your model data start in January
for i=1:12
    avg_clisccp(i,:,:,:,:)=squeeze(nanmean(clisccp(i:12:end,:,:,:,:),1));
    avg_albcs(i,:,:)=squeeze(nanmean(albcs(i:12:end,:,:),1));
end
clear clisccp albcs

%% Interpolate to a common grid (use the same grid as the kernels)
% first, a kluge to deal with weirdness of Matlab's interp function
x=[lon;lon(1)+360];
cat_albcs=cat(3,avg_albcs,avg_albcs(:,:,1)); 
cat_clisccp=cat(5,avg_clisccp,avg_clisccp(:,:,:,:,1)); 
clear avg_clisccp avgctl_albcs
    
[X1,Y1] = meshgrid(x,lat);
[X2,Y2] = meshgrid(kern_lon,kern_lat);
avg_clisccp_mod=nan*ones(12,7,7,90,144);
avg_albcs2=nan*ones(12,90,144);
for M=1:12
    avg_albcs2(M,:,:) = interp2(X1,Y1,squeeze(cat_albcs(M,:,:)),X2,Y2);  
    for T=1:7
        for P=1:7
            avg_clisccp_mod(M,T,P,:,:) = interp2(X1,Y1,squeeze(cat_clisccp(M,T,P,:,:)),X2,Y2);
        end
    end
end
clear cat_albcs cat_clisccp

%% load observed annual cycles
names(1)={'AC_clisccp'}; 
names(2)={'AC_clmodis'};
clisccp_obs=nan*ones(12,7,7,90,144,length(names));
for mm=1:length(names)
    variable=['/p/lscratche/zelinka1/cmip5/processed/yuying/',char(names(mm)),'.nc'];  % USER: MODIFY THIS LINE 
    avgclisccp = get_netcdf_data(variable,'avgclisccp');
    avgclisccp(avgclisccp>500)=NaN;
    clisccp_obs(:,:,:,:,:,mm)=100*avgclisccp; % units: percent
end
avg_clisccp_obs=flipdim(clisccp_obs,3); % SFC to TOA
clear avgclisccp clisccp_obs

% Note: both model and observed cloud fractions should be expressed in percent
% and both should be size (time=12,tau=7,CTP=7,lat=90,lon=144); 

%% TOTAL CLOUD AMOUNT METRIC
cutoff=2; % only clouds with tau>1.3
obs_pctau=squeeze(avg_clisccp_obs(:,cutoff+1:end,:,:,:,1));
LT=size(obs_pctau,2);

% The old (incorrect) way of computing Eq. 3 in Klein et al. (2013):
%cbar=count_wt_mean(squeeze(nanmean(nanmean(obs_pctau(:,:,:,pm60,:),1),5)),repmat(reshape(coslat(pm60),1,1,60),[LT,7,1]),3);
%anomobs=obs_pctau-repmat(reshape(cbar,1,LT,7,1,1),[12,1,1,90,144]); % anomaly of obs from its spatio-temporal mean
%ctot=squeeze(sum(sum(anomobs(:,:,:,pm60,:),2),3)); % size(12,Lreg,144)
%Z1denom=sqrt(count_wt_mean(squeeze(nanmean(nanmean(ctot.^2,1),3)),coslat(pm60)',2));

% The new (correct) way of computing Eq. 3 in Klein et al. (2013):
cbar=count_wt_mean(squeeze(nanmean(nanmean(sum(sum(obs_pctau(:,:,:,pm60,:),2),3),1),5)),coslat(pm60),1); % scalar
anomobs=squeeze(sum(sum(obs_pctau(:,:,:,pm60,:),2),3))-repmat(reshape(cbar,1,1,1),[12,60,144]); % anomaly of obs from its spatio-temporal mean
ctot=anomobs; % size(12,Lreg,144)
Z1denom=sqrt(count_wt_mean(squeeze(nanmean(nanmean(ctot.^2,1),3)),coslat(pm60)',2));

for dd=1:2 % Loop through possible differences       
    % Compute differences with respect to ISCCP
    if dd==2 % MODIS minus ISCCP
        anompctau=squeeze(nanmean(avg_clisccp_obs(:,cutoff+1:end,:,:,:,2),6))-obs_pctau;
    else
        anompctau=avg_clisccp_mod(:,cutoff+1:end,:,:,:)-obs_pctau;
    end    
    ctot=squeeze(sum(sum(anompctau(:,:,:,pm60,:),2),3)); % size(12,Lreg,144)
    Z1(dd)=sqrt(count_wt_mean(squeeze(nanmean(nanmean(ctot.^2,1),3)),coslat(pm60)',2));
end
     
%% CLOUD PROPERTY METRICS
cutoff=3; % only clouds with tau>3.6
obs_pctau=squeeze(avg_clisccp_obs(:,cutoff+1:end,:,:,:,1));
avg_clisccp_mod(:,1:cutoff,:,:,:,:)=[];
[a,b,c,d,e]=size(avg_clisccp_mod);
anompctau=avg_clisccp_mod-obs_pctau;
cbar=count_wt_mean(squeeze(nanmean(nanmean(obs_pctau(:,:,:,pm60,:),1),5)),repmat(reshape(coslat(pm60),1,1,60),[b,7,1]),3);
anomobs=obs_pctau-repmat(reshape(cbar,1,b,7,1,1),[12,1,1,90,144]); % anomaly of obs from its spatio-temporal mean
clear obs_pctau

% Reshape SW kernel to be the same size as pctau by mapping the albedos to the respective longitudes
SWkernel(:,1:cutoff,:,:,:)=[];
LWkernel(:,1:cutoff,:,:,:)=[];
SWkernelmod=NaN*ones(12,b,c,d,e);
for mm=1:12
    for i=1:90
        alon=squeeze(avg_albcs2(mm,i,:)); 
        if numel(alon(isnan(alon)))>0; continue; end
        for pp=1:7
            SWkernelmod(mm,:,pp,i,:)=interp2(albcs_midpt,1:b, squeeze(SWkernel(mm,:,pp,i,:)),alon,1:b);
        end
    end
end
clear alon avg_albcs2
LWkernelmod=repmat(LWkernel(:,:,:,:,1),[1,1,1,1,144]);

%% Multiply cloud fraction anomalies with kernels
mySW=anompctau.*SWkernelmod;
myLW=anompctau.*LWkernelmod;
obsSW=anomobs.*SWkernelmod;
obsLW=anomobs.*LWkernelmod;
clear LWkernelmod SWkernelmod

%% Aggregate high, mid, and low clouds over medium and thick ISCCP ranges
lo=1:2; mid=3:4; hi=5:7;
anompctau2=nan*ones(12,2,3,90,144); 
mySW2=nan*ones(12,2,3,90,144); 
myLW2=nan*ones(12,2,3,90,144);
anomobs2=nan*ones(12,2,3,90,144); 
obsSW2=nan*ones(12,2,3,90,144); 
obsLW2=nan*ones(12,2,3,90,144);
ttt=0;
for tt=[1,3]
    ttt=ttt+1;
    for pp=1:3
        if pp==1; PP=lo; elseif pp==2; PP=mid; elseif pp==3; PP=hi; end
        anompctau2(:,ttt,pp,:,:)=sum(sum(anompctau(:,tt:tt+1,PP,:,:),2),3);
        mySW2(:,ttt,pp,:,:)=sum(sum(mySW(:,tt:tt+1,PP,:,:),2),3);
        myLW2(:,ttt,pp,:,:)=sum(sum(myLW(:,tt:tt+1,PP,:,:),2),3);
        anomobs2(:,ttt,pp,:,:)=sum(sum(anomobs(:,tt:tt+1,PP,:,:),2),3);
        obsSW2(:,ttt,pp,:,:)=sum(sum(obsSW(:,tt:tt+1,PP,:,:),2),3);
        obsLW2(:,ttt,pp,:,:)=sum(sum(obsLW(:,tt:tt+1,PP,:,:),2),3);
    end
end

%% Compute sums and global averages
ctot1=squeeze(sum(sum(anompctau2(:,:,:,pm60,:).^2,2),3))/6;
Z2=sqrt(count_wt_mean(squeeze(nanmean(nanmean(ctot1,1),3)),coslat(pm60)',2));
ctot2=squeeze(sum(sum(anomobs2(:,:,:,pm60,:).^2,2),3))/6;
Z2denom=sqrt(count_wt_mean(squeeze(nanmean(nanmean(ctot2,1),3)),coslat(pm60)',2));

sumSW1=squeeze(sum(sum(mySW2(:,:,:,pm60,:).^2,2),3))/6;
Z3SW=sqrt(count_wt_mean(squeeze(nanmean(nanmean(sumSW1,1),3)),coslat(pm60)',2));
sumLW1=squeeze(sum(sum(myLW2(:,:,:,pm60,:).^2,2),3))/6;
Z3LW=sqrt(count_wt_mean(squeeze(nanmean(nanmean(sumLW1,1),3)),coslat(pm60)',2));
sumSW2=squeeze(sum(sum(obsSW2(:,:,:,pm60,:).^2,2),3))/6;
Z3SWdenom=sqrt(count_wt_mean(squeeze(nanmean(nanmean(sumSW2,1),3)),coslat(pm60)',2));
sumLW2=squeeze(sum(sum(obsLW2(:,:,:,pm60,:).^2,2),3))/6;
Z3LWdenom=sqrt(count_wt_mean(squeeze(nanmean(nanmean(sumLW2,1),3)),coslat(pm60)',2));
clear anompctau anompctau2 anomctp anomtau anomobs mySW myLW obsSW obsLW

%% CLOUD METRICS
% 1. MPI minus ISCCP 
% 2. MODIS minus ISCCP (for E_TCA only)

E_TCA=Z1/Z1denom % total cloud amount error 
E_CP=Z2/Z2denom % Cloud properties error 
E_LWCP=Z3LW/Z3LWdenom % LW-relevant Cloud properties error 
E_SWCP=Z3SW/Z3SWdenom % SW-relevant Cloud properties error 