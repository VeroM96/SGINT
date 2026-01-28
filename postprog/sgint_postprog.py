# %% 
import numpy as np
import matplotlib.pyplot as plt
import pickle 
import os 
import datetime
import yaml
from scipy.interpolate import CubicSpline

plotting = False
## Based on Zharova et al. 2001
##Changelog:
# mean_v1: Initial version, based on v11
#
#mean_v2: changes to dessication funciton
#mean_v4: changes to growth function
#         -no loss of leaf biomass due to aging
#         -dessication function changed based on temperature function in Voinov and Akhremenkov 1990, and moved from limiting growth to mortality
#mean_v5: no aging for single shoot biomass mortality
#mean_v6: leaf sloughing due to wind added

#STATEPARAMETERS:
#P : single leaf biomass [gc]
#N : shoot density [shoots/m²]
#R : root biomass [gC/m²]

def sgint(dt,S,N,R,temp,tdry,I_can,wind,doy, runpath,debug=False):
    # Load parameters from yaml file
    with open(f'{runpath}/sgint.yaml') as info:
        params = yaml.full_load(info)
    # Explicitly assign parameters from the YAML file to local variables
    grow_s = params['grow_s']            # maximum relative growth rate [d^-1]
    I_k20 = params['I_k20']              # saturation irradiance [W/m²]
    I_c20 = params['I_c20']              # compenstation irradiance at 20 degC [W/m²]
    t_Ic = params['t_Ic']                # temperature coefficient for compensation irradiance
    t_Ik = params['t_Ik']                # temperature coefficient for saturation irradiance
    
    k0_ph = params['k0_ph']              # photosynthesis function value at T = 0 degC [d^-1]
    km_ph = params['km_ph']              # photosynthesis function value at T = T_m [d^-1]
    t_opt_ph = params['t_opt_ph']        # optimal temperature for photosynthetic growth [degC]
    t_max_ph = params['t_max_ph']        # maximum temperature for photosynthetic growth [degC]
    stt_ph = params['stt_ph']            # controlled shape of function for relative growth rate
    loss_s = params['loss_s']            # leaf mortality [d^-1]
    t_loss = params['t_loss']            # temperature coefficient for leaf mortality
    
    tdrymin = params['tdrymin']          # minimum time for dry period [h] per 12 hours <= 12!
    tdrymax = params['tdrymax']          # maximum time for dry period [h] per 12 hours <= 12!
    tdrya = params['tdrya']              # rate coefficient for dry period t<tdrymin
    tdryb = params['tdryb']              # rate coefficient for dry period t>tdrymax
    sig_s = params['sig_s']              # maximum single leaf biomass
    
    grow_n = params['grow_n']            # maximum relative recruitment rate [d^-1]
    k0_pr = params['k0_pr']              # function value at T = 0 degC [d^-1]
    km_pr = params['km_pr']              # function value at T = T_m [d^-1]
    t_opt_pr = params['t_opt_pr']        # optimal temperature for recruitment [degC]
    t_max_pr = params['t_max_pr']        # maximum temperature for recruitment [degC]
    stt_pr = params['stt_pr']            # controlled shape of function for relative recruitment rate
    sig_ag = params['sig_ag']            # above ground biomass above which no new shoots are produced [gC m^-2]
    eps = params['eps']                  # half-saturated constant for below ground biomass [gC m^-2]
    loss_n = params['loss_n']            # shoot loss rate [d^-1]
    
    d_loss = params['d_loss']            # day of year at which aging of leaves increses shoot loss rate [d]
    d_grow = params['d_grow']            # day of year at which recruitment starts [d] 
    d_loss_r = params['d_loss_r']        # steepness of curve for loss rate of below ground biomass [d] (how fast loss rate increases)
    lmrv = params['lmrv']                # leaf mass removal via sloughing [d^-1 per (m/s) wind speed]
    
    k_rs = params['k_rs']                # fraction of growth that is allocated to below ground growth k_R:S
    r_upr = params['r_upr']              # below ground biomass of plant lost due to uprooting [gC]
    loss_r = params['loss_r']            # loss rate of below ground biomass [d^-1]
    
    AG = S*N
    #light limitation
    I_c = I_c20*(t_Ic**(temp-20))
    I_k = I_k20*(t_Ik**(temp-20))

    if I_can <= I_c:
        f_il = 0
    elif I_can >= I_k:
        f_il = 1
    else:
        f_il = (I_can-I_c)/(I_k-I_c)

    #temperature dependency
    if temp <= t_opt_ph:
        f_tphot = k0_ph**((t_opt_ph-temp)/t_opt_ph)**stt_ph
    elif temp > t_opt_ph:
        f_tphot = km_ph**((temp-t_opt_ph)/(t_max_ph-t_opt_ph))**stt_ph
 

    #maximum canopy height
    if S <= sig_s:
        f_slim = (1-(S/sig_s)**2) 
    else:
        f_slim = 0
    
    #---- biomass mortality ----------------------------------------------
    #temperature dependeny
    f_tloss = max(t_loss**(temp-20),0)

    #----shoot recruitment -----------------------------------------------
    #space limitation above ground biomass
    if AG <= sig_ag:
        f_aglim = 1-((AG/sig_ag)**2) 
    else:
        f_aglim = 0
    
    #space limitation below ground biomass
    f_rlim = max((R/(R+eps)),0)

    #temperature dependency
    if temp <= t_opt_pr:
        f_tprod = k0_pr**(((t_opt_pr-temp)/t_opt_pr)**stt_pr)
    elif temp > t_opt_pr:
        f_tprod = km_pr**((temp-t_opt_pr)/(t_max_pr-t_opt_pr))**stt_pr

    #---- aging ---------------------------------------------------------
    f_age = (np.tanh((((doy-d_grow)%366)-(d_loss-d_grow))/d_loss_r) +1 )/2

    #---- innundation ---------------------------------------------------

    if tdry < tdrymin:
        f_tdry = 1 + tdrya * (tdrymin-tdry)**2
    elif tdry > tdrymax:
        f_tdry = 1 + tdryb * (tdry-tdrymax)**2
    else:
        f_tdry = 1

    #---- leaf sloughing ------------------------------------------------
    f_wind = 1+(lmrv * wind)
        
    #---- root biomass growth -------------------------------------------
    R1 = (k_rs *grow_s *f_il *f_tphot * AG)
    R3 = (loss_n*r_upr*f_age*f_tdry*N) 
    R2 = (loss_r *R)

    
    #---- solve differential equations------------------------------------
    S_gain = (grow_s * f_il *f_tphot *f_slim) * dt *S
    S_loss = (loss_s*f_tloss*f_age*f_tdry) * dt *S 
    N_gain = (grow_n * f_il * f_tprod *f_aglim *f_rlim) * dt *N
    N_loss = (loss_n*f_age*f_tdry*f_wind) * dt *N
    R_gain = R1*dt
    R_loss = (R2 + R3 )*dt

    S += S * ((grow_s * f_il *f_tphot *f_slim) - (loss_s*f_tloss*f_age*f_tdry))  *dt
    N += N * ((grow_n * f_il * f_tprod *f_aglim *f_rlim) - (loss_n*f_age*f_tdry*f_wind))  *dt
    R += (R1 - R2 - R3)  *dt
    
    if R <= 0:
        R = 0
        S = 0
        N = 0
    if S <= 0:
        S = 0
        N = 0
    if N <= 0:
        N = 0
        S = 0
    if debug:
        
        return S,N,R, S_gain,S_loss,N_gain,N_loss,R_gain,R_loss
    else:
        return S,N,R