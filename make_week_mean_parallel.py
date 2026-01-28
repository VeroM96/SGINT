import os
from concurrent.futures import ProcessPoolExecutor
import sys
sys.path.insert(0, '/work/gg0877/g260204/tools/python_skripts/SchismUtils/')
from schism_utils import  read_data 
import numpy as np
import h5py
import glob
import yaml

#V. Mohr, 18.12.2025
#make weekly mean files from SCHISM output (unmerged, old ouput format)

#### INPUTS ############################
#check if arguments are provided, else throw an error and exit
if len(sys.argv) < 2:
    sys.exit("Error run_week_mean_parallel.py: Missing required argument 'runpat'. Usage: python make_week_mean_parallel.py <runpat>")

#read run pattern from command line argument
runpat = sys.argv[1]
ts = sys.argv[2] 

simpat = f'{runpat}/outputs'
if not os.path.exists(simpat):
    sys.exit(f"Error run_week_mean_parallel.py: The specified path '{simpat}' does not exist.")

#check number of files
nfiles = len(glob.glob(f'{simpat}/local_to_global*'))

outdir = f'{runpat}/outputs/week_mean'
if not os.path.exists(outdir):
    os.makedirs(outdir)
    print(f"Created output directory: {outdir}")

#### PARAMETERS ############################
with open(f'{runpat}/sgint.yaml') as info:
    params = yaml.full_load(info)
k_w = params['k_w'] #PAR attenuation coefficient [m-1] 
k_sed = params['k_sed'] #sediment attenuation coefficient [m-1]


#### FUNCTIONS #############################
def process_file(fnum):
    #check if file exists
    fname = f'{simpat}/schout_00{fnum:04d}_{ts}.nc'
    if not os.path.exists(fname):
        sys.exit(f"Error  run_week_mean_parallel.py: File {fname} does not exist")
    #read data from schism output file
    times, temp, wind, dry, rad, zcor, tsc = read_data(['time','temp','wind_speed','wetdry_node','solar_radiation','zcor','SED_TSC'], [fname])
    print(f'Processing file {fname} with {len(times)} time steps')

    #calculate mean over weekly period, limited to m2 tide cycles
    m2 = 12.42   #duration m2 tide [hours]
    dt = times[1] - times[0] #duration of SCHISM time step in seconds
    m2_seconds = m2 * 3600
    total_duration = times[-1] - times[0] #total duration of SCHISM output in seconds
    n_m2 = int(total_duration // m2_seconds) # number of 12.42 hour periods in the total duration
    n_steps = int(m2_seconds*n_m2 // dt) #number of steps for m2 periods in the total duration

    #take last n_steps for mean calculation
    mean_temp = np.mean(np.mean(temp[-n_steps:,:,:], axis=0), axis=1)
    mean_dry = np.sum(dry[-n_steps:]*dt, axis=0)/n_m2 /3600 # convert to hours

    #calculate light attenuation and store mean light intensity for lowest layer
    #light attenuation 

    dz = zcor[:,:,1:]-zcor[:,:,:-1]
    if dz.any()< 0:
        print(dz)
    i_n = rad #initial light at surface
    for i in range(dz.shape[2]-1,-1,-1):
        kw = k_w+k_sed*tsc[:,:,i] #total attenuation coefficient
        i_wc = i_n *np.exp(-kw*dz[:,:,i]/2)#calculate light at center of layer   
        i_n = i_wc * np.exp(-kw*dz[:,:,i]) #calculate light at bottom of layer
    mean_rad = np.mean(i_wc[-n_steps:,:], axis=0)
    #mean depth
    mean_depth = np.mean(np.sum(dz[-n_steps:,:], axis=0), axis=1)  # average depth over the last n_steps
    #wind extinction with depth
    k_wind = 1.2 #wind effect attenuation with depth
    wind_speed = np.sqrt(wind[:,:,1]**2+wind[::,:,0]**2)
    if np.any(zcor[:,:,-1]-zcor[:,:,1] < 0):
        print(zcor)
    wind_bot = wind_speed/10 * np.exp(-k_wind* (zcor[:,:,-1]-zcor[:,:,1]))#wind effect at bottom
    mean_wind = np.mean(wind_bot[-n_steps:,:], axis=0)


    with h5py.File(f'{outdir}/mean_{fnum:06d}_{ts}.h5', 'w') as f:
        f.create_dataset('dt', data=dt)
        f.create_dataset('n_steps', data=n_steps)
        f.create_dataset('mean_temp', data=mean_temp)
        f.create_dataset('mean_dry', data=mean_dry)
        f.create_dataset('mean_wind', data=mean_wind)
        f.create_dataset('mean_rad', data=mean_rad)
        f.create_dataset('mean_depth', data=mean_depth)
        print(f'Processed file {fname} and saved results to {outdir}/mean_{fnum:06d}_{ts}.h5')


#paralelizing
with ProcessPoolExecutor() as executor:
    executor.map(process_file, range(nfiles))
