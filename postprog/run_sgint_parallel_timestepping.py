import time
start = time.time()

import os
import yaml
import numpy as np
from sgint_postprog import sgint
import h5py
from concurrent.futures import ProcessPoolExecutor
import glob


run = 'sgtest014_2015_2'
outpat = f'/home/g/g260204/schism_v5.11/sgint_postprog/output/{run}/'
inpat = f'/home/g/g260204/schism_v5.11/test/{run}/'

plotting = False
debug = True

#check number of files
nfiles = len(glob.glob(f'{inpat}/outputs/local_to_global*'))
#### FUNCTIONS #############################
#write output to veg_ .gr3 files
def write_output_to_gr3(filename, updated_veg,nodes):
    with open(f'{outpat}/veg_N_0.gr3', 'r') as fin, open(filename, 'w') as fout:
        # Copy header
        fout.write(fin.readline())
        nElems, nNodes = map(int, fin.readline().split())
        fout.write(f"{nElems} {nNodes}\n")
        if len(nodes)==0:
            nodes = range(nNodes)
        
        for idx in range(len(nodes)):
            line = fin.readline()
            parts = line.split()
            # Replace the last value with updated_N
            fout.write(f"{parts[0]} {parts[1]} {parts[2]} {updated_veg[idx]}\n")

#read input veg_ .gr3 files
def read_input_veg(filename):
    with open(filename) as f:
        f.readline()  # skip title
        nElems, nNodes = map(int, f.readline().split())
        N = []
        for _ in range(nNodes):
            parts = f.readline().split()
            h = float(parts[3])
            N.append(h) 
    N = np.array(N)  # Convert to numpy array for easier manipulation
    return N
#### PARAMETERS ############################
with open(f'{inpat}/sgint.yaml') as info:
    params = yaml.full_load(info)
h0 = params['h0']               #minimum height of seagrass canopy [m]
a_ch = params['a_ch']           #conversion factor from gC to height [m gC^-1]
sgint_dt = params['sgint_dt']   #time step for sgint model [d]

#### OUTPUTS ############################
if not os.path.exists(outpat):
    os.makedirs(outpat)

#copy input veg_N.gr3 to output folder

os.system(f'cp {inpat}/veg_N_0.gr3 {outpat}/veg_N_0.gr3')
os.system(f'cp {inpat}/veg_RB_0.gr3 {outpat}/veg_RB_0.gr3')
os.system(f'cp {inpat}/veg_LB_0.gr3 {outpat}/veg_LB_0.gr3')
os.system(f'cp {inpat}/veg_h_0.gr3 {outpat}/veg_h_0.gr3')
#%%
step = 1
for step in range(1, 53):
    #### INITIALIZATION ############################  

    LB = read_input_veg(f'{outpat}/veg_LB_{step-1}.gr3')
    RB = read_input_veg(f'{outpat}/veg_RB_{step-1}.gr3')
    N = read_input_veg(f'{outpat}/veg_N_{step-1}.gr3')

    updated_LB = np.zeros_like(LB)
    updated_RB = np.zeros_like(RB)
    updated_N = np.zeros_like(N)
    updated_RB_gain = np.zeros_like(RB)
    updated_RB_loss = np.zeros_like(RB)
    updated_LB_gain = np.zeros_like(LB)
    updated_LB_loss = np.zeros_like(LB)
    updated_N_gain = np.zeros_like(N)
    updated_N_loss = np.zeros_like(N)



    def process_local(i):

        fpat = os.path.join(inpat, f'outputs/local_to_global_{i:06d}')
        with open(fpat, 'r') as f:
            lines = f.readlines()
        idx = 2
        n_faces = int(lines[idx].strip())
        idx += 1 + n_faces
        n_nodes = int(lines[idx].strip())
        idx += 1
        local_nodes = []
        global_nodes = []
        for j in range(n_nodes):
            local, global_ = map(int, lines[idx + j].split())
            local_nodes.append(local-1)
            global_nodes.append(global_-1)
        #print(f'Local to global mapping for local {i}: {len(local_nodes)} nodes')

        fpat = os.path.join(inpat, f'outputs/week_mean/mean_{i:06d}_{step}.h5')
        with h5py.File(fpat, 'r') as f:
            dt =  f['dt'][()]
            n_steps = f['n_steps'][()]
            mean_temp = f['mean_temp'][()]
            mean_tdry = f['mean_dry'][()]
            mean_wind = f['mean_wind'][()]
            mean_rad = f['mean_rad'][()]
        local_RB = RB[global_nodes]
        local_LB = LB[global_nodes]
        local_N = N[global_nodes]
        updated_RB_local = np.zeros_like(local_RB)
        updated_LB_local = np.zeros_like(local_LB)
        updated_N_local = np.zeros_like(local_N)
        if debug:
            debugthis = True
            updated_LB_gain_local = np.zeros_like(local_RB)
            updated_LB_loss_local = np.zeros_like(local_RB)
            updated_N_gain_local = np.zeros_like(local_RB)
            updated_N_loss_local = np.zeros_like(local_RB)
            updated_RB_gain_local = np.zeros_like(local_RB)
            updated_RB_loss_local = np.zeros_like(local_RB)
            
        else:
            debugthis = False

        #dt in days!
        for j in range(np.size(local_nodes)):
            if debugthis:
                updated_LB_local[j], updated_N_local[j], updated_RB_local[j], \
                updated_LB_gain_local[j], updated_LB_loss_local[j], updated_N_gain_local[j],\
                updated_N_loss_local[j], updated_RB_gain_local[j], updated_RB_loss_local[j] = sgint(
                    sgint_dt, local_LB[j], local_N[j], local_RB[j],
                    mean_temp[j], mean_tdry[j], mean_rad[j], mean_wind[j], 7*step,inpat,debug=debugthis
                )
            else:
            #dt in days!
                updated_LB_local[j], updated_N_local[j], updated_RB_local[j] = sgint(
                    sgint_dt, local_LB[j], local_N[j], local_RB[j], 
                    mean_temp[j], mean_tdry[j], mean_rad[j], mean_wind[j], 7*step,inpat,debug=debugthis
                )
        if debugthis:
            return global_nodes, updated_LB_local, updated_N_local, updated_RB_local, \
            updated_LB_gain_local, updated_LB_loss_local, updated_N_gain_local, updated_N_loss_local, \
            updated_RB_gain_local, updated_RB_loss_local
        else:
            return global_nodes, updated_LB_local, updated_N_local, updated_RB_local

    with ProcessPoolExecutor() as executor:
        results = list(executor.map(process_local, range(nfiles)))

    for global_nodes, updated_LB_local, updated_N_local, updated_RB_local, updated_LB_gain_local, updated_LB_loss_local, updated_N_gain_local, updated_N_loss_local, updated_RB_gain_local, updated_RB_loss_local in results:
        updated_LB[global_nodes] = updated_LB_local
        updated_RB[global_nodes] = updated_RB_local
        updated_N[global_nodes] = updated_N_local
        updated_LB_gain[global_nodes] = updated_LB_gain_local
        updated_LB_loss[global_nodes] = updated_LB_loss_local
        updated_N_gain[global_nodes] = updated_N_gain_local
        updated_N_loss[global_nodes] = updated_N_loss_local
        updated_RB_gain[global_nodes] = updated_RB_gain_local
        updated_RB_loss[global_nodes] = updated_RB_loss_local

    write_output_to_gr3(f'{outpat}/veg_LB_{step}.gr3', updated_LB,[])
    write_output_to_gr3(f'{outpat}/veg_RB_{step}.gr3', updated_RB,[])
    write_output_to_gr3(f'{outpat}/veg_N_{step}.gr3', updated_N,[])
    write_output_to_gr3(f'{outpat}/veg_h_{step}.gr3', updated_LB*a_ch,[])
    write_output_to_gr3(f'{outpat}/veg_LB_gain_{step}.gr3', updated_LB_gain,[])
    write_output_to_gr3(f'{outpat}/veg_LB_loss_{step}.gr3', updated_LB_loss,[])
    write_output_to_gr3(f'{outpat}/veg_N_gain_{step}.gr3', updated_N_gain,[])
    write_output_to_gr3(f'{outpat}/veg_N_loss_{step}.gr3', updated_N_loss,[])
    write_output_to_gr3(f'{outpat}/veg_RB_gain_{step}.gr3', updated_RB_gain,[])
    write_output_to_gr3(f'{outpat}/veg_RB_loss_{step}.gr3', updated_RB_loss,[])

print(f"Total runtime: {time.time() - start:.2f} seconds")
