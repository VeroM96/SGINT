import os
import numpy as np
from sgint import sgint
import h5py
from concurrent.futures import ProcessPoolExecutor
import sys
import glob
import yaml

#Author: Veronika Mohr, 18.12.2025
#run sgint using mean files from SCHISM output. Parallelizes the process over the number of local files created during SCHISM run.

#check if arguments are provided, else throw an error and exit
if len(sys.argv) < 3:
    sys.exit(f"Error run_sgint_parallel.py: Missing required arguments. Usage: python make_week_mean_parallel.py <runpat> <step>")

#read run pattern from command line argument
runpat = sys.argv[1]
step = int(sys.argv[2]) 

#check number of files from parallel jobs
nfiles = len(glob.glob(f'{runpat}/outputs/local_to_global*'))

#### FUNCTIONS #############################
#write output to veg_ .gr3 files
def write_output_to_gr3(filename, updated_veg):
    with open(f'{inpat}/sav_N.gr3', 'r') as fin, open(filename, 'w') as fout:
        # Copy header
        fin.readline()  # skip title
        print(f'Writing output to {filename.split("/")[-1]}')   
        fout.write(f'{filename.split("/")[-1]} at ts {step}\n')
        nElems, nNodes = map(int, fin.readline().split())
        fout.write(f"{nElems} {nNodes}\n")
        for idx in range(nNodes):
            line = fin.readline()
            parts = line.split()
            # Replace the last value with updated_N
            fout.write(f"{parts[0]} {parts[1]} {parts[2]} {updated_veg[idx]}\n")

#read input veg_ .gr3 files
def read_input_veg(filename):
    with open(filename) as f:
        f.readline()  # skip title
        _ , nNodes = map(int, f.readline().split())
        N = []
        for _ in range(nNodes):
            parts = f.readline().split()
            h = float(parts[3])
            N.append(h) 
    N = np.array(N)  # Convert to numpy array for easier manipulation
    return N
#### PARAMETERS ############################
with open(f'{runpat}/sgint.yaml') as info:
    params = yaml.full_load(info)
h0 = params['h0']               #minimum height of seagrass canopy [m]
a_ch = params['a_ch']           #conversion factor from gC to height [m gC^-1]
sgint_dt = params['sgint_dt']   #time step for sgint model [d]

#### OUTPUTS ############################
outpat = f'{runpat}/outputs/week_mean'

#### PATHS ############################
datapat = f'{runpat}/outputs/week_mean'
if not os.path.exists(datapat):
    sys.exit(f"Error run_sgint_parallel.py: The specified path '{datapat}' does not exist.")
inpat =  runpat

#copy input veg_N.gr3 to output folder
if step == 1:
    os.system(f'cp {inpat}/sav_N.gr3 {outpat}/sav_N_0.gr3')
    os.system(f'cp {inpat}/sav_RB.gr3 {outpat}/sav_RB_0.gr3')
    os.system(f'cp {inpat}/sav_LB.gr3 {outpat}/sav_LB_0.gr3')
    os.system(f'cp {inpat}/sav_h.gr3 {outpat}/sav_h_0.gr3')

#### INITIALIZATION ############################  

LB = read_input_veg(f'{outpat}/sav_LB_{step-1}.gr3')
RB = read_input_veg(f'{outpat}/sav_RB_{step-1}.gr3')
N = read_input_veg(f'{outpat}/sav_N_{step-1}.gr3')

updated_LB = np.zeros_like(LB)
updated_RB = np.zeros_like(RB)
updated_N = np.zeros_like(N)

#parallelizing
def process_local(i):
    fpat = f'{runpat}/outputs/local_to_global_{i:06d}'
    if not os.path.exists(fpat):
        sys.exit(f"Error run_sgint_parallel.py: The specified path '{fpat}' does not exist.")
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

    fpat = f'{datapat}/mean_{i:06d}_{step}.h5'
    if not os.path.exists(fpat):
        sys.exit(f"Error run_sgint_parallel.py: The specified path '{fpat}' does not exist.")
    with h5py.File(fpat, 'r') as f:
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

    for j in range(np.size(local_nodes)):
        if (local_LB[j] + local_RB[j] +local_N[j]) <=0:
            updated_LB_local[j] = 0
            updated_RB_local[j] = 0
            updated_N_local[j] = 0
            continue
        # Call the sgint function for each local node !!!dt in days!!!!
        updated_LB_local[j], updated_N_local[j], updated_RB_local[j] = sgint(
            sgint_dt, local_LB[j], local_N[j], local_RB[j],
            mean_temp[j], mean_tdry[j], mean_rad[j], mean_wind[j],7*step, runpat
        )
    return global_nodes, updated_LB_local, updated_N_local, updated_RB_local

with ProcessPoolExecutor() as executor:
    results = list(executor.map(process_local, range(nfiles)))

for global_nodes, updated_LB_local, updated_N_local, updated_RB_local in results:
    updated_LB[global_nodes] = updated_LB_local
    updated_RB[global_nodes] = updated_RB_local
    updated_N[global_nodes] = updated_N_local

write_output_to_gr3(f'{outpat}/sav_LB_{step}.gr3', updated_LB)
write_output_to_gr3(f'{outpat}/sav_RB_{step}.gr3', updated_RB)
write_output_to_gr3(f'{outpat}/sav_N_{step}.gr3', updated_N)
write_output_to_gr3(f'{outpat}/sav_h_{step}.gr3', updated_LB*a_ch)

os.system(f'cp {outpat}/sav_N_{step}.gr3 {inpat}/sav_N.gr3')
os.system(f'cp {outpat}/sav_RB_{step}.gr3 {inpat}/sav_RB.gr3')
os.system(f'cp {outpat}/sav_LB_{step}.gr3 {inpat}/sav_LB.gr3')
os.system(f'cp {outpat}/sav_h_{step}.gr3 {inpat}/sav_h.gr3')