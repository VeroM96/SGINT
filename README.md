The intertidal seagrass model was developped for simulating the seasonality of the Zostera notlii meadows in the Wadden Sea (Mohr, 2026)\
It is written to be coupled with the 3D-hydrodynamic model SCHISM using the old output format. The workflow for the model is the following:

SCHISM -> calculate mean values for wind stress, water temperature, innundation time, radiation -> run SGINT -> hotstart SCHISM

if no feedback of the seagrass model on the hydrodynamics is needed, SGINT can be run offline after completing the SCHISM simulation. 

### Functions
**sgint.py**\
core SGINT function for single timestep & node

**run_sgint_parallel.py**\
function looping over SCHISM nodes

**make_week_mean_parallel.py**\
function creating weekly averaged values for wind stress, water temperature, innundation time and radiation from SCHISM output, 

### Folder structure in SCHISM run:
```
runpat
├─ sgint.yaml
├─ veg_*.gr3 (used for SCHISM input)
├─ restart_run.sh (for submitting the model workflow to the cluster)
├─ merge_hotstart.sh (batch script calling sgint functions during model iterations)
├─ ...
└─ outputs
  ├─ schout_**.nc
  ├─ local_to_global*
  ├─ hotstart*
  └─ week_mean
    ├─ mean_*.h5 (mean values for each timestep used as SGINT input)
    └─ veg_*_*.gr3 (output of SGINT for each timestep)
```

### additional code
**sed_source_sep.py:**\
converts below ground biomass lost between two SGINT timesteps to POC and adds it to the hotstart file as a sediment fraction in the bottom sediments

**create_source_nc_sep.py**\
converts above ground biomass lost between two SGINT timesteps to POC and adds it to the source file as a sediment fraction to be released during the next SCHISM run

**postprog:**
- sgint version that can be run when all mean values from SCHISM are already availiable. Includes timestepping and sgint version returns files with gain and loss during each timestep additionally to S, N, R

run_sgint_parallel_timestepping.py\
sgint_portprog.py
