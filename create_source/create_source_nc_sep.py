"""
This script calculates the seagrass leaf biomass that is input into the water column as a source term.
Better version: depth is calculated from mean values, source nc only has two time steps: one at start and one at end.
_sep: This version is for separate sediment classes for seagrass AG and BG biomass.
BG (root and rhizome material) is added to sediment class 6.
AG (leaf material) is added to sediment class 5.

Author: Veronika Mohr
Date: August 2025
"""

import netCDF4 as nc
import numpy as np
import pandas as pd
import sys
import yaml

runpat = sys.argv[1]
ts = int(sys.argv[2]) 

#### PARAMETERS ############################
with open(f'{runpat}/sgint.yaml') as info:
    params = yaml.full_load(info)
idsed = params['id_agPOC']  # index for water column sediment fraction (AG biomass material)
f_exp = params['f_exp_ag']  # fraction of above ground biomass lost to POC

#### FUNCTIONS #############################
def create_source_nc(filename='source.nc', 
                     nsources=1, 
                     ntracers=2, 
                     time_msource=None, 
                     time_vsource=None, 
                     time_step_seconds=120,
                     source_elements=None, 
                     vsource_data=None, 
                     msource_data=None, 
                     time_msource_data=None, 
                     time_vsource_data=None):
    """
    Create a NetCDF file with SCHISM source structure.
    
    Parameters:
    -----------
    filename : str
        Output filename for the NetCDF file
    nsources : int
        Number of sources
    ntracers : int
        Number of tracers
    time_msource : int
        Number of time steps for mass source
    time_vsource : int
        Number of time steps for volume source
    time_step_seconds : float
        Time step in seconds
    source_elements : array-like, optional
        Source element indices (1-based)
    vsource_data : array-like, optional
        Volume source data (time_vsource, nsources)
    msource_data : array-like, optional
        Mass source data (time_msource, ntracers, nsources)
    time_msource_data : array-like, optional
        Time array for mass source
    time_vsource_data : array-like, optional
        Time array for volume source
    
    Returns:
    --------
    str : Success message with file details
    """
    
    # Convert to integers if needed
    time_msource = int(time_msource) if time_msource is not None else None
    time_vsource = int(time_vsource) if time_vsource is not None else None
    
    # Create the NetCDF file
    dataset = nc.Dataset(filename, 'w', format='NETCDF4')
    
    try:
        # Create dimensions
        nsinks = 0  # UNLIMITED, currently 0
        time_vsink = 0  # UNLIMITED, currently 0
        one = 1

        # Define dimensions - make time dimensions unlimited for concatenation
        dataset.createDimension('nsources', nsources)
        dataset.createDimension('nsinks', None)  # UNLIMITED
        dataset.createDimension('ntracers', ntracers)
        dataset.createDimension('time_msource', None)  # UNLIMITED for concatenation
        dataset.createDimension('time_vsource', None)  # UNLIMITED for concatenation
        dataset.createDimension('time_vsink', None)  # UNLIMITED
        dataset.createDimension('one', one)

        # Create variables
        source_elem = dataset.createVariable('source_elem', 'i4', ('nsources',))
        vsource = dataset.createVariable('vsource', 'f4', ('time_vsource', 'nsources'))
        msource = dataset.createVariable('msource', 'f4', ('time_msource', 'ntracers', 'nsources'))
        time_msource_var = dataset.createVariable('time_msource', 'f8', ('time_msource',))
        time_vsource_var = dataset.createVariable('time_vsource', 'f8', ('time_vsource',))
        time_step_vsource = dataset.createVariable('time_step_vsource', 'f4', ('one',))
        time_step_msource = dataset.createVariable('time_step_msource', 'f4', ('one',))
        time_step_vsink = dataset.createVariable('time_step_vsink', 'f4', ('one',))


        # Initialize variables with provided data or defaults
        # Source elements (1-based indexing for SCHISM)
        if source_elements is not None:
            source_elem[:] = source_elements
        else:
            source_elem[:] = np.arange(1, nsources + 1)

        # Time arrays
        if time_msource_data is not None:
            time_msource_var[:time_msource] = time_msource_data
        else:
            # Default time array starting from 0
            time_msource_var[:time_msource] = np.arange(time_msource) * time_step_seconds
            
        if time_vsource_data is not None:
            time_vsource_var[:time_vsource] = time_vsource_data
        else:
            # Default time array starting from 0
            time_vsource_var[:time_vsource] = np.arange(time_vsource) * time_step_seconds

        # Time steps
        time_step_vsource[:] = time_step_seconds
        time_step_msource[:] = time_step_seconds
        time_step_vsink[:] = time_step_seconds

        # Initialize vsource and msource with provided data or zeros
        if vsource_data is not None:
            vsource[:time_vsource, :] = vsource_data
        else:
            vsource[:time_vsource, :] = 0.0
            
        if msource_data is not None:
            msource[:time_msource, :, :] = msource_data
        else:
            msource[:time_msource, :, :] = 0.0
        
        success_msg = f"NetCDF file '{filename}' created successfully!\n"
        success_msg += f"Dimensions:\n"
        success_msg += f"  nsources: {nsources}\n"
        success_msg += f"  nsinks: UNLIMITED (currently {nsinks})\n"
        success_msg += f"  ntracers: {ntracers}\n"
        success_msg += f"  time_msource: UNLIMITED (currently {time_msource})\n"
        success_msg += f"  time_vsource: UNLIMITED (currently {time_vsource})\n"
        success_msg += f"  time_vsink: UNLIMITED (currently {time_vsink})\n"
        success_msg += f"  one: {one}"
        
        return success_msg
        
    finally:
        # Close the dataset
        dataset.close()

def read_gr3_file(filepath, read_tri=True):
    """
    Read x, y, depth from gr3 file.
    
    Parameters:
    -----------
    filepath : str
        Path to the .gr3 file
    read_tri : bool, optional
        Whether to read triangle connectivity data (default: True)
    
    Returns:
    --------
    tuple : (x, y, z) if read_tri=False, (x, y, z, tri) if read_tri=True
        x, y : array-like
            Node coordinates
        z : array-like
            Depth values (positive upwards)
        tri : array-like
            Triangle connectivity (zero-based indexing)
    """
    with open(filepath) as f:
        f.readline()  # skip title
        nElems, nNodes = map(int, f.readline().split())
        coords = []
        tri = []
        for _ in range(nNodes):
            parts = f.readline().split()
            coords.append([int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])])
        if read_tri:
            for _ in range(nElems):
                parts = f.readline().split()
                tri.append([int(parts[2]), int(parts[3]), int(parts[4])])
    
    grid = pd.DataFrame(coords, columns=['node', 'x', 'y', 'depth'])
    x = grid['x'].values
    y = grid['y'].values
    z = grid['depth'].values  # positive upwards
    
    if read_tri:
        tri = np.array(tri) - 1  # Convert to zero-based indexing
        return x, y, z, tri
    else:
        return x, y, z
   
def ll2xy(lat, lon, midy, midx):
    '''
    convert x,y from lat, lon to meters. use mean(lon),mean(lat) as reference x -> lon, y = lat
    approximation is enough: 1deg lat = 111.319 km; 1 deg lon = 111.319 km * cos(phi)
    '''
    r_earth = 40075000 #m 
    x = (lon-midx)*r_earth/360*np.cos(midy*np.pi/180)
    y = (lat-midy)*r_earth/360
    return x, y

# Parse SAND_SRHO from sediment.nml file
with open(f'{runpat}/sediment.nml', 'r') as f:
    for line in f:
        if 'Srho =' in line:
            values = line.split('=')[1].split('!')[0].strip().split(',')
            srho = np.array([float(v.replace('d0', '').replace('D0', '')) for v in values])
            break

outpat = f'{runpat}/outputs/week_mean'
xll,yll,_,tri=  read_gr3_file(runpat + '/hgrid.gr3')
ilonlat = 1  # 1 if lon/lat, 0 if x/y
nElems = len(tri)  # number of elements, same as number of sources
ntracers = 2+len(srho) #T,S, n*SED
time_step_seconds = 31536000 #seconds in a year
time_msource = 2  # 1 year in terms of time steps
time_hotstart = 604800   # length of hotstart in seconds (7 days)


#create sample msource data
msource_data = np.zeros([int(time_msource), ntracers, nElems])
vsource_data = np.zeros([int(time_msource), nElems])

#calculate area of elements
if ilonlat == 1:
    midx = np.mean(xll)
    midy = np.mean(yll)
    # Note: You'll need to implement ll2xy function for lat/lon conversion
    x, y = ll2xy(yll, xll, midy, midx)

# Calculate triangle areas using the shoelace formula
triArea = np.abs((x[tri[:, 0]] * (y[tri[:, 1]] - y[tri[:, 2]]) +
                    x[tri[:, 1]] * (y[tri[:, 2]] - y[tri[:, 0]]) +
                    x[tri[:, 2]] * (y[tri[:, 0]] - y[tri[:, 1]])) / 2)


if ts > 0:

    #read seagrass leaf data from file
    _,_,LB0 = read_gr3_file(runpat + f'/outputs/week_mean/veg_LB_{ts-1}.gr3',read_tri=False) #gC
    _,_,N0 = read_gr3_file(runpat + f'/outputs/week_mean/veg_N_{ts-1}.gr3',read_tri=False) #shoots/m^2
    _,_,LB1 = read_gr3_file(runpat + f'/outputs/week_mean/veg_LB_{ts}.gr3',read_tri=False)
    _,_,N1 = read_gr3_file(runpat + f'/outputs/week_mean/veg_N_{ts}.gr3',read_tri=False)
    diff = LB0*N0 - LB1*N1 # positive -> biomass lost
    #reshape to element shape
    elem_diff = np.mean(diff[tri], axis=1)*f_exp #gC/m^2, f_exp of leaf NPP is exported from meadow as DOC & POC
    elem_diff = elem_diff / 1000 * triArea  # convert to kgC
    velem = elem_diff / srho[idsed]  # get volume in m^3
    #if values are positive, calculate the loss per time step
    velem = np.where(velem > 0, velem, 0)  # m^3
    tmp = np.where(elem_diff > 0, elem_diff, 0)  # SED4 #kgC/m^3
    print(f"Total biomass loss in kgC/m^3: {np.sum(tmp)}")

    msource_data[:,idsed+2,:] = srho[idsed]  # SED4
    vsource_data[:, :] = np.tile(velem/time_hotstart, (int(time_msource), 1)) # SED4, volume source data in m^3/s


# create sample vsource data
msg = create_source_nc(filename=outpat+f'/source_{ts}.nc', 
                     nsources=nElems,
                     ntracers=ntracers, 
                     time_msource=time_msource,
                     time_vsource=time_msource,
                     time_step_seconds=time_step_seconds,
                     msource_data=msource_data,
                     vsource_data=vsource_data)
print(msg)
