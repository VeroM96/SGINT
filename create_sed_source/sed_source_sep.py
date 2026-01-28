import netCDF4 as nc
import numpy as np
import pandas as pd
import sys
import yaml

runpat = sys.argv[1]
ts = int(sys.argv[2]) 
ilonlat = 1  # 1 if lon/lat, 0 if x/y (in hgrid.gr3)

#### PARAMETERS ############################
with open(f'{runpat}/sgint.yaml') as info:
    params = yaml.full_load(info)
idsed = params['id_bgPOC']  # index for water column sediment fraction (AG biomass material)
f_exp = params['f_exp_bg']  # fraction of above ground biomass lost to POC

#### FUNCTIONS #############################
idsed = 6  # index for bottom sediment fraction
## DO NOT RUN THIS SCRIPT DIRECTLY, IT IS PART OF A WORKFLOW
#If you run this script directly, it will overwrite the sediment fractions in the hotstart file, whithout keeping the original values.
# It is used to create an updated sediment source file for the sediment model.
# _sep: this version is for seperate sediment classes for seagrass AG and BG biomass
# BG (root and rhizome material) is added to sediment class 7
# AG (leaf material) is added to sediment class 6
###1. Calculate bed mass from hotstart file
#get info from hotstart file
fpat = f'{runpat}/hotstart.nc'
with nc.Dataset(fpat, 'r') as f:
    bed = f.variables['SED3D_bed'][:]
    bedfrac = f.variables['SED3D_bedfrac'][:]
bed_thck = bed[:,:,0]
bed_poro = bed[:,:,2]
nlev = bed.shape[1]

# Parse SAND_SRHO from sediment.in file
with open(f'{runpat}/sediment.nml', 'r') as f:
    for line in f:
        if 'Srho =' in line:
            values = line.split('=')[1].split('!')[0].strip().split(',')
            srho = np.array([float(v.replace('d0', '').replace('D0', '')) for v in values])
            print(srho)
            break

bed_mass = np.zeros_like(bedfrac)
for ised in range(len(srho)):
    bed_mass[:,:,ised] = bed_thck[:,0:nlev] * srho[ised] * (1 - bed_poro[:,0:nlev]) * bedfrac[:, 0:nlev, ised]


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
    z = grid['depth'].values  # positive downwards
    
    if read_tri:
        tri = np.array(tri) - 1  # Convert to zero-based indexing
        return x, y, z, tri
    else:
        return x, y, z


### get grid file for tri & Area
xll,yll,_,tri = read_gr3_file(runpat + f'/hgrid.gr3')
def ll2xy(lat, lon, midy, midx):
    #convert x,y from lat, lon to meters. use mean(lon),mean(lat) as reference x -> lon, y = lat
    #approximation is enough: 1deg lat = 111.319 km; 1 deg lon = 111.319 km * cos(phi)
    r_earth = 40075000 #m 
    x = (lon-midx)*r_earth/360*np.cos(midy*np.pi/180)
    y = (lat-midy)*r_earth/360
    return x, y

#calculate area of elements:
if ilonlat == 1:
    midx = np.mean(xll)
    midy = np.mean(yll)
    # Note: You'll need to implement ll2xy function for lat/lon conversion
    x, y = ll2xy(yll, xll, midy, midx)

# Calculate triangle areas using the shoelace formula
triArea = np.abs((x[tri[:, 0]] * (y[tri[:, 1]] - y[tri[:, 2]]) +
                    x[tri[:, 1]] * (y[tri[:, 2]] - y[tri[:, 0]]) +
                    x[tri[:, 2]] * (y[tri[:, 0]] - y[tri[:, 1]])) / 2)

###2. calculate POC from roots to input into sediment source
if ts > 0:
    #read seagrass leaf data from file
    _,_,RB0 = read_gr3_file(runpat + f'/outputs/week_mean/veg_RB_{ts-1}.gr3',read_tri=False) #gC
    _,_,RB1 = read_gr3_file(runpat + f'/outputs/week_mean/veg_RB_{ts}.gr3',read_tri=False)

    diff = RB0 - RB1 # positive -> biomass lost
    
    diff_elem = np.mean(diff[tri], axis=1)  # gC

POC = np.where(diff_elem > 0, diff_elem / 1000 * f_exp, 0)  # convert to kgC/m², 3% of root NPP is stored in sediment as POC long term
print(f'total change in bed POC: {np.sum(POC* triArea)} kgC')

###3. add BG biomass to sed fraction
bed_mass[:,0,idsed] += POC  # add POC to the POC sediment fraction (index idsed) in top most layer (?) check index!

#change in bedthickness due to POC
bed_thck_new = np.zeros_like(bed_thck)
for lev in range(nlev):
    for nsed in range(len(srho)):
        bed_thck_new[:,lev] += bed_mass[:,lev,nsed]/ (srho[nsed] * (1 - bed_poro[:,lev]))  # kgC/m^2 / (kg/m^3) = m

#convert back to sediment fraction
bedfrac_new = np.zeros_like(bedfrac)
for ised in range(len(srho)):
    bedfrac_new[:, 0:nlev, ised] = bed_mass[:,:,ised] /np.sum(bed_mass, axis=2)

# Open the hotstart file in write mode to update the SED3D_bedfrac variable
with nc.Dataset(fpat, 'r+') as f:
    # Update the SED3D_bedfrac variable
    f.variables['SED3D_bedfrac'][:] = bedfrac_new
    f.variables['SED3D_bed'][:,:,0] = bed_thck_new
    # Ensure data is written to disk
    f.sync()

print(f"Successfully updated SED3D_bedfrac in {fpat}")
print(f"Added POC mass to sediment fraction {idsed} (POC class)")