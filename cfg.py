import math

### I/O Filenames

inputmetafile = "meta.sys"

inputsinglemoleculefile = ["data.dopcsingle", "data.dopcsingle_cpy", "data.dopcsingle_cpy", "data.dopcsingle_cpy", "data.dopcsingle_cpy"]
inputsinglesolventfile = "data.watersingle"
outputsystemfile = "data.sys"

inputforcefield = "forcefield.dopc" 
outputforcefield = "forcefield.dopc_scaled"

### Box dimensions
### LAMMPS goes from -x to +x where x = xL/2

xL = yL = 68.0    # Set for a square system
zL = 63.0         #
zoffset = -1.5    # Controls interdigitation of lipids (+/- digitates/separates)

### Molecules
### For best initial packing use numbers of the form: 2*N*N
### So 128 would be 2*8*8

molecules = [64, 66, 2, 2, 2]
totalmol = reduce(lambda a,b: a+b, molecules)
alternatephase = True	
randomphase = False


### Solvent

noofwaters = 1000   	  # Set to zero for water free   
optimvolume = 8.0	  # Initial idealized water volume
empericalvolume = 30.0    # Theoretical water volume
molrep = 0                # What we name the molecules

### BondScaling - Bonded NN pair scaling
### Performs scaling on NN bonds according the the dimensions (sigma) 
### of those constituent atoms bondscaling*(sigma_i + sigma_j)*(0.5)

bondscaling = 0.90


### End of user definable constants ###


### Lateral density of molecules

if totalmol > 0:
	
	squarelipid = math.ceil(math.sqrt(float(totalmol/2.0)))
	delx = dely = xL/float(squarelipid)
 
### Data Formats

ROWDATA  = 7 ## This is the number of lumps of data we take from LAMMPS data file
ROWBOND  = 4 ## Bond data from the LAMMPS data file
ROWANGLE = 5 ## Angle data from the LAMMPS data file
ROWDIHEDRAL = ROWIMPROPER = 6 ## Dihedral data from the LAMMPS data file

