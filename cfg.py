import math

### I/O Filenames

inputmetafile = "meta.sys"
inputsinglemoleculefile = ["data.dopcsingle"]
inputsinglesolventfile = "data.watersingle"
outputsystemfile = "data.sys"
inputforcefield = "forcefield.dopc" 
outputforcefield = "forcefield.sys"

### Box dimensions
### LAMMPS dimensions from -x to +x where x = xL/2

xL = yL = 70.0    # Set for a square system
zL = 90.0         #
zoffset = [-2.5]  # Controls interdigitation of lipids (+/- digitates/separates)

### Solute Molecules
### ================
### For best initial packing use numbers of the form: 2*N*N
### So 128 would be 2*8*8

molecules = [128]
alternatephase = False	# Set Multi-species systems to checkerboard until one exhausted.
randomphase = False     #  Randomize the species - giving a new config every time.
invertmoleculetype = ["s"]

### For non-bilayer geometries

MixPhaseRandom = False
approxvolmol = [2000] ## User specified input to help give some idea of the initial
		      ## system start-up time (limited to single component for reporting).

### We divide the whole volume up and use it all while randomizing the phases
### of the solute and water. There is no rotational element as yet - which can
### be added with a rotation matrix.

## End non-bilayer geometries

## "s" means simple invert where the chains are not swapped
## "c" means invert and swap the chains to cope with interdigitation

totalmol = reduce(lambda a,b: a+b, molecules)

### Solvent
### -------

noofwaters = 2048         # Set to zero for water free   
optimvolume = 30.0	  # Initial idealized water volume - this is what you want the system to try and
			  # use to pack the water molecules. Tune this to help with system building.
empericalvolume = 30.0    # Theoretical water volume (vol waters would like to live at)
molrep = 0                # What we name the molecules

### BondScaling - Bonded NN pair scaling
### Performs scaling on NN bonds according the the dimensions (sigma) 
### of those constituent atoms bondscaling*(sigma_i + sigma_j)*(0.5)

bondscaling = 0.90
excludescaling = [] ## keeps the bond length as in the input file


### ---------------------------------------------------------------------------------------------------------

### End of user definable constants ###


if totalmol > 0:
	squarelipid = math.ceil(math.sqrt(float(totalmol/2.0)))
	delx = dely = xL/float(squarelipid)
 
### Data Formats

ROWDATA  = 7 ## This is the number of lumps of data we take from LAMMPS data file
ROWBOND  = 4 ## Bond data from the LAMMPS data file
ROWANGLE = 5 ## Angle data from the LAMMPS data file
ROWDIHEDRAL = ROWIMPROPER = 6 ## Dihedral data from the LAMMPS data file



