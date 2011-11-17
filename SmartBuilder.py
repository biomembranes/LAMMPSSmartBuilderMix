import math
import random
import cfg
import os
import shutil
import re
import sys
from copy import copy
from copy import deepcopy

### Building of the forcefield from META file

def checkpath(filename):
	try:
		f = open(filename)
		return f
	except IOError:
    		print "File not found: error on %s " % filename
		sys.exit()
	except:
    		print "Unexpected error:", sys.exc_info()[0]
    		raise
	

def getatommasses():
	massdict = {}
	f = checkpath(cfg.inputforcefield)	
	foundmasses = reading = False
	for line in f:
		
		if ("mass" in line) and not foundmasses:
			foundmasses = True
		if foundmasses:
			if (re.search('\d+', line)) and not "#" in line.split()[0]: 	### End reading on blank line
				reading = True
				masstype = int(line.split()[1])
				massdict[masstype] = float(line.split()[2])
			if line.strip() == '':						### Check for blank line and end reading
				if reading:
					foundmasses = False
					reading = False	
	f.close()
	return massdict

def findoutermostmol(atomoffset, systemmoldictbyatom):
	## Get the outermost position + sigma (generous)
	maxz = sigma = 0.0		
	for atom in range(1, atomoffset+1):
		if systemmoldictbyatom[atom][1][2] > maxz:
			maxz = systemmoldictbyatom[atom][1][2]
			myatomtype = systemmoldictbyatom[atom][0]
			atomtypes = getatomtypes()
			LJsinglearray = readLJsingles(atomtypes)
			sigma = LJsinglearray[myatomtype][2]
	return (abs(maxz)+abs(sigma))
 		
	
def createwaters(noofwaters, systemdict, atomoffset, solarray, systemmoldictbyatom):
	allsoldict = {}
	area = float(cfg.xL)*float(cfg.yL)
	waterdim = math.pow(cfg.optimvolume, 1.0/3.0)
	squarewater = int (float(cfg.xL)/float(waterdim))
	xstart = -float(cfg.xL)/2.0; ystart = -float(cfg.yL)/2.0
	if atomoffset > 0:
		zstart = abs(findoutermostmol(atomoffset, systemmoldictbyatom))
	else:
		zstart = 0.0
	atom = atomoffset+1
	xint = yint = maxz_new = 0
	## Create atom positions first
	while atom < noofwaters+atomoffset+1:
		xint = (xint + 1) % int(squarewater) 
		yint = ((int((atom-atomoffset)/squarewater)) % int(squarewater))
		zint = int((atom-atomoffset)/squarewater**2)*(1-2*(atom % 2))
		xnew = xstart + waterdim*xint
		ynew = ystart + waterdim*yint
		znew = zstart*(1-2*(atom % 2)) + waterdim*zint	
		tempsol = deepcopy(solarray)
		tempsol[1][1] = [xnew, ynew, znew]
		tempsol[1][2] = int(cfg.molrep)
		allsoldict[atom-atomoffset] = tempsol[1]
		if znew > maxz_new:
			maxz_new = znew
		atom+=1
	if maxz_new > float(cfg.zL)/2.0:
		print "Please revise your zL, the water pushes zL to: %f\n" % (maxz_new*2)
	return allsoldict, atom-1


### Building of the forcefield from META file

def createnewFF(forcefieldLJ, forcefieldbonds):
	f = checkpath(cfg.inputforcefield)	
	w = open(cfg.outputforcefield, "w")
	fftext = ""
	foundpairLJ = False; foundpairbond = False
	for line in f:
		if ("pair_coeff" in line) and not foundpairLJ:
			foundpairLJ = True
			fftext += forcefieldLJ
		if ("bond_coeff" in line) and not foundpairbond:
			foundpairbond = True
			fftext += forcefieldbonds
		if not ("pair_coeff" in line) and not ("bond_coeff" in line):
			fftext += line

	if not foundpairLJ:
		### dump anyhow to make sure
		fftext += forcefieldLJ
	
	w.write(fftext)
	f.close()
	w.close()
	return 0

def getatomtypes():
	f = checkpath(cfg.inputmetafile)
	collect = reading = False
	type_no = 0
	for line in f:
		if "LJData" in line:
			collect = True
		if collect: 
			if (re.search('\d+', line)):
				reading = True
				atom_data = line.split()
			else:
				if reading:
					collect = False
	f.close()
	return int(atom_data[0])

def getmaxbondtypes(systemdict, species):
	f = checkpath(cfg.inputforcefield)	
	foundpairbond = reading = False
	maxbond = 0
	for line in f:
		if ("bond_coeff" in line) and not foundpairbond:
			foundpairbond = True
		if foundpairbond:
			if (re.search('\d+', line)) and not "#" in line.split()[0]: 	### End reading on blank line
				reading = True
				if int(line.split()[1]) > maxbond:
					maxbond = int(line.split()[1])
			if line.strip() == '':						### Check for blank line and end reading
				if reading:
					foundpairbond = False
					reading = False	
	f.close()
	return maxbond

def getmaxbendtypes(systemdict, species):
	f = checkpath(cfg.inputforcefield)	
	foundpairbend = reading = False
	maxbend = 0
	for line in f:
		if ("angle_coeff" in line) and not foundpairbend:
			foundpairbend = True
		if foundpairbend:
			if (re.search('\d+', line)) and not "#" in line.split()[0]: 	### End reading on blank line
				reading = True
				if int(line.split()[1]) > maxbend:
					maxbend = int(line.split()[1])
			if line.strip() == '':						### Check for blank line and end reading
				if reading:
					foundpairbend = False
					reading = False	
	f.close()
	return maxbend

		
def readLJsingles(atomtypes):
	LJdict = {}
	f = checkpath(cfg.inputmetafile)
	collect = False
	type_no = 0
	for line in f:
		if "LJData" in line:
			collect = True
		if collect and (re.search('\d+', line)) and (int(line.split()[0]) <= int(atomtypes)):
			atom_data = line.split()
			type_no = int(atom_data[0])
			LJdict[type_no] = [int(atom_data[0]), float(atom_data[1]), float(atom_data[2])]
			if int(line.split()[0]) == int(atomtypes):
				collect = False
	f.close()
	return LJdict

def readLJpairs(atomtypes):
	LJscaledict = {}
	f = checkpath(cfg.inputmetafile)
	pair_no = 0
	collect = False
	for line in f:
		if "Scaling" in line:
			collect = True
		if collect and (re.search('\d+', line)) and (int(line.split()[0]) <= int(atomtypes)):
			atom_data = line.split()
			LJscaledict[pair_no] = [int(atom_data[0]), int(atom_data[1]), float(atom_data[2])]                               
			pair_no += 1
			if line == "":
				collect = False
	f.close()
	#print LJscaledict
	#print "======================================="
	return LJscaledict, pair_no

def getscalefactor(i, j, LJpairs, no_pairs):
	scale = 1.0
	for pair in range(0, no_pairs):
		if int(LJpairs[pair][0]) == int(i) and int(LJpairs[pair][1]) == int(j):
			scale = LJpairs[pair][2]
	return scale	
	

def buildLJpairs(atomtypes, LJarray, LJpairs, no_pairs):
	pairdata = ""
	paircount = 0
	#print no_pairs
	for i in range(0, atomtypes+1):
		for j in range(i,atomtypes+1):
			scale = getscalefactor(i, j, LJpairs, no_pairs)
			if i in LJarray and j in LJarray:
				epsilon = scale * math.sqrt(float(LJarray[i][1])*float(LJarray[j][1]))
				sigma = (float(LJarray[i][2]) + float(LJarray[j][2]))/2.0
				pairdata += 'pair_coeff %7s %7s   %7f   %7f\n' % (i, j, epsilon, sigma)
			
	#print pairdata
	return pairdata	

	
def getatomno(atomspermol, molarray, atom_type):
	for atom_no in range(1, atomspermol+1):
		if (int(molarray[atom_no][0]) == int(atom_type)):
			return int(atom_no)

def geti_LJarray(atomtypes, LJmetaarray, atom_type):
	for pos in range(0, atomtypes):
		if (int(LJmetaarray[pos][0]) == int(atom_type)):
			return int(pos)
		
def getatomids_bondarray(bondspermol, bondarray, bondtype):
	i = j = None
	for bond in range(1, bondspermol+1):
		if int(bondarray[bond][0]) == int(bondtype):
			i = int(bondarray[bond][1])
			j = int(bondarray[bond][2])	
	return i, j

def getrealsigma_meta(atomtypes, atom_type):
	LJmetaarray = readLJsingles(atomtypes)
	posi = geti_LJarray(atomtypes, LJmetaarray, atom_type)
	sigma = float(LJmetaarray[posi][2])
	return sigma
	
def readcomputebondpairs(atomspermol, bondspermol, atomtypes, bondtypes, bondarray, molarray, scale):
	bondscaledsigmaarray = [None]*(int(bondtypes)+1)
	for bondtype in range(1,bondtypes+1):
		atomi, atomj = getatomids_bondarray(bondspermol, bondarray, bondtype)
		if atomi != None:
			LJmetaarray = readLJsingles(atomtypes)
			typei = molarray[atomi][0]; typej = molarray[atomj][0]
			#pos_i = geti_LJarray(atomtypes, LJmetaarray, typei)
			#pos_j = geti_LJarray(atomtypes, LJmetaarray, typej)
			sigmai = float(LJmetaarray[typei][2])
			sigmaj = float(LJmetaarray[typej][2])
 			scaledsigma = (sigmai+sigmaj)*(scale/2.0)
			bondscaledsigmaarray[bondtype] = scaledsigma
	return bondscaledsigmaarray

def createLAMMPSforcefieldscaledbond(bondscaledsigmaarray):
	f = checkpath(cfg.inputforcefield)
	bondtypedata = ""
	for line in f:
		if "bond_coeff" in line:
			pair_data = line.split()
			bond_lookup = int(pair_data[1])
			if bondscaledsigmaarray[bond_lookup] != None:
				bondtypedata += 'bond_coeff %7d   %7f   %7f\n' % (int(pair_data[1]), float(pair_data[2]), bondscaledsigmaarray[bond_lookup])
			else:
				bondtypedata += 'bond_coeff %7d   %7f   %7f\n' % (int(pair_data[1]), float(pair_data[2]), float(pair_data[3]))
	return bondtypedata	
	f.close()

### Create a data file based on N molecules



def simpleinvertmolecule(atomspermol, singlemolarray):
	newmolarray = copymol(atomspermol, singlemolarray)
	for atom in range(1, atomspermol+1):	
		newmolarray[atom][1][2] *= -1	# z positions
		newmolarray[atom][4][2] *= -1	# z positions (dipole)
	return newmolarray

def invertmolecule(atomspermol, singlemolarray, systemmoldict, intspecies):
	## To get interdigitation we swap the rotate about the xy plane
	## we also rotate lipid to allow for closer packing of asymmetric
	## tails.
	## We typically use this for Ceramide moeties since we have two heads.
	## See <simpleinvertmolecule> for simple version that flips and does not
	## worry about chain flipping etc.
	#print systemmoldict
	newmolarray = copymol(atomspermol, systemmoldict[intspecies])
	heads = findlipidheadbeads(atomspermol, newmolarray)
	#print heads
	heada = heads[0]; headb = heads[1]
	lena = headb - heada; lenb = atomspermol - lena
	chainaxyz = singlemolarray[heada][1]; chainbxyz = singlemolarray[headb][1]
	for atom in range(heada, headb):
		newmolarray[atom][1][0] = chainbxyz[0]
		newmolarray[atom][1][1] = chainbxyz[1]

	for atom in range(headb, atomspermol+1):
		newmolarray[atom][1][0] = chainaxyz[0]
		newmolarray[atom][1][1] = chainaxyz[1]
	
	for atom in range(1, atomspermol+1):	
		newmolarray[atom][1][2] *= -1	# z positions
		newmolarray[atom][4][2] *= -1	# z positions (dipole)
	return newmolarray
	
def findlipidheadbeads(atomspermol, singlemolarray):
	## This attempts to locate the 'lipid' two head beads
	## The head beads can then used to help swap the chains over 
	## and manipulate the lipids.
	# print atomspermol
	for atom in range(2, atomspermol):
		if ((singlemolarray[atom+1][1][2] - singlemolarray[atom][1][2]) > 0.0):
			currentatom = atom 
	return 1, currentatom + 1
	
			
def getsinglelipidprops():
	f = checkpath(cfg.inputsinglemoleculefile)
	for line in f:
		if "atoms" in line:
			atoms = re.search('\d+', line).group()
		if "bonds" in line:
			bonds = re.search('\d+', line).group()	
		if "angles" in line:
			angles = re.search('\d+', line).group()
	
	print "We have the following structure: " + " " + str(angles) + " " + str(bonds) + " " + str(angles) 
	f.close()
	return int(atoms), int(bonds), int(angles)

def readfullheader(filename):
	f = checkpath(filename)
	atoms = bonds = angles = dihedrals = impropers = atomtypes = bondtypes = angletypes = dihedraltypes = impropertypes = 0

	for line in f:
		if "atoms" in line:
			atoms = int(re.search('\d+', line).group())
		if "bonds" in line:
			bonds = int(re.search('\d+', line).group())	
		if "angles" in line:
			angles = int(re.search('\d+', line).group())
		if "dihedrals" in line:
			dihedrals = int(re.search('\d+', line).group())
		if "impropers" in line:
			impropers = int(re.search('\d+', line).group())
		if "atom types" in line:
			atomtypes = int(re.search('\d+', line).group())	
		if "bond types" in line:
			bondtypes = int(re.search('\d+', line).group())
		if "angle types" in line:
			angletypes = int(re.search('\d+', line).group())
		if "dihedral types" in line:
			dihedraltypes = int(re.search('\d+', line).group())
		if "improper types" in line:
			impropertypes = int(re.search('\d+', line).group())
	f.close()
	return atoms, bonds, angles, dihedrals, impropers, atomtypes, bondtypes, angletypes, dihedraltypes, impropertypes

def createfullheader(molecules, atoms, bonds, angles, dihedrals, impropers, atomtypes, bondtypes, angletypes, dihedraltypes, impropertypes):
	headdata = ""
	headdata += "%12s atoms\n" % atoms
	headdata += "%12s bonds\n" % bonds
	headdata += "%12s angles\n" % angles
	headdata += "%12s dihedrals\n" % dihedrals
	headdata += "%12s impropers\n\n" % impropers
	headdata += "%12s atom types\n" % atomtypes
	headdata += "%12s bond types\n" % bondtypes
	headdata += "%12s angle types\n" % angletypes
	headdata += "%12s dihedral types\n" % dihedraltypes
	headdata += "%12s improper types\n" % impropertypes
	return str(headdata)

def createboxdimensions(xL, yL, zL):
	xlo = -abs(xL/2.0); ylo = -abs(yL/2.0); zlo = -abs(zL/2.0)
	xhi = -xlo; yhi = -ylo; zhi = -zlo
	boxdata = ""
	boxdata += '%7.3f %7.3f  xlo xhi\n' % (xlo,xhi)
	boxdata += '%7.3f %7.3f  ylo yhi\n' % (ylo,yhi)
	boxdata += '%7.3f %7.3f  zlo zhi\n' % (zlo,zhi)
	return str(boxdata)

def readsinglemolecule(atomspermol, filename):
	moldict = {}
	f = checkpath(filename)
	collect = False
	for line in f:
		if "Atoms" in line:
			collect = True
		if collect and (re.search('\d+', line)) and (int(line.split()[0]) <= int(atomspermol)):
			atom_data = line.split()
			atom_no = int(atom_data[0])
			moldict[atom_no] = [int(atom_data[1]), [float(atom_data[2]), float(atom_data[3]), float(atom_data[4])], int(atom_data[5]), float(atom_data[6]), [float(atom_data[7]), float(atom_data[8]), float(atom_data[9])], float(atom_data[10]), float(atom_data[11])]
			if int(line.split()[0]) == int(atomspermol):
				collect = False
	f.close()
	return moldict

def readsinglemoleculebonds(bondspermol, filename):
	bondarray = [[0 for col in range(cfg.ROWBOND) ] for row in range(1, bondspermol+2)]
	f = checkpath(filename)
	collect = False
	for line in f:
		if "Bonds" in line:
			collect = True
		if collect and (re.search('\d+', line)) and (int(line.split()[0]) <= int(bondspermol)):
			bond_data = line.split()
			bond_no = int(bond_data[0])
			for row in range(0,cfg.ROWBOND):
				bondarray[bond_no][row] = bond_data[row]
			if int(line.split()[0]) == int(bondspermol):
				collect = False
	f.close()
	return bondarray

def readsinglemoleculeangles(anglespermol, filename):
	anglesarray = [[0 for col in range(cfg.ROWANGLE) ] for row in range(1, anglespermol+2)]
	f = checkpath(filename)
	collect = False
	for line in f:
		if "Angles" in line:
			collect = True
		if collect and (re.search('\d+', line)) and (int(line.split()[0]) <= int(anglespermol)):
			angles_data = line.split()
			angles_no = int(angles_data[0])
			for row in range(0,cfg.ROWANGLE):
				anglesarray[angles_no][row] = angles_data[row]
			if int(line.split()[0]) == int(anglespermol):
				collect = False
	f.close()
	return anglesarray

def readsinglemoleculedihedrals(dihedralspermol, filename):
	dihedarray = [[0 for col in range(cfg.ROWDIHEDRAL) ] for row in range(1, dihedralspermol+2)]
	f = checkpath(filename)
	collect = False
	for line in f:
		if "Dihedrals" in line:
			collect = True
		if collect and (re.search('\d+', line)) and (int(line.split()[0]) <= int(dihedralspermol)):
			dihed_data = line.split()
			dihed_no = int(dihed_data[0])
			for row in range(0,cfg.ROWDIHEDRAL):
				dihedarray[dihed_no][row] = dihed_data[row]
			if int(line.split()[0]) == int(dihedralspermol):
				collect = False
	f.close()
	return dihedarray

def readsinglemoleculeimpropers(improperspermol, filename):
	impropsarray = [[0 for col in range(cfg.ROWIMPROPERS) ] for row in range(1, improperspermol+2)]
	f = checkpath(filename)
	collect = False
	for line in f:
		if "Impropers" in line:
			collect = True
		if collect and (re.search('\d+', line)) and (int(line.split()[0]) <= int(improperspermol)):
			improp_data = line.split()
			improp_no = int(improp_data[0])
			for row in range(0,cfg.ROWIMPROPERS):
				improparray[dihed_no][row] = improp_data[row]
			if int(line.split()[0]) == int(improperspermol):
				collect = False
	f.close()
	return impropsarray

def copymol(atomspermol, molarray):
	newmolarray = deepcopy(molarray)
	return newmolarray

def createspeciesseq(molecules, nospecies):
	currentcount = deepcopy(molecules)
	monolayer = []; molcount = 1; speciesseq = []
	for item in currentcount:
		monolayer.append(item/2)
	while (molcount <= int(cfg.totalmol)/2):
		for spec in range(0, (nospecies-1)):
			if (monolayer[spec] > 0):
				molcount += 1
				monolayer[spec] -= 1
				speciesseq.append(spec)
	
	if (cfg.randomphase):
		random.shuffle(speciesseq)
	#print speciesseq	
	#print speciesseq+speciesseq
	return speciesseq+speciesseq
		
		
def createbilayerproperties(systemdict, systemmoldict, speciesseq, LJdictbytype, atommassesall):
	newsystemmoldictbyatom = {}
	xstart = -float(cfg.xL)/2.0; ystart = -float(cfg.yL)/2.0
	zoffset = -float(cfg.zoffset)
	molcurrent = 1; atomcurrent = 1
	## Create atom positions first
	for monolayer in range (2):
		for xlip in range(0, int(cfg.squarelipid)):
			for ylip in range(0, int(cfg.squarelipid)):
				if (molcurrent <= cfg.totalmol):
					speciesint = speciesseq[molcurrent-1]
					atomspermol = systemdict[speciesint][0]
					monolayeroffset = zoffset*(1-2*(monolayer % 2))
					cpymolarray = deepcopy(systemmoldict[speciesint])
					#print systemmoldict
					#print "=========="
					#print cpymolarray
					#print "=========="

					if (monolayer % 2 == 1):
						cpymolarray = invertmolecule(atomspermol, cpymolarray, systemmoldict, speciesint) 
						# Can replace with ch.inversion <invertmolecule>	
					for atom in range(1, atomspermol+1):
						#print atom, atomtype
						atomtype = int(cpymolarray[atom][0])
					
						newmass = atommassesall[atomtype]
						newsigma = LJdictbytype[atomtype][2]
						newdensity = newmass/(4.0/3.0 * math.pi * (newsigma/2.0)**3)
						xnew = float(cpymolarray[atom][1][0]) + float(xstart) + float(cfg.delx)*xlip
						ynew = float(cpymolarray[atom][1][1]) + float(ystart) + float(cfg.dely)*ylip
						znew = float(cpymolarray[atom][1][2]) + monolayeroffset
						cpymolarray[atom][1] = [xnew, ynew, znew]
						cpymolarray[atom][2] = molcurrent
						cpymolarray[atom][5] = newsigma
						cpymolarray[atom][6] = newdensity 
						newsystemmoldictbyatom[atomcurrent] = cpymolarray[atom]
						atomcurrent+=1
					molcurrent+=1
	return newsystemmoldictbyatom, atomcurrent-1		

def createbilayerpropertiesbond(systemdict, systembonddict, speciesseq, molecules):
	newsystembonddictbybond = {}
	molcurrent = 1; bondcurrent = 1
	atomoffset = 0
	for molecule in range(1, int(molecules)+1):
		speciesint = speciesseq[molcurrent-1]
		atomspermol = systemdict[speciesint][0]
		bondspermol = systemdict[speciesint][1]
		cpybondarray = deepcopy(systembonddict[speciesint])
		
		for bond in range(1, bondspermol+1):
			newsystembonddictbybond[bondcurrent] = [int(cpybondarray[bond][1]), int(cpybondarray[bond][2])+atomoffset, int(cpybondarray[bond][3])+atomoffset]
			bondcurrent+=1
		molcurrent+=1
		atomoffset+=atomspermol
	return newsystembonddictbybond, bondcurrent-1


def createbilayerpropertiesbend(systemdict, systembenddict, speciesseq, molecules):
	newsystembenddictbybend = {}
	molcurrent = 1; bendcurrent = 1
	atomoffset = 0
	## Create atom positions first
	for molecule in range(1, int(molecules)+1):
		speciesint = speciesseq[molcurrent-1]
		atomspermol = systemdict[speciesint][0]
		bendspermol = systemdict[speciesint][2]
		cpybendarray = deepcopy(systembenddict[speciesint])
		for bend in range(1, bendspermol+1):
			newsystembenddictbybend[bendcurrent] = [int(cpybendarray[bend][1]), int(cpybendarray[bend][2]) + atomoffset, int(cpybendarray[bend][3]) + atomoffset, int(cpybendarray[bend][4]) + atomoffset]
			bendcurrent+=1
		molcurrent+=1
		atomoffset+=atomspermol
	return newsystembenddictbybend, bendcurrent-1

def createbilayerpropertiesdihedrals(systemdict, systemdiheddict, speciesseq, molecules):
	newsystemdiheddictbydihed = {}
	molcurrent = 1; dihedcurrent = 1
	atomoffset = 0
	## Create atom positions first
	for molecule in range(1, int(molecules)+1):
		speciesint = speciesseq[molcurrent-1]
		atomspermol = systemdict[speciesint][0]
		dihedspermol = systemdict[speciesint][3]
		cpydihedarray = deepcopy(systemdiheddict[speciesint])
		for dihed in range(1, dihedspermol+1):
			newsystemdiheddictbydihed[dihedcurrent] = [int(cpydihedarray[dihed][1]), int(cpydihedarray[dihed][2]) + atomoffset, int(cpydihedarray[dihed][3]) + atomoffset, int(cpydihedarray[dihed][4]) + atomoffset, int(cpydihedarray[dihed][5]) + atomoffset]
			dihedcurrent+=1
		molcurrent+=1
		atomoffset+=atomspermol
	return newsystemdiheddictbydihed, dihedcurrent-1

def createbilayerpropertiesimpropers(systemdict, systemimpropdict, speciesseq, molecules):
	newsystemimpropdictbyimprop = {}
	molcurrent = 1; impropcurrent = 1
	atomoffset = 0
	## Create atom positions first
	for molecule in range(1, int(molecules)+1):
		speciesint = speciesseq[molcurrent-1]
		atomspermol = systemdict[speciesint][0]
		impropspermol = systemdict[speciesint][4]
		cpydihedarray = deepcopy(systemimpropdict[speciesint])
		for improp in range(1, impropspermol+1):
			newsystemimpropdictbyimprop[impropcurrent] = [int(cpydihedarray[improp][1]), int(cpydihedarray[improp][2]) + atomoffset, int(cpydihedarray[improp][3]) + atomoffset, int(cpydihedarray[improp][4]) + atomoffset, int(cpydihedarray[improp][5]) + atomoffset]
			impropcurrent+=1
		molcurrent+=1
		atomoffset+=atomspermol
	return newsystemimpropdictbyimprop, impropcurrent-1


def createatomdataLAMMPS(totalatoms, molarray, offset):
	atomdata = ""
	for atom in range(1, totalatoms+1):
		atomdata += '%7d %7d %10.3f %10.3f %10.3f %7d %7d %10.3f %10.3f %10.3f %7.3f %7.3f \n' % (atom+offset, molarray[atom][0], molarray[atom][1][0], molarray[atom][1][1], molarray[atom][1][2], molarray[atom][2], molarray[atom][3], molarray[atom][4][0], molarray[atom][4][1], molarray[atom][4][2], molarray[atom][5], molarray[atom][6])
	return atomdata


def createbonddataLAMMPS(totalbonds, molarray, offset):
	bonddata = "Bonds\n\n"
	for bond in range(1, totalbonds+1):
		bonddata += '%7d %7d %7d %7d\n' % (bond+offset, molarray[bond][0], molarray[bond][1], molarray[bond][2])
	return bonddata

def createbenddataLAMMPS(totalbends, molarray, offset):
	benddata = "Angles\n\n"
	for bend in range(1, totalbends+1):
		benddata += '%7d %7d %7d %7d %7d\n' % (int(bend)+int(offset), int(molarray[bend][0]), int(molarray[bend][1]), int(molarray[bend][2]), int(molarray[bend][3]))	
	return benddata

def createdihedralsdataLAMMPS(totalbends, molarray, offset):
	benddata = ""	
	if totalbends > 0:
		benddata = "Dihedrals\n\n"
	for bend in range(1, totalbends+1):
		benddata += '%7d %7d %7d %7d %7d %7d\n' % (int(bend)+int(offset), int(molarray[bend][0]), int(molarray[bend][1]), int(molarray[bend][2]), int(molarray[bend][3]), int(molarray[bend][4]))	
	return benddata

def createimpropersdataLAMMPS(totalbends, molarray, offset):
	benddata = ""
	if totalbends > 0:
		benddata = "Impropers\n\n"
	for bend in range(1, totalbends+1):
		benddata += '%7d %7d %7d %7d %7d %7d\n' % (int(bend)+int(offset), int(molarray[bend][0]), int(molarray[bend][1]), int(molarray[bend][2]), int(molarray[bend][3]), int(molarray[bend][4]))	
	return benddata


def createnewfileLAMMPS(buildtxt):
	f = open(cfg.outputsystemfile, "w")
	f.write(buildtxt)
	f.close()

def membranestats():
	areaperlipid = float(cfg.xL)*float(cfg.yL)/(float(cfg.squarelipid)**2)
	return areaperlipid

def main():
	intspecies = 0; systemdict = {}; systemmoldict = {}; systembonddict = {}; systembenddict = {}; systemdiheddict = {}; systemimprop = {};
	atomtypesall = getatomtypes()
	atommassesall = getatommasses()
	LJarray = readLJsingles(atomtypesall)
	for strspecies in (cfg.inputsinglemoleculefile):
		atoms, bonds, angles, dihedrals, impropers, atomtypes, bondtypes, angletypes, dihedraltypes, impropertypes = readfullheader(strspecies)
		systemdict[intspecies] = readfullheader(strspecies)
		systemmoldict[intspecies] = readsinglemolecule(atoms, strspecies)
		systembonddict[intspecies] = readsinglemoleculebonds(bonds, strspecies)
		systembenddict[intspecies] = readsinglemoleculeangles(angles, strspecies)
		systemdiheddict[intspecies] = readsinglemoleculedihedrals(dihedrals, strspecies)
		systemimprop[intspecies] = readsinglemoleculedihedrals(impropers, strspecies)
		intspecies+=1
	speciesseq = createspeciesseq(cfg.molecules, intspecies+1)
	newsystemmoldictbyatom, totalatomssolute = createbilayerproperties(systemdict, systemmoldict, speciesseq, LJarray, atommassesall)
	newsystembonddictbybond, totalbondssolute = createbilayerpropertiesbond(systemdict, systembonddict, speciesseq, cfg.totalmol)
        newsystembenddictbybend, totalbendssolute = createbilayerpropertiesbend(systemdict, systembenddict, speciesseq, cfg.totalmol)
        newsystemdiheddictbydihed, totaldihedssolute = createbilayerpropertiesdihedrals(systemdict, systemdiheddict, speciesseq, cfg.totalmol)
        newsystemimpropdictbyimprop, totalimpropssolute = createbilayerpropertiesimpropers(systemdict, systemimprop, speciesseq, cfg.totalmol)
	atomoffset = totalatomssolute
	formatatommol = createatomdataLAMMPS(totalatomssolute, newsystemmoldictbyatom, 0)
	formatbondmol = createbonddataLAMMPS(totalbondssolute, newsystembonddictbybond, 0)
	formatbendmol = createbenddataLAMMPS(totalbendssolute, newsystembenddictbybend, 0)
	formatdihedmol = createdihedralsdataLAMMPS(totaldihedssolute, newsystemdiheddictbydihed, 0)
	formatimpropmol = createimpropersdataLAMMPS(totalimpropssolute, newsystemimpropdictbyimprop, 0)
	if int(cfg.noofwaters) > 0:
		solarray = readsinglemolecule(1, cfg.inputsinglesolventfile)
		solfield, atomoffset = createwaters(int(cfg.noofwaters), systemdict, totalatomssolute, solarray, newsystemmoldictbyatom)
		formatsolmols = createatomdataLAMMPS(cfg.noofwaters, solfield, totalatomssolute)
	
	bondtypessolute = getmaxbondtypes(systemdict, intspecies)
	bendtypessolute = getmaxbendtypes(systemdict, intspecies)

	formatheader = createfullheader(cfg.totalmol, atomoffset, totalbondssolute, totalbendssolute, totaldihedssolute, totalimpropssolute, atomtypesall, bondtypessolute, bendtypessolute, dihedraltypes, impropertypes)
	formatboxhead = createboxdimensions(cfg.xL, cfg.yL, cfg.zL)
	systemfilecontent = "LAMMPS Bilayer comprising " + str(cfg.totalmol) + " created from a single molecule\n"
	systemfilecontent += formatheader + "\n"
	systemfilecontent += formatboxhead + "\n"
	systemfilecontent += "Atoms\n\n"
	systemfilecontent += formatatommol
	if int(cfg.noofwaters) > 0:
		systemfilecontent += formatsolmols 
	systemfilecontent += "\n\n\n\n"
	systemfilecontent += formatbondmol + "\n\n\n\n"
	systemfilecontent += formatbendmol + "\n\n\n\n"
	systemfilecontent += formatdihedmol + "\n\n\n\n"
	systemfilecontent += formatimpropmol + "\n\n\n\n"
	createnewfileLAMMPS(systemfilecontent)
	if int(cfg.totalmol) > 0:
		print 'Creating %d molecules' % (cfg.totalmol)
		print 'Initial area per lipid: %f Angstroms sq.' % membranestats()
		volumesys = float(cfg.xL * cfg.yL * cfg.zL)
		volumesol = float(cfg.empericalvolume*cfg.noofwaters)
		print 'Initial volume per lipid: %s Angstroms cbd.' % ((volumesys - volumesol)/(cfg.totalmol))

	### Handle the forcefield creation
	### Read in existing FF and ensure we account for LJ scaling and create i, j pairs
	### Additionaly we scale LJ pairs and scale bonded interactions with scaling factor 
	### pulled in from 'cfg' file.

	
	
	LJpairs, no_pairs = readLJpairs(atomtypesall)
	forcefieldLJ = buildLJpairs(atomtypesall, LJarray, LJpairs, no_pairs)
	bondscaledsigmaarray = readcomputebondpairs(totalatomssolute, totalbondssolute, atomtypesall, bondtypessolute, newsystembonddictbybond, newsystemmoldictbyatom, float(cfg.bondscaling))
	forcefieldbonds = createLAMMPSforcefieldscaledbond(bondscaledsigmaarray)
	createnewFF(forcefieldLJ, forcefieldbonds)

if __name__ == "__main__":
	main()

	
