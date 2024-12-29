#! /usr/bin/env python

import sys
from pathlib import Path
import sys
DATADIR = str(Path(__file__).resolve().parent.parent.parent )+ "/data/"
CODEDIR = str(Path(__file__).resolve().parent.parent.parent )+ ""
sys.path.append(CODEDIR)
from maptial.map import mapsmanager as mman
from maptial.xyz import vectorthree as v3
import maptial.map.mapsmanager as mman
import maptial.map.mapfunctions as mfun
import maptial.map.mapplothelp as mph

# And example of rings in disulfide bonds
###### Edit the pdb and the atoms ######
pdb_code = "1ejg"
central_atoms = ["A:3@SG.A","A:4@SG.A","A:16@SG.A"]
linear_atoms = ["A:40@SG.A","A:32@SG.A","A:26@SG.A"]
planar_atoms = ["A:3@O.A","A:4@O.A","A:16@O.A"]
interpolation = "bspline"
width = 10
samples = 200

mman.MapsManager().set_dir(DATADIR)
ml = mman.MapsManager().get_or_create(pdb_code,file=1,header=1,values=1,cif=False)
first_3 = ml.pobj.get_first_three()
print(first_3)



  