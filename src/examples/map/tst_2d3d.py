#! /usr/bin/env python

import sys
import os
from pathlib import Path
import sys
DATADIR = str(Path(__file__).resolve().parent.parent.parent )+ "/data/"
CODEDIR = str(Path(__file__).resolve().parent.parent.parent )+ ""
sys.path.append(CODEDIR)

from maptial.geo import pdbloader as pl
from maptial.map import mapsmanager as mman

pdb = "6eex"
print(pdb,DATADIR)

mman.MapsManager().set_dir(DATADIR)

print("Loading",pdb)
pla = pl.PdbLoader(pdb,DATADIR,cif=False)    
po = pla.load_pdb()









