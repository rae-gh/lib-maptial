#! /usr/bin/env python

import sys
import os
from pathlib import Path
import sys
DATADIR = str(Path(__file__).resolve().parent.parent )+ "/data/"
CODEDIR = str(Path(__file__).resolve().parent.parent )+ ""
sys.path.append(CODEDIR)
print(sys.path)

from prometry import pdbloader as pl
from prometry import pdbgeometry as pg
from maptial.map import mapsmanager as mman
from maptial.map import mapfunctions as mfun
from maptial.xyz import vectorthree as v3
from maptial.map import mapplothelp as mph
import numpy as np



pdb = "6eex"
print(pdb,DATADIR)

mman.MapsManager().set_dir(DATADIR)

print("Loading",pdb)
pla = pl.PdbLoader(pdb,DATADIR,cif=False)    
po = pla.load_pdb()









