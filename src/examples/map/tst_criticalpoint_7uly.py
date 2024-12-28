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

  # An example laplacian in peptide bonds in 7uly an electron crystallography structure
###### Edit the pdb and the atoms ######
pdb_code = "7uly"
central_atoms = ["A:5@C.A","A:6@C.A","A:7@C.A","A:8@C.A"]
linear_atoms = ["A:5@CA.A","A:6@CA.A","A:7@CA.A","A:8@CA.A"]
planar_atoms = ["A:5@O.A","A:6@O.A","A:7@O.A","A:8@O.A"]
interpolation = "bspline"
width = 4.5 #Angstrom
samples = 50
#### ----------- ###
mman.MapsManager().set_dir(DATADIR)
ml = mman.MapsManager().get_or_create(pdb_code,file=1,header=1,values=1)
mf = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj,interpolation)
slice_vectors = []
for i in range(len(central_atoms)):
  central_atom = central_atoms[i]
  linear_atom = linear_atoms[i]
  planar_atom = planar_atoms[i]  
  cc = v3.VectorThree().from_coords(ml.pobj.get_coords_key(central_atom))
  ll = v3.VectorThree().from_coords(ml.pobj.get_coords_key(linear_atom))
  pp = v3.VectorThree().from_coords(ml.pobj.get_coords_key(planar_atom))
  slice_vectors.append((cc,ll,pp))  

filename = "SHOW"

for cc,ll,pp in slice_vectors:
  vals2d1 = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=1,ret_type="2d")
  vals2d2 = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=2,ret_type="2d")
  vals2d3 = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=3,ret_type="2d")
  mplot = mph.MapPlotHelp(filename)
  mplot.make_plot_slice_2d(vals2d1,min_percent=1,max_percent=1,samples=samples,width=width,points=[cc,ll,pp],title=pdb_code, hue="BW",plot_type="heatmap")
  mplot.make_plot_slice_2d(vals2d2,min_percent=1,max_percent=1,samples=samples,width=width,points=[cc,ll,pp],title=pdb_code, hue="RB")
  mplot.make_plot_slice_2d(vals2d3,min_percent=1,max_percent=1,samples=samples,width=width,points=[cc,ll,pp],title=pdb_code, hue="CP",plot_type="contour")

  
