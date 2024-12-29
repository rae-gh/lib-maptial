import sys
from pathlib import Path
import sys
DATADIR = str(Path(__file__).resolve().parent.parent.parent )+ "/data/"
print(DATADIR)
CODEDIR = str(Path(__file__).resolve().parent.parent.parent )+ ""
sys.path.append(CODEDIR)

###### Enter the pdb #######################
pdb_code = "6q53"
plane = "xy"
##########################################
import maptial.map.mapsmanager as mman
import maptial.map.mapfunctions as mfun
import matplotlib.pyplot as plt
mman.MapsManager().set_dir(DATADIR)
ml = mman.MapsManager().get_or_create(pdb_code,file=1,header=1,values=1,)
mf = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj,"nearest")
valsxy = mf.get_map_projection("xy")
xs,ys,zs,vs,xx,yy,zz = mf.get_atoms_projection("linear",log_level=1)
mf.get_map_cross_section("xy",0)
mf.get_map_projection("xy",0,0,0,0)
plt.scatter(ys,xs, c=vs,cmap="rainbow")