import glob
import os, sys
from yt.mods import *

try:
    from mpi4py import MPI
    parallel = True
    comm = MPI.COMM_WORLD
    nprocs = comm.size
except:
    parallel = False
    nprocs = 1

if len(sys.argv) > 1:
    fname = [sys.argv[-1]]
    if not os.path.exists(fname[0]):
        fname = glob.glob("DD*/DD*.hierarchy")
        fname += glob.glob("RD*/RedshiftOutput*.hierarchy")
else:
    fname = glob.glob("DD*/DD*.hierarchy")
    fname += glob.glob("RD*/RedshiftOutput*.hierarchy")
    
if not os.path.exists("frames") and ytcfg.getint("yt","__global_parallel_rank") == 0:
    os.mkdir("frames")

fields = ['Density', 'Temperature', 'HII_Fraction', 'HI_Fraction', 'HI_kph']

zlim = {}
#zlim['Density'] = (1e-20, 5e-18)
#zlim['Temperature'] = (1e2, 3e4)
#zlim['HI_Fraction'] = (1e-5, 1)
#zlim['HII_Fraction'] = (1e-5, 1)
#zlim['HI_kph'] = (1e-6, 30)
#zlim['TotalMetallicity'] = (1e-6, 1)
#zlim[dm] = (5e-3, 5.0)

cmap = {}
cmap['Density'] = 'algae'
cmap['Temperature'] = 'hot'
cmap['HI_kph'] = 'idl05'
#cmap['TotalMetallicity'] = 'gist_stern'
#cmap[dm] = 'bone'

ts = TimeSeriesData.from_filenames(fname)
for pf in ts.piter():
    test_pic_name = "frames/%s_Slice_x_%s.png" % (pf, fields[0])
    if os.path.exists(test_pic_name): continue
    if 'HI_kph' not in pf.h.field_list:
        ff = fields[:]
        ff.remove('HI_kph')
    else:
        ff = fields
    p = SlicePlot(pf, 'x', ff, center=[0.49]*3)
    for f in fields:
        if f in zlim.keys():
            if f in fields:
                p.set_zlim(f, zlim[f][0], zlim[f][1])
        if f in cmap.keys():
            if f in fields:
                p.set_cmap(f, cmap[f])
    #p.annotate_velocity()
    p.annotate_grids()
    p.set_width((0.1, "pc"))
    p.save("frames/%s" % pf)

