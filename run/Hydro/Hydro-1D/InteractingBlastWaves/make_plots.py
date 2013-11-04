# See Figure 9 of Stone et al. 2008 ApJS 178 137

from yt.mods import *
import pylab

### define problem name
problem_name = 'InteractingBlastWaves'
#problem_name = 'SodShockTubeAMR'


### define simulation output directory and filename base
output_dir_base = 'DD'
datafile_base = 'data'


### define output to be plotted
dumpid = '0039'


### construct input filename
filename = './' + output_dir_base + dumpid + '/' + datafile_base + dumpid
print "Plotting output file %s\n" % filename


png_filename = './' + problem_name + '.png'

### load data
pf = load(filename)


### define InternalEnergy field
def _InternalEnergy(field, data):
    return data['TotalEnergy'] - 0.5*data['x-velocity']*data['x-velocity']

add_field('InternalEnergy', function=_InternalEnergy, units=r'\rm{erg}/\rm{g}')

### define Pressure field
def _Pressure(field, data):
    return (data.pf['Gamma'] - 1.0) * data['Density'] * data['InternalEnergy']

add_field('Pressure', function=_Pressure, units=r'\rm{dyne}/\rm{cm}^{2}')

### extract an ortho_ray (1D solution vector)
ray = pf.h.ortho_ray(0, [0.5, 0.5])

### define fields vector
fields = ('Density', 'x-velocity', 'InternalEnergy', 'Pressure' )

### make plot

pylab.figure()

# Density Plot
pylab.plot(ray['x'],ray['Density'], 'ro', ms=4)
pylab.xlim((0.4,1.0))
pylab.ylim((0.0,7.0))
pylab.xlabel('Position')
pylab.ylabel('Density')

### Save plot
pylab.savefig(png_filename)
