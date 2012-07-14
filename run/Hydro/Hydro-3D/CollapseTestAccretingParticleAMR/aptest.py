from yt.mods import *
import h5py
import asciitable

fns = sorted(glob.glob('./DD[0-9][0-9][0-9][0-9]/DD[0-9][0-9][0-9][0-9]'))

ParticleData = []
TotalMass = []
ParticleMass = []
Time = []
MinTemperature = []
MaxDensity = []

msun = 1.98892012e33

for i,fn in enumerate(fns):
    print '---------------------------------------------------------'
    print '---------------------------------------------------------'
    print fn
    pf = load(fn)
    dd = pf.h.all_data()
    TotalMass.append(dd.quantities["TotalQuantity"](["CellMassCode"])[0])
    MinTemperature.append(dd.quantities["Extrema"](["Temperature"])[0][0])
    MaxDensity.append(dd.quantities["Extrema"](["Density"])[0][1])
    cpufns = glob.glob(fn+'.cpu*')
    ParticleData.append([])
    Time.append(pf.current_time)
    dx = pf.domain_right_edge[0].astype('float64')/pf['TopGridDimensions'][0].astype('float64')/ \
        2.**pf['MaximumRefinementLevel']
    ParticleMass.append(0)
    for cpufn in cpufns:
        handle = h5py.File(cpufn,'r')
        keys = sorted(handle.keys())
        for key in keys:
            if key[0:4] != u'Grid':
                continue
            Grid = handle[key]
            gkeys = Grid.keys()
            if u'ActiveParticles' in gkeys:
                print ''
                AP = Grid['ActiveParticles']['AccretingParticle']
                np = AP[AP.keys()[0]].shape[0]
                print 'Found %(n)3d Active Particles in %(key)s' % {'n': np, 'key': key} 
                for j in na.arange(np):
                    ParticleDatum = {}
                    for apkey in AP.keys():
                        print apkey,AP[apkey][:][j]
                        ParticleDatum.setdefault(apkey,AP[apkey][:][j])
                    ParticleData[i].append(ParticleDatum)
                    ParticleMass[i]+=AP['mass'][:][j]*dx**3
        handle.close()
    c = [0.49,0.49,0.49]
    pc = PlotCollection(pf,c)
    slc = SlicePlot(pf, 0, 'Density', center=c, width=0.3)
    slc.set_zlim('Density',1e-20,3.165577e-17)
    slc.save('./frames/'+fn[9:])

TotalMass = na.array(TotalMass)
ParticleMass = na.array(ParticleMass)

print TotalMass
print ParticleMass
print TotalMass+ParticleMass
print zip(fns,ParticleMass)
print MinTemperature
print MaxDensity

asciitable.write({"TotalMass":TotalMass,
                  "ParticleMass":ParticleMass,
                  "Total+Particle":TotalMass+ParticleMass,
                  "MinTemperature":MinTemperature,
                  "Time":Time},
                 "masses.dat")
