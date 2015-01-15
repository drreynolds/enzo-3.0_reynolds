from yt.mods import *
from matplotlib import use; use('Agg')

#Test that the correct number of Dark Matter particles are produced and 
#that their masses are as expected. 

CHILDRENPERPARENT = 12
File1 = "RD0000/RD0000"
File2 = "RD0001/RD0001"

centre = [0.5, 0.5, 0.5]
leftedge = [0.25, 0.25, 0.25]
rightedge = [0.75, 0.75, 0.75]
domainedgeleft = [0.0,0.0,0.0]
domainedgeright = [1.0,1.0,1.0]
pf = load(File1)
reg = pf.h.region(center=centre, left_edge=leftedge, right_edge=rightedge)

NumberofParents = len(reg["particle_mass"])
print "Number of Parents = ", NumberofParents
ParentMass = reg["ParticleMassMsun"][0]

pf = load(File2)
reg = pf.h.region(center=centre, left_edge=domainedgeleft, right_edge=domainedgeright)
sp = pf.h.sphere(centre, 0.1)
NumberofChildlessParents = 32**3  - NumberofParents 
ExpectedParticleNumber = NumberofChildlessParents + NumberofParents*(CHILDRENPERPARENT + 1)
print "Expected Number of Particles = ", ExpectedParticleNumber
TotalNumberofParticles =  len(reg["particle_mass"])
print "Number of Particles = ", TotalNumberofParticles
ChildMass = sp["ParticleMassMsun"][0]

Delta = TotalNumberofParticles - ExpectedParticleNumber
assert(Delta == 0.0)
print "Parent Mass = %e Msun\n" % (ParentMass)
print "Child Mass = %e Msun\n" % (ChildMass)
print "Mass Difference = %.2f\n" % (ParentMass/ChildMass)
