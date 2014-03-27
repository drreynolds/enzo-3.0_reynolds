"""
The necessary imports for accessing Enzo data

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://enzo-project.org/
License: This file is covered under the Enzo license
"""

cdef int E_ENPY_INT = ENPY_INT
cdef int E_ENPY_PFLOAT = ENPY_PFLOAT
cdef int E_ENPY_BFLOAT = ENPY_BFLOAT

cdef int E_COMMUNICATION_SEND_RECEIVE = COMMUNICATION_SEND_RECEIVE
cdef int E_COMMUNICATION_POST_RECEIVE = COMMUNICATION_POST_RECEIVE
cdef int E_COMMUNICATION_SEND = COMMUNICATION_SEND
cdef int E_COMMUNICATION_RECEIVE = COMMUNICATION_RECEIVE

cdef int E_MAX_DEPTH_OF_HIERARCHY = MAX_DEPTH_OF_HIERARCHY
cdef int E_MAX_DIMENSION = MAX_DIMENSION
cdef int E_MAX_NUMBER_OF_BARYON_FIELDS = MAX_NUMBER_OF_BARYON_FIELDS
cdef int E_MAX_NUMBER_OF_PARTICLE_ATTRIBUTES = MAX_NUMBER_OF_PARTICLE_ATTRIBUTES
cdef int E_MAX_COUNTERS = MAX_COUNTERS
cdef int E_MAX_FLAGGING_METHODS = MAX_FLAGGING_METHODS
cdef int E_MAX_MOVIE_FIELDS = MAX_MOVIE_FIELDS
cdef int E_MAX_NAME_LENGTH = MAX_NAME_LENGTH
cdef int E_MAX_NUMBER_OF_NODES = MAX_NUMBER_OF_NODES
cdef int E_MAX_NUMBER_OF_OUTPUT_REDSHIFTS = MAX_NUMBER_OF_OUTPUT_REDSHIFTS
cdef int E_MAX_NUMBER_OF_TASKS = MAX_NUMBER_OF_TASKS
cdef int E_MAX_RECEIVE_BUFFERS = MAX_RECEIVE_BUFFERS
cdef int E_MAX_STATIC_REGIONS = MAX_STATIC_REGIONS
cdef int E_MAX_TIME_ACTIONS = MAX_TIME_ACTIONS

cdef int E_ZERO_ALL_FIELDS = ZERO_ALL_FIELDS
cdef int E_ZERO_UNDER_SUBGRID_FIELD = ZERO_UNDER_SUBGRID_FIELD

cdef int E_PARTICLE_TYPE_GAS = PARTICLE_TYPE_GAS
cdef int E_PARTICLE_TYPE_DARK_MATTER = PARTICLE_TYPE_DARK_MATTER
cdef int E_PARTICLE_TYPE_STAR = PARTICLE_TYPE_STAR
cdef int E_PARTICLE_TYPE_TRACER = PARTICLE_TYPE_TRACER
cdef int E_PARTICLE_TYPE_MUST_REFINE = PARTICLE_TYPE_MUST_REFINE
cdef int E_PARTICLE_TYPE_SINGLE_STAR = PARTICLE_TYPE_SINGLE_STAR
cdef int E_PARTICLE_TYPE_BLACK_HOLE = PARTICLE_TYPE_BLACK_HOLE
cdef int E_PARTICLE_TYPE_CLUSTER = PARTICLE_TYPE_CLUSTER
cdef int E_PARTICLE_TYPE_MBH = PARTICLE_TYPE_MBH

cdef int E_NORMAL_STAR = NORMAL_STAR
cdef int E_UNIGRID_STAR = UNIGRID_STAR
cdef int E_KRAVTSOV_STAR = KRAVTSOV_STAR
cdef int E_POP3_STAR = POP3_STAR
cdef int E_SINK_PARTICLE = SINK_PARTICLE
cdef int E_STAR_CLUSTER = STAR_CLUSTER
cdef int E_INSTANT_STAR = INSTANT_STAR
cdef int E_SPRINGEL_HERNQUIST_STAR = SPRINGEL_HERNQUIST_STAR
cdef int E_MBH_PARTICLE = MBH_PARTICLE

cdef int E_TO_DELETE = TO_DELETE
cdef int E_NO_FEEDBACK = NO_FEEDBACK
cdef int E_ACCRETION = ACCRETION
cdef int E_SUPERNOVA = SUPERNOVA
cdef int E_CONT_SUPERNOVA = CONT_SUPERNOVA
cdef int E_FORMATION = FORMATION
cdef int E_STROEMGREN = STROEMGREN
cdef int E_DEATH = DEATH
cdef int E_MBH_THERMAL = MBH_THERMAL

cdef int E_FAIL = FAIL
cdef int E_SUCCESS = SUCCESS
cdef int E_TRUE = TRUE
cdef int E_FALSE = FALSE

cdef field_enums = {
    "Density" : Density,
    "TotalEnergy" : TotalEnergy,
    "GasEnergy" : InternalEnergy,
    "Pressure" : Pressure,
    "x-velocity" : Velocity1,
    "y-velocity" : Velocity2,
    "z-velocity" : Velocity3,
    "Electron_Density" : ElectronDensity,
    "HI_Density" : HIDensity,
    "HII_Density" : HIIDensity,
    "HeI_Density" : HeIDensity,
    "HeII_Density" : HeIIDensity,
    "HeIII_Density" : HeIIIDensity,
    "HM_Density" : HMDensity,
    "H2I_Density" : H2IDensity,
    "H2II_Density" : H2IIDensity,
    "DI_Density" : DIDensity,
    "DII_Density" : DIIDensity,
    "HDI_Density" : HDIDensity,
    "SNColour" : SNColour,
    "Metallicity" : Metallicity,
    "ExtraType0" : ExtraType0,
    "ExtraType1" : ExtraType1,
    "kphHI" : kphHI,
    "PhotoGamma" : PhotoGamma,
    "kphHeI" : kphHeI,
    "gammaHeI" : gammaHeI,
    "kphHeII" : kphHeII,
    "gammaHeII" : gammaHeII,
    "kdissH2I" : kdissH2I,
    "GravPotential" : GravPotential,
    "Acceleration0" : Acceleration0,
    "Acceleration1" : Acceleration1,
    "Acceleration2" : Acceleration2,
    "RadPressure0" : RadPressure0,
    "RadPressure1" : RadPressure1,
    "RadPressure2" : RadPressure2,
    "Emissivity0" : Emissivity0,
    "Dark_Matter_Density" : Dark_Matter_Density,
  
    "gParticlePosition" : gParticlePosition,
    "gParticleVelocity" : gParticleVelocity,
    "gParticleMass" : gParticleMass,
    "gParticleAcceleration" : gParticleAcceleration,
    "gParticleNumber" : gParticleNumber,
    "gParticleType" : gParticleType,
    "gParticleAttribute" : gParticleAttribute,
    "gPotentialField" : gPotentialField,
    "gAccelerationField" : gAccelerationField,
    "gGravitatingMassField" : gGravitatingMassField,
    "gFlaggingField" : gFlaggingField,
    "gVelocity" : gVelocity,
  
    "Bfield1" : Bfield1,
    "Bfield2" : Bfield2,
    "Bfield3" : Bfield3,
    "PhiField" : PhiField,
    "Phi_pField" : Phi_pField,
    "DebugField" : DebugField,
  
    "DrivingField1" : DrivingField1,
    "DrivingField2" : DrivingField2,
    "DrivingField3" : DrivingField3,
  
    "AccelerationField1" : AccelerationField1,
    "AccelerationField2" : AccelerationField2,
    "AccelerationField3" : AccelerationField3,
  
    "Galaxy1Colour" : Galaxy1Colour,
    "Galaxy2Colour" : Galaxy2Colour,
  
    "Mach" : Mach,
    "PreShockTemperature" : PreShockTemperature,
    "PreShockDensity" : PreShockDensity,
    "CRDensity" : CRDensity,
  
    "CI_Density" : CIDensity,
    "CII_Density" : CIIDensity,
    "OI_Density" : OIDensity,
    "OII_Density" : OIIDensity,
    "SiI_Density" : SiIDensity,
    "SiII_Density" : SiIIDensity,
    "SiIII_Density" : SiIIIDensity,
    "CHI_Density" : CHIDensity,
    "CH2I_Density" : CH2IDensity,
    "CH3II_Density" : CH3IIDensity,
    "C2I_Density" : C2IDensity,
    "COI_Density" : COIDensity,
    "HCOII_Density" : HCOIIDensity,
    "OHI_Density" : OHIDensity,
    "H2OI_Density" : H2OIDensity,
    "O2I_Density" : O2IDensity,
  
    "MBHColour" : MBHColour,
    "ForbiddenRefinement" : ForbiddenRefinement,
  
  
    "RadiationFreq0" : RadiationFreq0,
    "RadiationFreq1" : RadiationFreq1,
    "RadiationFreq2" : RadiationFreq2,
    "RadiationFreq3" : RadiationFreq3,
    "RadiationFreq4" : RadiationFreq4,
    "RadiationFreq5" : RadiationFreq5,
    "RadiationFreq6" : RadiationFreq6,
    "RadiationFreq7" : RadiationFreq7,
    "RadiationFreq8" : RadiationFreq8,
    "RadiationFreq9" : RadiationFreq9,
  
    "RaySegments" : RaySegments,
  
    "FieldUndefined" : FieldUndefined
}

cdef extern from "fix_enzo_defs.h":
    pass
