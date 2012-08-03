"""
The necessary imports for accessing Enzo data

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://enzo-project.org/
License: This file is covered under the Enzo license
"""

ctypedef double Eflt
ctypedef double FLOAT

# Now the business with the ints
ctypedef long long long_int
ctypedef long_int Eint
ctypedef int Eint32
ctypedef long_int Eint64
ctypedef long_int PINT

cdef extern from "preincludes.h":
    pass

cdef extern from "complex":
    pass

cdef extern from "ErrorExceptions.h":
    pass

cdef extern from "macros_and_parameters.h":
    enum: ENPY_INT
    enum: ENPY_PFLOAT
    enum: ENPY_BFLOAT

    enum: COMMUNICATION_POST_RECEIVE
    enum: COMMUNICATION_SEND
    enum: COMMUNICATION_RECEIVE
    enum: COMMUNICATION_SEND_RECEIVE 

    enum: DEFAULT_GHOST_ZONES
    enum: MAX_DEPTH_OF_HIERARCHY
    enum: MAX_DIMENSION
    enum: MAX_NUMBER_OF_BARYON_FIELDS
    enum: MAX_NUMBER_OF_PARTICLE_ATTRIBUTES
    enum: MAX_BARYON_FIELDS
    enum: MAX_COUNTERS
    enum: MAX_FLAGGING_METHODS
    enum: MAX_MOVIE_FIELDS
    enum: MAX_NAME_LENGTH
    enum: MAX_NUMBER_OF_NODES
    enum: MAX_NUMBER_OF_OUTPUT_REDSHIFTS
    enum: MAX_NUMBER_OF_TASKS
    enum: MAX_RECEIVE_BUFFERS
    enum: MAX_STATIC_REGIONS
    enum: MAX_TIME_ACTIONS
    
    enum: ZERO_ALL_FIELDS
    enum: ZERO_UNDER_SUBGRID_FIELD

    enum: PARTICLE_TYPE_GAS
    enum: PARTICLE_TYPE_DARK_MATTER
    enum: PARTICLE_TYPE_STAR
    enum: PARTICLE_TYPE_TRACER
    enum: PARTICLE_TYPE_MUST_REFINE
    enum: PARTICLE_TYPE_SINGLE_STAR
    enum: PARTICLE_TYPE_BLACK_HOLE
    enum: PARTICLE_TYPE_CLUSTER
    enum: PARTICLE_TYPE_MBH

    enum: NORMAL_STAR
    enum: UNIGRID_STAR
    enum: KRAVTSOV_STAR
    enum: POP3_STAR
    enum: SINK_PARTICLE
    enum: STAR_CLUSTER
    enum: INSTANT_STAR
    enum: SPRINGEL_HERNQUIST_STAR
    enum: MBH_PARTICLE

    enum: TO_DELETE
    enum: NO_FEEDBACK
    enum: ACCRETION
    enum: SUPERNOVA
    enum: CONT_SUPERNOVA
    enum: FORMATION
    enum: STROEMGREN
    enum: DEATH
    enum: MBH_THERMAL

    enum: FAIL
    enum: SUCCESS
    enum: TRUE
    enum: FALSE

cdef extern from "typedefs.h":
    pass

cdef extern from "global_data.h":
    pass

cdef extern from "Fluxes.h":
    pass

cdef extern from "GridList.h":
    pass

cdef extern from "ExternalBoundary.h":
    pass

cdef extern from "Grid.h":
    cdef cppclass grid:
        pass

cdef extern from "Hierarchy.h":
    pass

cdef extern from "communication.h":
    pass

cdef extern from "CommunicationUtilities.h":
    pass

cdef extern from "TopGridData.h":
    struct c_TopGridData "TopGridData":
        Eint   CycleNumber
        Eint   SubcycleNumber
        FLOAT Time
        double CPUTime
        double StartCPUTime
        double LastCycleCPUTime
        Eint ResubmitOn

        # Script names for resubmission to queues and restarting to reduce
        # memory fragmentation.

        char *ResubmitCommand

        # Stopping criteria for TopGrid.

        FLOAT StopTime
        Eint   StopCycle
        Eint   StopSteps
        Eflt StopCPUTime

        # Parameters governing when output is done.

        Eflt TimeLastRestartDump
        Eflt dtRestartDump

        FLOAT TimeLastDataDump
        FLOAT dtDataDump

        FLOAT TimeLastHistoryDump
        FLOAT dtHistoryDump

        FLOAT TimeLastMovieDump
        FLOAT dtMovieDump

        FLOAT TimeLastTracerParticleDump
        FLOAT dtTracerParticleDump

        FLOAT MovieRegionLeftEdge[MAX_DIMENSION]
        FLOAT MovieRegionRightEdge[MAX_DIMENSION]

        FLOAT NewMovieLeftEdge[MAX_DIMENSION]
        FLOAT NewMovieRightEdge[MAX_DIMENSION]

        Eint CycleLastRestartDump
        Eint CycleSkipRestartDump

        Eint CycleLastDataDump
        Eint CycleSkipDataDump

        Eint SubcycleLastDataDump
        Eint SubcycleSkipDataDump

        Eint CycleLastHistoryDump
        Eint CycleSkipHistoryDump
        Eint CycleSkipGlobalDataDump

        Eint OutputFirstTimeAtLevel
        Eint StopFirstTimeAtLevel

        # Parameters governing output names.

        Eint RestartDumpNumber
        Eint DataDumpNumber
        Eint HistoryDumpNumber
        Eint MovieDumpNumber
        Eint TracerParticleDumpNumber

        char *RestartDumpName
        char *DataDumpName
        char *HistoryDumpName
        char *MovieDumpName
        char *TracerParticleDumpName
        char *RedshiftDumpName

        char *RestartDumpDir
        char *DataDumpDir
        char *HistoryDumpDir
        char *MovieDumpDir
        char *TracerParticleDumpDir
        char *RedshiftDumpDir

        char *LocalDir
        char *GlobalDir

        # TopGrid Parameters governing hierarchy

        Eint StaticHierarchy

        # Some grid defining data
        # These are here out of convenience, the real ones are in the grids.

        Eint TopGridRank
        Eint TopGridDims[MAX_DIMENSION]
        #boundary_type  LeftFaceBoundaryCondition[MAX_DIMENSION],
        #              RightFaceBoundaryCondition[MAX_DIMENSION]
        char *BoundaryConditionName

        # Gravity data -- used only for top grid potential field solve

        #gravity_boundary_type GravityBoundary

        # Particle and Particle boundary data. (real one in ExternalBoundary).

        #boundary_type ParticleBoundaryType
        Eint           NumberOfParticles

        Eflt  CourantSafetyNumber
        Eint    PPMFlatteningParameter
        Eint    PPMDiffusionParameter
        Eint    PPMSteepeningParameter

        #AMRHDF5Writer AmiraGrid
        Eint FirstTimestepAfterRestart

E_ENPY_INT = ENPY_INT
E_ENPY_PFLOAT = ENPY_PFLOAT
E_ENPY_BFLOAT = ENPY_BFLOAT

E_COMMUNICATION_SEND_RECEIVE = COMMUNICATION_SEND_RECEIVE
E_COMMUNICATION_POST_RECEIVE = COMMUNICATION_POST_RECEIVE
E_COMMUNICATION_SEND = COMMUNICATION_SEND
E_COMMUNICATION_RECEIVE = COMMUNICATION_RECEIVE

E_DEFAULT_GHOST_ZONES = DEFAULT_GHOST_ZONES
E_MAX_DEPTH_OF_HIERARCHY = MAX_DEPTH_OF_HIERARCHY
E_MAX_DIMENSION = MAX_DIMENSION
E_MAX_NUMBER_OF_BARYON_FIELDS = MAX_NUMBER_OF_BARYON_FIELDS
E_MAX_NUMBER_OF_PARTICLE_ATTRIBUTES = MAX_NUMBER_OF_PARTICLE_ATTRIBUTES
E_MAX_COUNTERS = MAX_COUNTERS
E_MAX_FLAGGING_METHODS = MAX_FLAGGING_METHODS
E_MAX_MOVIE_FIELDS = MAX_MOVIE_FIELDS
E_MAX_NAME_LENGTH = MAX_NAME_LENGTH
E_MAX_NUMBER_OF_NODES = MAX_NUMBER_OF_NODES
E_MAX_NUMBER_OF_OUTPUT_REDSHIFTS = MAX_NUMBER_OF_OUTPUT_REDSHIFTS
E_MAX_NUMBER_OF_TASKS = MAX_NUMBER_OF_TASKS
E_MAX_RECEIVE_BUFFERS = MAX_RECEIVE_BUFFERS
E_MAX_STATIC_REGIONS = MAX_STATIC_REGIONS
E_MAX_TIME_ACTIONS = MAX_TIME_ACTIONS

E_ZERO_ALL_FIELDS = ZERO_ALL_FIELDS
E_ZERO_UNDER_SUBGRID_FIELD = ZERO_UNDER_SUBGRID_FIELD

E_PARTICLE_TYPE_GAS = PARTICLE_TYPE_GAS
E_PARTICLE_TYPE_DARK_MATTER = PARTICLE_TYPE_DARK_MATTER
E_PARTICLE_TYPE_STAR = PARTICLE_TYPE_STAR
E_PARTICLE_TYPE_TRACER = PARTICLE_TYPE_TRACER
E_PARTICLE_TYPE_MUST_REFINE = PARTICLE_TYPE_MUST_REFINE
E_PARTICLE_TYPE_SINGLE_STAR = PARTICLE_TYPE_SINGLE_STAR
E_PARTICLE_TYPE_BLACK_HOLE = PARTICLE_TYPE_BLACK_HOLE
E_PARTICLE_TYPE_CLUSTER = PARTICLE_TYPE_CLUSTER
E_PARTICLE_TYPE_MBH = PARTICLE_TYPE_MBH

E_NORMAL_STAR = NORMAL_STAR
E_UNIGRID_STAR = UNIGRID_STAR
E_KRAVTSOV_STAR = KRAVTSOV_STAR
E_POP3_STAR = POP3_STAR
E_SINK_PARTICLE = SINK_PARTICLE
E_STAR_CLUSTER = STAR_CLUSTER
E_INSTANT_STAR = INSTANT_STAR
E_SPRINGEL_HERNQUIST_STAR = SPRINGEL_HERNQUIST_STAR
E_MBH_PARTICLE = MBH_PARTICLE

E_TO_DELETE = TO_DELETE
E_NO_FEEDBACK = NO_FEEDBACK
E_ACCRETION = ACCRETION
E_SUPERNOVA = SUPERNOVA
E_CONT_SUPERNOVA = CONT_SUPERNOVA
E_FORMATION = FORMATION
E_STROEMGREN = STROEMGREN
E_DEATH = DEATH
E_MBH_THERMAL = MBH_THERMAL

E_FAIL = FAIL
E_SUCCESS = SUCCESS
E_TRUE = TRUE
E_FALSE = FALSE

cdef extern from "typedefs.h" nogil:

    enum: Density
    enum: TotalEnergy
    enum: InternalEnergy
    enum: Pressure
    enum: Velocity1
    enum: Velocity2
    enum: Velocity3
    enum: ElectronDensity
    enum: HIDensity
    enum: HIIDensity
    enum: HeIDensity
    enum: HeIIDensity
    enum: HeIIIDensity
    enum: HMDensity
    enum: H2IDensity
    enum: H2IIDensity
    enum: DIDensity
    enum: DIIDensity
    enum: HDIDensity
    enum: SNColour
    enum: Metallicity
    enum: ExtraType0
    enum: ExtraType1
    enum: kphHI
    enum: PhotoGamma
    enum: kphHeI
    enum: gammaHeI
    enum: kphHeII
    enum: gammaHeII
    enum: kdissH2I
    enum: GravPotential
    enum: Acceleration0
    enum: Acceleration1
    enum: Acceleration2
    enum: RadPressure0
    enum: RadPressure1
    enum: RadPressure2
    enum: Emissivity0
  
    enum: gParticlePosition
    enum: gParticleVelocity
    enum: gParticleMass
    enum: gParticleAcceleration
    enum: gParticleNumber
    enum: gParticleType
    enum: gParticleAttribute
    enum: gPotentialField
    enum: gAccelerationField
    enum: gGravitatingMassField
    enum: gFlaggingField
    enum: gVelocity
  
    enum: Bfield1
    enum: Bfield2
    enum: Bfield3
    enum: PhiField
    enum: Phi_pField
    enum: DebugField
  
    enum: DrivingField1
    enum: DrivingField2
    enum: DrivingField3
  
    enum: AccelerationField1
    enum: AccelerationField2
    enum: AccelerationField3
  
    enum: Galaxy1Colour
    enum: Galaxy2Colour
  
    enum: Mach
    enum: PreShockTemperature
    enum: PreShockDensity
    enum: CRDensity
  
    enum: CIDensity
    enum: CIIDensity
    enum: OIDensity
    enum: OIIDensity
    enum: SiIDensity
    enum: SiIIDensity
    enum: SiIIIDensity
    enum: CHIDensity
    enum: CH2IDensity
    enum: CH3IIDensity
    enum: C2IDensity
    enum: COIDensity
    enum: HCOIIDensity
    enum: OHIDensity
    enum: H2OIDensity
    enum: O2IDensity
  
    enum: MBHColour
    enum: ForbiddenRefinement
  
  
    enum: RadiationFreq0
    enum: RadiationFreq1
    enum: RadiationFreq2
    enum: RadiationFreq3
    enum: RadiationFreq4
    enum: RadiationFreq5
    enum: RadiationFreq6
    enum: RadiationFreq7
    enum: RadiationFreq8
    enum: RadiationFreq9
  
    enum: RaySegments
  
    enum: FieldUndefined

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
