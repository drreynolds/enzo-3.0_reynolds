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

cdef extern from "Hierarchy.h":
    struct c_HierarchyEntry "HierarchyEntry":
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
    enum: Dark_Matter_Density

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
