class EventDataContainer {

};

class EvolveLevelEventDataContainer : public EventDataContainer {
public:
    int level; 
    int grid_ind; 
    SiblingGridList *SiblingList; 
    int NumberOfGrids; 
    LevelHierarchyEntry **LevelArray; 
    ExternalBoundary *Exterior; 
    fluxes ***SubgridFluxesEstimate; 
    int *NumberOfSubgrids; 
    Star *AllStars; 
    int *TotalStarParticleCountPrevious; 
    float dtLevelAbove; 
    void *ImplicitSolver; 
    static float TopGridTimeStep;
    static float norm;
    static int StaticLevelZero;
};

class TestEventDataContainer : public EventDataContainer {
public:
    int SomethingToPrint;
};
