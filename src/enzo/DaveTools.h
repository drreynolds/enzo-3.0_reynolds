void Pout( char * string, int i1 = -12345, int i2 = -12345,
           int i3 = -12345, int i4 = -12345, int i5 = -12345);

void PoutF( char * string, float f1 = -12345.6,float f2 = -12345.6,float f3 = -12345.6,
	    float f4 = -12345.6,float f5 = -12345.6,float f6 = -12345.6 );

void dump(float *A, int nx, int ny, int nz, int nb, char * filename);
//<dcc>
void WriteSingleCube(float * array, int Dims[], char* string, int dNum, int gNum, char * label);
//</dcc>

int TracerParticlesAddToRestart_DoIt(char * filename, HierarchyEntry *TopGrid,
				    TopGridData *MetaData);

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif


double wall_time (char * string);
EXTERN FILE * wall_ptr;
void wall_time_start();
void wall_time_stop();
void wall_time_flush();

float wall_time();
