#ifndef LATTICE_TYPE
#define LATTICE_TYPE

typedef struct lattice_type{
  char comment[100];
  
  double ax, ay, az;		/* lattice vectors */
  double bx, by, bz;
  double cx, cy, cz;

  int tot;			/* total number of lattice points */
  double *x, *y, *z;		/* directions */
  double *charge;		/* charge */
  double *moments;		/* local moments */
  double *vx, *vy, *vz;		/* vibration modes */
  
  int species;			/* atomic info (POSCAR)*/
  char atom_sym[100][4];
  int atom_num[100];
  double atom_radius[100];
  double atom_mass[100];
} lattice_t;

#endif

/*
  Read lattice information from a XYZ file or a POSCAR file.
 */
int LoadLatticeFromXYZ(char * fname, lattice_t *lattice);
int LoadLatticeFromPOSCAR(char *fname, lattice_t *lattice);
int LoadLatticeFromXDATCAR(char *fname, lattice_t *lattice);

/*
  Print Latttice information to standard output;
  PrintLattice: print one lattice;
  Print2Lattices: print 2 lattices;
 */
int PrintLattice(lattice_t *lattice);
int print2Lattices(lattice_t *lattice1, lattice_t *lattice2);

int CopyLattice(lattice_t *dest, lattice_t *source);
int FreeLattice(lattice_t *lattice);

int LoadRWIGS(lattice_t *lattice);
int LoadAtomMass(lattice_t * lattice);
