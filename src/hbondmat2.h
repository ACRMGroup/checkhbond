#define DIV_PER_ANGSTROM 2
#define MAXCAHBDIST      15

#define MAXSIZE (2*MAXCAHBDIST*DIV_PER_ANGSTROM)
#define DIV ((REAL)1/DIV_PER_ANGSTROM)
#define OFFSET (MAXSIZE/2)
#define MAXBUFF 300
#define PDBEXT ".ent"
#define DOMAINLOC "/acrm/data/dompdb/"
#define PDBLOC "/data/pdb/pdb"
#define PDBSTART ""
#define PGPFILE "Explicit.pgp"
#define COORD_2_GRID(x,y) (x) = ((int)(y/DIV) + OFFSET)
#define GRID_2_COORD(x,y) (y) = ((REAL)(x - OFFSET) * DIV)
#define VALIDGRIDCOORDS(x, y, z)  \
 (((x) < MAXSIZE) && ((x) >= 0) && \
  ((y) < MAXSIZE) && ((y) >= 0) &&  \
  ((z) < MAXSIZE) && ((z) >= 0))
#define GRIDSPACING ((REAL)1/DIV_PER_ANGSTROM)
#define RESOL 2.5
#define DEFAULT_CUTOFF_VALUE 0.25 /* default cutoff for off-grid matching */

