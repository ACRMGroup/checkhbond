#ifndef ORIENTATE_H
#define ORIENTATE_H

BOOL OrientatePDB(PDB *pdb, PDB *res1_start, PDB *res1_next);
BOOL OrientateN_PDB(PDB *pdb, PDB *prev, PDB *start, PDB *next);
BOOL OrientateCO_PDB(PDB *pdb, PDB *start, PDB *next);
void RotateToXZ(PDB *pdb, VEC3F *n);
void RotateToX(PDB *pdb, VEC3F *n);
void RotateToXY(PDB *pdb, VEC3F *c_beta);
REAL TrueAngle(REAL opp, REAL adj);
BOOL FindNCACBAtoms(PDB *res1_start, PDB *stop, VEC3F *c_alpha, VEC3F *c_beta, VEC3F *n);
BOOL FindCACOAtoms(PDB *res1_start, PDB *stop, VEC3F *c_alpha, VEC3F *c_beta, VEC3F *n);
BOOL FindCNCAAtoms(PDB *res0_start, PDB *res1_start, PDB *stop, VEC3F *c, VEC3F *n, VEC3F *c_alpha);
BOOL ResiduesBonded(PDB *pdb1, PDB *pdb2);
void PrintError(FILE *out, char *text);

/* Max distance for C...N peptide bond  */
#define CNBONDSQ 2.25    /* 1.5A */

#endif
