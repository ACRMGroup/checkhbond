#ifndef ORIENTATE_H
#define ORIENTATE_H

BOOL OrientatePDB(PDB *pdb, PDB *res1_start, PDB *res1_next);
void RotateToXZ(PDB *pdb, VEC3F *n);
void RotateToX(PDB *pdb, VEC3F *n);
void RotateToXY(PDB *pdb, VEC3F *c_beta);
REAL TrueAngle(REAL opp, REAL adj);
BOOL FindNCACAtoms(PDB *res1_start, PDB *stop, VEC3F *c_alpha, VEC3F *c_beta, VEC3F *n);

#endif
