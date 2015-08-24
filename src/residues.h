#ifndef RESIDUES_H
#define RESIDUES_H

void WriteResidue(FILE *fp, PDB  *start, PDB *stop);
void WriteResidueRecord(FILE *fp, PDB *pdb);
PDB *SelectAtomsResidue(PDB *start, PDB *stop, int nsel, char **sel, int *natom);
void ApplyMatrixResidue(PDB  *start, PDB *stop, REAL matrix[3][3]);
void TranslateResidue(PDB *start, PDB *stop, VEC3F tvect);

#endif
