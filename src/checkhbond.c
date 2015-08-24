/*************************************************************************

   Program:    checkhbond
   File:       checkhbond.c
   
   Version:    V2.0
   Date:       24.01.06
   Function:   Generate matrices of hydrogen bond information for use
               by checkhbond
   
   Copyright:  (c) University of Reading / Alison L. Cuff 2002-2006
   Author:     Alison L. Cuff
   Address:    School of Animal and Microbial Sciences,
               The University of Reading,
               Whiteknights,
               P.O. Box 228,
               Reading RG6 6AJ.
               England.
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   If a hydrogen-capable residue is substituted for another will the 
   hydrogen-bond between the replacement residue and the partner residue 
   of the original be maintained?

   Steps:

   Orientate protein so that res 1 has its CA at the origin.
   Cull matrices to ensure that there are no steric clashes vector is 
      calculated from CA of res 1 to CA of res 2
   res 2 is moved with CA at the origin
   matfit is then used to obtain a rotation matrix (with res 2 as the 
      'mobile' residue)

   Run though grid for res 2 for all non-zero positions
   for each position:
      convert to real x, y, z co-ordinates
      multiply the resulting vector by the rotation matrix -- res 2 grid 
         points should be at same orientation now as res 1 grid points
      move the vector by the CatoCa vector so that both grids are with 
         CA at the origin
      Convert real co-ordinates back to grid locations and match with 
         locations in other grid

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  17.06.03 Original version as used for thesis
   V1.1  19.08.05 Various bug fixes By: ACRM
   V2.0  24.01.05 Modified to allow mc/sc matrices to be used

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/matrix.h"
#include "bioplib/hbond.h"
#include "bioplib/fit.h"
#include "residues.h"
#include "orientate.h"
#include "hbondmat2.h"

/************************************************************************/
/* Defines and macros
*/
/* default matrix files */
#define MATRIXFILE            "/acrm/home/alison/hydrogen_bonding/matrices_05new.txt"
#define MATRIXFILE_MCDONOR    "/acrm/home/alison/hydrogen_bonding/matrices_05new.txt"
#define MATRIXFILE_MCACCEPTOR "/acrm/home/alison/hydrogen_bonding/matrices_05new.txt"
/* radius of atom */
#define RAD 25
/* calculates distance in angstroms between atoms */
#define DISTSQ_HATOMS(a,b) (a.x - b.x) * (a.x - b.x) + \
                       (a.y - b.y) * (a.y - b.y) + \
                       (a.z - b.z) * (a.z - b.z)
/* Error codes for GetResidues() */
#define ERR_NOMEM      0
#define ERR_NOPREVRES1 1
#define ERR_NOPREVRES2 2

/* Matrix reading styles for ReadInMatrices */
#define MAT_READ_DONOR1    1
#define MAT_READ_DONOR2    2
#define MAT_READ_ACCEPTOR1 4
#define MAT_READ_ACCEPTOR2 8
#define MAT_READ_BOTH      9
#define MAT_RES_1          1
#define MAT_RES_2          2
#define MAT_RES_BOTH       3

/* Atom sets for CreateRotationMatrix */
#define ATOMS_NCAC   1
#define ATOMS_CNCA   2
#define ATOMS_CACO   3
#define ATOMS_NCACB  4

/* Max distance for C...N peptide bond  */
#define CNBONDSQ 2.25    /* 1.5A */

/* for debugging purposes */

/*
#define DEBUG
#define DEBUG2
#define DEBUG3
*/

/************************************************************************/
/* Globals
*/
/* matrices that store all acceptor and donor atoms */
int gDonate[MAXSIZE][MAXSIZE][MAXSIZE];
int gAccept[MAXSIZE][MAXSIZE][MAXSIZE];
/* matrix storing partner atoms to hydrogen accepting atoms */
int gPartnertoAccept[MAXSIZE][MAXSIZE][MAXSIZE];
int gPartnertoDonate[MAXSIZE][MAXSIZE][MAXSIZE];
/* rotation matrix */
REAL gRotation_matrix[3][3];
#if defined(DEBUG1) || defined(DEBUG2)
char gChain = 'Z';
#endif
   
/************************************************************************/
/* Prototypes
*/
int main (int argc, char *argv[]);
void Usage(void);
BOOL ReadInMatrices(char *res1, char *res2, FILE *matrix, int type, int whichres);
void CullArrays(PDB *pdb, PDB *res1, PDB *res2,
                int donate_array[MAXSIZE][MAXSIZE][MAXSIZE],
                int accept_array[MAXSIZE][MAXSIZE][MAXSIZE]);
void CalculateCaToCaVector(PDB *res1_start, PDB *res1_stop,
                           PDB *res2_start, PDB *res2_stop,
                           VEC3F *CAtoCAVector);
PDB *FindAtom(PDB *start, PDB *stop, char *atnam, VEC3F *c_alpha);
BOOL CreateRotationMatrix(PDB *pdb, PDB *res1_start,
                          PDB *res1_stop, PDB *res2_start,
                          PDB *res2_stop, VEC3F CAtoCAVector, int atomset1, int atomset2, PDB *prevres1);
BOOL CheckValidHBond(VEC3F CAtoCAVector, REAL cutoff,
                     FILE *out,
                     int keyarray[MAXSIZE][MAXSIZE][MAXSIZE],
                     int partnerarray[MAXSIZE][MAXSIZE][MAXSIZE]);
void ClearArrays();
BOOL ParseCmdLine(int argc, char **argv, REAL *cutoff,
                  BOOL *hbplus, char *hatom1,
                  char *hatom2, char *matrix_file, char *matrix_file2, char *pdbfile,
                  char *locres1, char *locres2, char *res2,
                  char *outputfile);
BOOL CalculateHBondEnergy(PDB *pdb, PDB *res1_start, PDB *res1_stop,
                          PDB *res2_start, PDB *res2_stop, VEC3F CAtoCAVector,
                          REAL cutoff,
                          char *hatom1, char *hatom2, FILE *out);
void OrientateMatrix(VEC3F CAtoCAVector,int x, int y, int z,
                     VEC3F *rotated_coord);
int CalculateTotalCounts(int array[MAXSIZE][MAXSIZE][MAXSIZE]);
REAL CalcEnergy(int count1, int totalcount1, int count2, int totalcount2);
REAL Distance_squared (int i, int j, int k, int x_coord, int y_coord,
                       int z_coord);
REAL DoCheckHBond(int x, int y, int z, VEC3F CAtoCAVector, VEC3F *partner_coord,
                  int totalcount1, int totalcount2,
                  int keyarray[MAXSIZE][MAXSIZE][MAXSIZE],
                  int partnerarray[MAXSIZE][MAXSIZE][MAXSIZE],
                  REAL cutoff, FILE *out);
FILE *OpenMatrixFile(char *matrix_file, char *def_matrix_file);
BOOL Open_Std_Files(char *infile, char *outfile, FILE **in, FILE **out);
BOOL PrepareHBondingPair(int resnum1, int resnum2, PDB *pdb, FILE *matrix, char *chain1, 
                         char *chain2, char *insert1, char *insert2, BOOL *hbplus,
                         char *hatom1, char *hatom2, REAL cutoff, char *res1, char *res2,
                         FILE *OUT);
PDB *GetResidues(PDB *pdb, char *chain1, int resnum1, char *insert1, 
                 char *chain2, int resnum2, char *insert2, int *errorcode);
void FindRes1Type(PDB *pdb, char *chain, int resnum, char *insert, char *res);
BOOL ResiduesBonded(PDB *pdb1, PDB *pdb2);
BOOL AnalyzeMCDonorPair(int resnum1, int resnum2, PDB *pdb, FILE *matrix, FILE *matrix2,
                        char *chain1, 
                        char *chain2, char *insert1, char *insert2, BOOL *hbplus,
                        char *hatom1, char *hatom2, REAL cutoff, char *res1, char *res2,
                        FILE *OUT);
void CalculateNToCaVector(PDB *res1_start, PDB *res1_stop,
                          PDB *res2_start, PDB *res2_stop, VEC3F *NtoCAVector);
BOOL AnalyzeMCAcceptorPair(int resnum1, int resnum2, PDB *pdb, FILE *matrix, FILE *matrix2,
                           char *chain1, 
                           char *chain2, char *insert1, char *insert2, BOOL *hbplus,
                           char *hatom1, char *hatom2, REAL cutoff, char *res1, char *res2,
                           FILE *OUT);
void CalculateCToCaVector(PDB *res1_start, PDB *res1_stop,
                          PDB *res2_start, PDB *res2_stop, VEC3F *CtoCAVector);
PDB *SwapCarbon(PDB *pdb, PDB *prev);


/************************************************************************/
int main (int argc, char *argv[])
{
   FILE *PDBFILE = stdin,
      *OUT = stdout,
      *matrix = NULL;

   char locres1[6], locres2[6],res1[4], res2[6],
      pdbfile[MAXBUFF], outputfile[MAXBUFF];
   char chain1[6], insert1[6], chain2[6], insert2[6], hatom1[6], hatom2[6];
   char matrix_file[MAXBUFF];
   char matrix_file2[MAXBUFF];
   PDB *pdb;
   int natoms, resnum1, resnum2, errorcode;
   REAL cutoff;
   BOOL hbplus = FALSE;

#if defined(MCDONOR) || defined(MCACCEPTOR)
   FILE *matrix2 = NULL;
#endif

   
   /* set array elements to 0 */
   ClearArrays();
   
   if(ParseCmdLine(argc, argv,  &cutoff, &hbplus,  hatom1, hatom2,
                   matrix_file, matrix_file2, pdbfile, locres1, locres2, res2, outputfile))
   {
      if((matrix = OpenMatrixFile(matrix_file, MATRIXFILE)))
      {
         if(ParseResSpec(locres1, chain1, &resnum1, insert1))
         {
            if(ParseResSpec(locres2, chain2, &resnum2, insert2))
            {               
               if(Open_Std_Files(pdbfile, outputfile, &PDBFILE, &OUT))
               { 
                  /* create linked list of pdb file */ 
                  if((pdb = ReadPDB(PDBFILE, &natoms)) !=NULL)
                  {
                     /* ACRM 08.09.05 Get only the residues of interest */
                     if((pdb = GetResidues(pdb, chain1, resnum1, insert1, 
                                           chain2, resnum2, insert2,
                                           &errorcode))==NULL)
                     {
                        if(errorcode == ERR_NOMEM)
                        {
                           fprintf(stderr,"No memory for storing residues of interest\n");
                           return(1);
                        }
                        else if(errorcode == ERR_NOPREVRES1)
                        {
                           fprintf(stderr,"No preceeding residue for residue %c%d%c\n",
                                   chain1[0], resnum1, insert1[0]);
                           return(1);
                        }
                        else if(errorcode == ERR_NOPREVRES2)
                        {
                           fprintf(stderr,"No preceeding residue for residue %c%d%c\n",
                                   chain2[0], resnum2, insert2[0]);
                           return(1);
                        }
                        else
                        {
                           fprintf(stderr,"Undefined error in getting residues\n");
                           return(1);
                        }
                     }
                     FindRes1Type(pdb, chain1, resnum1, insert1, res1);


#ifdef MCDONOR
                     if((matrix2 = OpenMatrixFile(matrix_file2, MATRIXFILE_MCDONOR)))
                     {
                        if(!AnalyzeMCDonorPair(resnum1, resnum2, pdb, matrix, matrix2, chain1, chain2, 
                                               insert1, insert2, &hbplus, hatom1, hatom2,
                                               cutoff, res1, res2, OUT))
                        {
                           fprintf(stderr, "No hydrogen bonds\n");
                        }
                     }
#elif  MCACCEPTOR
                     if((matrix2 = OpenMatrixFile(matrix_file2, MATRIXFILE_MCDONOR)))
                     {
                        if(!AnalyzeMCAcceptorPair(resnum1, resnum2, pdb, matrix, matrix2, chain1, chain2, 
                                                  insert1, insert2, &hbplus, hatom1, hatom2,
                                                  cutoff, res1, res2, OUT))
                        {
                           fprintf(stderr, "No hydrogen bonds\n");
                        }
                     }
#else
                     if(!PrepareHBondingPair(resnum1, resnum2, pdb, matrix, chain1, chain2, 
                                             insert1, insert2, &hbplus, hatom1, hatom2,
                                             cutoff, res1, res2, OUT))
                     {
                        /* ACRM 08.09.05 Swap chain, inserts and hatom as well! */
                        if(!PrepareHBondingPair(resnum2, resnum1, pdb, matrix, chain2, chain1, 
                                                insert2, insert1, &hbplus, hatom2, hatom1,
                                                cutoff, res2, res1, OUT))
                        {
                           fprintf(stderr, "No hydrogen bonds\n");
                        }
                     }
#endif
                  }
                  else
                  {
                     fprintf(stderr, "Sorry, cannot read pdb file\n");   
                     return(1);
                  }
               }
               else
               {
                  fprintf(stderr, "Sorry, unable to open PDB file or output file\n");
                  return(1);
               }
            }
            else
            {
               fprintf(stderr, "Sorry, unable to parse second residue correctly\n");
               return(1);
            }
         }
         else
         {
            fprintf(stderr, "Sorry, unable to parse first residue correctly\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr, "Sorry, unable to open matrix file\n");
      }
   }
   else
   {
      Usage();
   }
  
   return(0);
}


/************************************************************************/
BOOL PrepareHBondingPair(int resnum1, int resnum2, PDB *pdb, FILE *matrix,
                         char *chain1, char *chain2, 
                         char *insert1, char *insert2, BOOL *hbplus,
                         char *hatom1, char *hatom2, REAL cutoff, 
                         char *res1, char *res2, FILE *OUT)
{
   VEC3F CAtoCAVector;

   PDB *res1_start, *res1_stop, *res2_start, *res2_stop;
   
   for(res1_start = pdb; res1_start !=NULL;
       res1_start=res1_stop)
   {
      res1_stop = FindNextResidue(res1_start);
                        
      /* if *key* residue of choice */
      if((resnum1 == res1_start->resnum)
         && (chain1[0]==res1_start->chain[0])
         && (insert1[0] == res1_start->insert[0]))
      {
         break;
      }
   }
   if(res1_start==NULL)
   {
      fprintf(stderr, "Error: Can't find residue 1  in PDB file\n");
      return(FALSE);
   }
   
   /*step though each residue in list
     is it *partner* residue? */
   for(res2_start = pdb; res2_start !=NULL;
       res2_start=res2_stop)
   {
      res2_stop = FindNextResidue(res2_start);
      
      /* if *partner* residue of choice */
      if((resnum2 == res2_start->resnum)
         && (chain1[0]==res1_start->chain[0])
         && (insert1[0] == res1_start->insert[0]))
      {
         break;
      }
   }
   if(res2_start==NULL)
   {
      fprintf(stderr, "Error: Can't find residue 2 in PDB file\n");
      return(FALSE);
   }
   
   if(!ReadInMatrices(res1, res2, matrix, MAT_READ_BOTH, MAT_RES_BOTH))
   {
      /* ACRM 08.09.05 - error message */
/*      fprintf(stderr, "Unable to find one of the residues in the matrix file: %s %s\n", res1, res2); */
      return(FALSE);
   }
   
   OrientatePDB(pdb, res2_start, res2_stop);
#ifndef NOCULL
   CullArrays(pdb, res1_start, res2_start, gPartnertoDonate,
              gPartnertoAccept);
#endif
   OrientatePDB(pdb, res1_start, res1_stop);
#ifndef NOCULL
   CullArrays(pdb, res1_start, res2_start,
              gDonate, gAccept);
#endif   
   CalculateCaToCaVector(res1_start, res1_stop, res2_start,
                         res2_stop, &CAtoCAVector);
   
   /* 19.01.06 Now fits on N,CA,CB rather than N,CA,C for consistency with
      the frame of reference
   */
   CreateRotationMatrix(pdb, res1_start, res1_stop, res2_start,
                        res2_stop, CAtoCAVector, ATOMS_NCACB, ATOMS_NCACB, NULL);

   if(*hbplus == TRUE)    /* ACRM 23.07.04 Dereference pointer! */
   {
      OrientatePDB(pdb, res2_start, res2_stop);
      if(!CalculateHBondEnergy(pdb, res1_start, res1_stop,
                               res2_start,
                               res2_stop, CAtoCAVector, cutoff,
                               hatom1, hatom2, OUT))
      {
         fprintf(stderr, "No hydrogen bonds\n");
      }
   }
   else
   {
      if(!CheckValidHBond(CAtoCAVector, cutoff, OUT,
                          gDonate, gPartnertoAccept))
      {
         if(!CheckValidHBond(CAtoCAVector, cutoff, OUT,
                             gAccept, gPartnertoDonate))
         {
            return(FALSE);
         }
      }
   }
   return(TRUE);
}

/************************************************************************/
BOOL CalculateHBondEnergy(PDB *pdb, PDB *res1_start, PDB *res1_stop,
                          PDB *res2_start,PDB *res2_stop,
                          VEC3F CAtoCAVector, REAL cutoff,
                          char *hatom1, char *hatom2, FILE *out)
{
   int x_coord1, y_coord1, z_coord1,
      totalcount1 = 0, 
      totalcount2 = 0;
   FILE *OUT = out;
   VEC3F partner_coord;
   REAL pseudoenergy = -1;
   PDB *p;
   BOOL OK = FALSE;
      
   /* calculate total number of hydrogen bonding atoms for residue 1 and 2 */
   totalcount1 = CalculateTotalCounts(gDonate);
   totalcount2 = CalculateTotalCounts(gPartnertoAccept);
   
   /* obtaining x,y,z co-ordinates for donor atom */
   for(p = res1_start; p !=res1_stop; NEXT(p))
   {
      if(strstr(p->atnam, hatom1))
      {
         COORD_2_GRID(x_coord1,p->x);
         COORD_2_GRID(y_coord1,p->y);
         COORD_2_GRID(z_coord1,p->z);
         
         OK = TRUE;
         break;
      }
   }

   if(!OK)
   {
      fprintf(stderr, "Can't find %s in protein structure\n", hatom1);
      return(FALSE);
   }
      
   pseudoenergy = DoCheckHBond(x_coord1, y_coord1, z_coord1, CAtoCAVector,&partner_coord,
                               totalcount1, totalcount2,
                               gDonate, gPartnertoAccept, cutoff, OUT); 
   if((pseudoenergy !=-1)&& (pseudoenergy !=9999.9999))
   {
      fprintf(out, "Pseudoenergy of best quality hydrogen bond: %.2f\n",
              pseudoenergy);
      printf("Pseudoenergy of best quality hydrogen bond: %.2f\n",
             pseudoenergy);
      return(TRUE);
   }
   else
   {
      fprintf(out, "No hydrogen bonds\n");
      printf("No hydrogen bonds\n");
      return(FALSE);
   }
}

/************************************************************************/
int CalculateTotalCounts(int array[MAXSIZE][MAXSIZE][MAXSIZE])
{
   int totalcount = 0, x, y, z;
   
   for(x=0; x <MAXSIZE; x++)
   {
      for(y=0; y <MAXSIZE; y++)
      {
         for(z=0; z <MAXSIZE; z++)
         {
            totalcount = totalcount + array[x][y][z];
         }
      }
   }
   return(totalcount);
}

/************************************************************************/
/* function that calculates the vector from the CA of res1 to CA of res2 
*/
void CalculateCaToCaVector(PDB *res1_start, PDB *res1_stop,
                           PDB *res2_start, PDB *res2_stop, 
                           VEC3F *CAtoCAVector)
{
   VEC3F res1_calpha, 
      res2_calpha;

   /*find co-ordinates of CA atom */   
   if(!FindAtom(res1_start, res1_stop, "CA  ", &res1_calpha))
   {
      fprintf(stderr, "Error: Can't find c-alpha atoms of key residue\n");
   }
   
   if(!FindAtom(res2_start, res2_stop, "CA  ", &res2_calpha))
   {
      fprintf(stderr, "Error: Can't find c-alpha atoms of partner residue\n");
   }
   
   CAtoCAVector->x =  res2_calpha.x - res1_calpha.x;
   CAtoCAVector->y =  res2_calpha.y - res1_calpha.y;
   CAtoCAVector->z =  res2_calpha.z - res1_calpha.z;
}

/************************************************************************/
/* function that creates a rotation matrix (global) to fit res 1
   (*key* residue) onto res 2 (*partner* residue). A weight of 1.0
   is added to atoms CA and N and 0.1 to atom C. 
   Includes option to print out residues at set stages for debugging
   purposes.

   If we are ever to do backbone-backbone we will need to support
   atomset2 being C,N,CA
*/
BOOL CreateRotationMatrix(PDB *pdb, PDB *res1_start, PDB *res1_stop,
                          PDB *res2_start, PDB *res2_stop,VEC3F Vector,
                          int atomset1, int atomset2, PDB *prevres1)
{
   VEC3F tempv;   
   
   PDB  *keyres1_pdb = NULL, 
        *partnerres2_pdb = NULL;
   char *sel1[4], *sel2[4];
   
   REAL *weight = NULL;
   int natoms,
       NumCo_ord1 = 0,
       NumCo_ord2 = 0,
       natom      = 0,
       resnum     = 0,
       i          = 0;
   char insert = ' ';
   PDB *p         = NULL,
       *q         = NULL,
       *start     = NULL;
   COOR *keyres1_coor = NULL,
      *partnerres2_coor = NULL;
    
   switch(atomset1)
   {
   case ATOMS_NCAC:
      SELECT(sel1[0], "N   ");
      SELECT(sel1[1], "CA  ");
      SELECT(sel1[2], "C   ");
      break;
   case ATOMS_NCACB:
      SELECT(sel1[0], "N   ");
      SELECT(sel1[1], "CA  ");
      SELECT(sel1[2], "CB  ");
      break;
   case ATOMS_CNCA:
      SELECT(sel1[0], "C   ");
      SELECT(sel1[1], "N   ");
      SELECT(sel1[2], "CA  ");
      break;
   case ATOMS_CACO:
      SELECT(sel1[0], "CA  ");
      SELECT(sel1[1], "C   ");
      SELECT(sel1[2], "O   ");
      break;
   default:
      return(FALSE);
   }
   
   switch(atomset2)
   {
   case ATOMS_NCAC:
      SELECT(sel2[0], "N   ");
      SELECT(sel2[1], "CA  ");
      SELECT(sel2[2], "C   ");
      break;
   case ATOMS_NCACB:
      SELECT(sel2[0], "N   ");
      SELECT(sel2[1], "CA  ");
      SELECT(sel2[2], "CB  ");
      break;
   case ATOMS_CNCA:
      fprintf(stderr,"INTERNAL ERROR: Code must be modified to support C,N,CA in second position\n");
      exit(1);
      SELECT(sel2[0], "C   ");
      SELECT(sel2[1], "N   ");
      SELECT(sel2[2], "CA  ");
      break;
   case ATOMS_CACO:
      SELECT(sel2[0], "CA  ");
      SELECT(sel2[1], "C   ");
      SELECT(sel2[2], "O   ");
      break;
   default:
      return(FALSE);
   }
   
   if((sel1[0] ==NULL)||(sel2[0]==NULL))
   {
      return(FALSE);
   }
   /* create linked list containing just C, N and CA atoms of residue 1 */
   if((keyres1_pdb = SelectAtomsResidue(res1_start, res1_stop,
                                        3, sel1, &natoms)) == NULL)
   {
      return(FALSE);
   }  

   /* If we are doing C,N,CA, we want to swap the carbon from res1 for the
      carbon from the previous residue
   */
   if(atomset1 == ATOMS_CNCA)
   {
      if((keyres1_pdb = SwapCarbon(keyres1_pdb, prevres1))==NULL)
      {
         return(FALSE);
      }
   }
   /* create linked list containing just C, N and CA atoms of residue 2 */
   if((partnerres2_pdb = SelectAtomsResidue(res2_start, res2_stop,
                                            3, sel2, &natoms)) == NULL)
   {
      return(FALSE);
   }

#ifdef DEBUG
   printf("original\n");
   WritePDB(stdout, keyres1_pdb);
   WritePDB(stdout, partnerres2_pdb);
#endif
   
   /*moving partnerres2_pdb  CA atom to the origin */
   tempv.x = -Vector.x;
   tempv.y = -Vector.y;
   tempv.z = -Vector.z;
   
   TranslatePDB(partnerres2_pdb, tempv);
    
#ifdef DEBUG
   printf("translate\n");
   WritePDB(stdout, keyres1_pdb);
   WritePDB(stdout, partnerres2_pdb);
#endif

   /* create coordinate arrays for res1 and res2 */
   NumCo_ord2 = GetPDBCoor(partnerres2_pdb, &partnerres2_coor);
   NumCo_ord1 = GetPDBCoor(keyres1_pdb, &keyres1_coor);

   /* create the weight array */
   if((weight = (REAL *)malloc(NumCo_ord1 * sizeof(REAL))) == NULL)
   {
      if(partnerres2_coor)
      {
         free(partnerres2_coor);
         partnerres2_coor = NULL;
      }
      if(keyres1_coor)
      {
         free(keyres1_coor);
         keyres1_coor = NULL;
      }
      return(FALSE);
   }

   /* Set up the weight array */
   resnum = partnerres2_pdb->resnum;
   insert = partnerres2_pdb->insert[0];
   start = partnerres2_pdb;
   natom = 0;
   i = 0;
   for(q = start; q!=NULL; NEXT(q))
   {
      if(!strncmp(q->atnam, sel2[0], 4))
         weight[i] = (REAL)1.0;
      if(!strncmp(q->atnam, sel2[1], 4))
         weight[i] = (REAL)1.0;
      if(!strncmp(q->atnam, sel2[2], 4))
         weight[i] = (REAL)(0.1);
      i++;
   }
   
   /* create rotation matrix. */
   if(!matfit(partnerres2_coor, keyres1_coor, gRotation_matrix,
              NumCo_ord2, weight, FALSE))
   {
      fprintf(stderr,"Fitting failed!\n");
      return(FALSE);
   }

#ifdef DEBUG   
   ApplyMatrixPDB(keyres1_pdb, gRotation_matrix);
   printf("rotate\n");
   WritePDB(stdout, keyres1_pdb);
   WritePDB(stdout, partnerres2_pdb);
#endif
   
   FREELIST(partnerres2_pdb, PDB);
   partnerres2_pdb = NULL;
   
   FREELIST(keyres1_pdb, PDB);
   keyres1_pdb = NULL;
   
   free(weight);
   weight = NULL;
   
   free(sel1[0]);
   sel1[0] = NULL;
   
   free(sel1[1]);
   sel1[1] = NULL;
   
   free(sel1[2]);
   sel1[2] = NULL;
   
   free(sel2[0]);
   sel2[0] = NULL;
   
   free(sel2[1]);
   sel2[1] = NULL;
   
   free(sel2[2]);
   sel2[2] = NULL;
   
   return(TRUE);
}

/************************************************************************/
/* function that assesses whether a hydrogen-bond between the *key* 
   residue and *partner* residue is valid. It compares the location of 
   the *key* residue hydrogen-capable atoms  with that the *partner* 
   hydrogen atoms (i.e. where the *partner* residue likes its partner 
   hydrogen atoms to be).
   Fits *partner* residue matrix on top of *key* residue matrix and 
   calculates the difference between their respective hydrogen-capable 
   atoms.
   A valid hydrogen bond is said to exist if the distance between the 
   two atoms are within a certain cutoff distance 
*/
BOOL CheckValidHBond(VEC3F CAtoCAVector, REAL cutoff, FILE *out,
                     int keyarray[MAXSIZE][MAXSIZE][MAXSIZE],
                     int partnerarray[MAXSIZE][MAXSIZE][MAXSIZE])
{
   int x, y, z, totalcount1 = 0, totalcount2 = 0, atnum = 0, resnum = 0;
   FILE *OUT = out;
   REAL pseudoenergy, final_penergy  = 9999.9999;
   VEC3F partner_coord;

#if defined(DEBUG1) || defined(DEBUG2) || defined(DEBUG3)
         gChain--;
#endif
   
#ifdef DEBUG3
{
   int x2, y2, z2;
   int x_coord, y_coord, z_coord;
   VEC3F rotated_coord;

   fprintf(stderr, "REMARK Key array: chain %c\n", gChain);

   for(x2=0; x2<MAXSIZE; x2++)
   {
      for(y2=0; y2<MAXSIZE; y2++)
      {
         for(z2=0; z2<MAXSIZE; z2++)
         {
            OrientateMatrix(CAtoCAVector, x2, y2, z2, &rotated_coord);
            COORD_2_GRID(x_coord,rotated_coord.x);
            COORD_2_GRID(y_coord,rotated_coord.y);
            COORD_2_GRID(z_coord,rotated_coord.z);
            if(VALIDGRIDCOORDS(x_coord, y_coord, z_coord))
            {
               if(keyarray[x_coord][y_coord][z_coord])
               {
                  fprintf(stdout, "ATOM  %5d  CA  ALA %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                          atnum++, gChain, resnum++, rotated_coord.x, rotated_coord.y,  rotated_coord.z, 1.00, 2.00);
               }
            }
         }
      }
   }
   fprintf(stdout, "TER   \n");
   gChain--;
}
#endif
         
   /* first ... compare *key* residue hydrogen-donor atoms
      with *partner* residue *partner-to-accept* donor atoms */ 
   totalcount1 = CalculateTotalCounts(keyarray);
   totalcount2 = CalculateTotalCounts(partnerarray);
      
#ifdef DEBUG3
   fprintf(stderr, "REMARK Partner array: chain %c\n", gChain);
#endif
   /* go though grid for *partner* residue */
   for(x=0; x < MAXSIZE; x++)
   {
      for(y=0; y < MAXSIZE; y++)
      {
         for(z = 0; z < MAXSIZE; z++)
         {
            pseudoenergy = DoCheckHBond(x, y, z, CAtoCAVector, &partner_coord, totalcount1,
                                        totalcount2, keyarray, partnerarray,
                                        cutoff, OUT);           
            
            if((pseudoenergy != -1) && (final_penergy > pseudoenergy))
            {
               final_penergy = pseudoenergy;
#ifdef DEBUG
   fprintf(stdout, "ATOM  %5d  C   THR  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
           atnum++, resnum++, partner_coord.x, partner_coord.y,  partner_coord.z, 1.00, 2.00);
#endif
            }
         }
      }
   }
   
   if(final_penergy != 9999.9999)
   {
      fprintf(out, "Pseudoenergy of best quality hydrogen bond: %.2f\n", 
              final_penergy);
      return(TRUE);
   }

   return(FALSE);
}

/************************************************************************/
REAL DoCheckHBond(int x, int y, int z, VEC3F CAtoCAVector, 
                  VEC3F *partner_coord, int totalcount1,
                  int totalcount2,
                  int keyarray[MAXSIZE][MAXSIZE][MAXSIZE],
                  int partnerarray[MAXSIZE][MAXSIZE][MAXSIZE],
                  REAL cutoff, FILE *out)
{
   VEC3F rotated_coord;
   int x_coord, y_coord, z_coord,
       count1 = 0, count2 = 0;
   int number_of_cells = 1+(cutoff / GRIDSPACING);
   int i, j, k;
   REAL cutoff_squared = cutoff * cutoff,
        pseudoenergy = -1, final_penergy = 9999.9999,
        dist_squared;
   int final_x,final_y,final_z;
   static int atnum = 0, resnum = 0;

   /* call routine to rotate matrix */
   if(partnerarray[x][y][z] > 0)
   {  
      count2 = partnerarray[x][y][z];
      OrientateMatrix(CAtoCAVector, x, y, z, &rotated_coord);
      
      /* convert real co-ordinates back into grid locations */
      COORD_2_GRID(x_coord,rotated_coord.x);
      COORD_2_GRID(y_coord,rotated_coord.y);
      COORD_2_GRID(z_coord,rotated_coord.z);

      if(VALIDGRIDCOORDS(x_coord, y_coord, z_coord))
      { 
#ifdef DEBUG2
   fprintf(stdout, "ATOM  %5d  C   THR %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
           atnum++, gChain, resnum++, rotated_coord.x, rotated_coord.y,  rotated_coord.z, 1.00, 2.00);
#endif
         count1 = keyarray[x_coord][y_coord][z_coord];
             
         /* if *key* and *partner* hydrogen atoms match exactly */
         if(count1 > 0)
         { 
            final_penergy = CalcEnergy(count1, totalcount1,
                                      count2, totalcount2);

            final_x = x_coord;
            final_y = y_coord;
            final_z = z_coord;

            GRID_2_COORD(final_x, partner_coord->x);
            GRID_2_COORD(final_y, partner_coord->y);
            GRID_2_COORD(final_z, partner_coord->z);

#ifdef DEBUG1
   fprintf(stdout, "ATOM  %5d  C   THR %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
           atnum++, gChain, resnum++, rotated_coord.x, rotated_coord.y,  rotated_coord.z, 1.00, 2.00);
#endif
         }
         else if(cutoff > 0.0)
         {
            for(i=x_coord-number_of_cells; i<=x_coord+number_of_cells; i++)
            {
               for(j=y_coord-number_of_cells; j<=y_coord+number_of_cells; j++)
               {
                  for(k=z_coord-number_of_cells; k<=z_coord+number_of_cells; k++)
                  {
                     if(VALIDGRIDCOORDS(i, j, k))
                     {
                        count1 = keyarray[i][j][k];
                     
                        if(count1 > 0)
                        {
                           dist_squared = Distance_squared(i, j, k, x_coord, y_coord, z_coord);
                           
                           if(dist_squared <=cutoff_squared)     
                           {
                              /* ACRM 08.09.02 Missing this VITAL line!!!! */
                              pseudoenergy = CalcEnergy(count1, totalcount1,
                                                        count2, totalcount2);

                              if(pseudoenergy < final_penergy)
                              {
                                 final_penergy = pseudoenergy;
                                 
                                 final_x = x_coord;
                                 final_y = y_coord;
                                 final_z = z_coord;
                                 
                                 GRID_2_COORD(final_x, partner_coord->x);
                                 GRID_2_COORD(final_y, partner_coord->y);
                                 GRID_2_COORD(final_z, partner_coord->z);
                              }
#ifdef DEBUG1
   fprintf(stdout, "ATOM  %5d  C   THR %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
           atnum++, gChain, resnum++, rotated_coord.x, rotated_coord.y,  rotated_coord.z, 1.00, 3.00);
#endif
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   return(final_penergy);
}
        
/************************************************************************/
REAL CalcEnergy(int count1, int totalcount1, int count2, int totalcount2)
{
   REAL potential1, potential2, log_potential;
   
   potential1 = ((REAL)count1/(REAL)totalcount1);
   potential2 = ((REAL)count2/(REAL)totalcount2);  
   log_potential = -log(potential1)-log(potential2);
   return(log_potential);
}

/************************************************************************/
void OrientateMatrix(VEC3F CAtoCAVector,int x, int y, int z, 
                     VEC3F *rotated_coord)
{
   VEC3F partner_real_coord, partner_real_coord_rotated;

   /* convert location of hydrogen atoms in *partner* residue
      matrix to real coordinates */
   GRID_2_COORD(x, partner_real_coord.x);
   GRID_2_COORD(y, partner_real_coord.y);
   GRID_2_COORD(z, partner_real_coord.z);

   /* multiply vector by rotation matrix */
   MatMult3_33(partner_real_coord, gRotation_matrix,
               &partner_real_coord_rotated);
  
   /* add CA to CA vector */
   partner_real_coord_rotated.x += CAtoCAVector.x;
   partner_real_coord_rotated.y += CAtoCAVector.y;
   partner_real_coord_rotated.z += CAtoCAVector.z;
  
   *rotated_coord = partner_real_coord_rotated;
}

/************************************************************************/
/* function that calculates distance (in angstroms) between two grid 
   points 
*/
REAL Distance_squared (int i, int j, int k, int x_coord,
                       int y_coord, int z_coord)
{
   REAL squared_dist;
   VEC3F test_value, c_value;
   
   GRID_2_COORD(i, test_value.x);
   GRID_2_COORD(j, test_value.y);
   GRID_2_COORD(k, test_value.z);

   GRID_2_COORD(x_coord, c_value.x);
   GRID_2_COORD(y_coord, c_value.y);
   GRID_2_COORD(z_coord, c_value.z);

   squared_dist = DISTSQ_HATOMS(test_value, c_value);
      
   return(squared_dist);   
}

/************************************************************************/
/* function that finds location of an atom in a residue */
PDB *FindAtom(PDB *res1_start, PDB *res1_next, char *atnam, VEC3F *atm)
{
   PDB *p,
      *found = NULL;
        
   for(p = res1_start; p !=res1_next; NEXT(p))
   {
      if(!strncmp(p->atnam, atnam, 4))
      {
         if(atm != NULL)
         {
            atm->x = p->x;
            atm->y = p->y;
            atm->z = p->z;
         }
         found = p;
         break;
      }
   }

   return(found);
}
   
/************************************************************************/
void CullArrays(PDB *pdb, PDB *res1, PDB *res2,
                int donate_array[MAXSIZE][MAXSIZE][MAXSIZE],
                int accept_array[MAXSIZE][MAXSIZE][MAXSIZE])
{
   PDB *p;
   int x_coord, y_coord, z_coord, x, y, z, total;
   
   for(p = pdb; p!=NULL; NEXT(p))
   {
      COORD_2_GRID(x_coord,p->x);
      COORD_2_GRID(y_coord,p->y);
      COORD_2_GRID(z_coord,p->z);    
      
      for(x = -RAD; x <=RAD; x++)
      {
         for(y = -RAD; y <=RAD; y++)
         {
            for(z = -RAD; z <=RAD; z++) 
            {
               total = ((x*x) + (y*y) + (z*z));
               if(total  < (RAD*RAD))
               {
                  x_coord = x + x_coord;
                  y_coord = y + y_coord;
                  z_coord = z + z_coord;

                  if(VALIDGRIDCOORDS(x_coord, y_coord, z_coord))
                  {
                     if((p->resnum !=res1->resnum)
                        && (p->chain != res1->chain)
                        && (p->insert != res1->insert)
                        && (p->resnum !=res2->resnum)
                        && (p->chain != res2->chain)
                        && (p->insert != res2->insert))
                     {
                        donate_array[x_coord][y_coord][z_coord] = 0; 
                        accept_array[x_coord][y_coord][z_coord] = 0; 
                     }
                  }
               }
            }
         }
      }
   }
}
            
/************************************************************************/
/* function that populates matrices from text file (created in 
   hydrogen_matrices.c program 
*/
BOOL ReadInMatrices(char *res1, char *res2, FILE *matrix, 
                    int type, int whichres)
{
   char buffer[MAXBUFF], junk[15], minibuffer[6];
   int x, y, z, count;
   BOOL inResidue1, inResidue2,
        found_residue1 = FALSE, 
        found_residue2 = FALSE;

   rewind(matrix);
   while(fgets(buffer, MAXBUFF, matrix))
   {     
      TERMINATE(buffer);
      if(!strncmp(buffer, "residue", 7))
      {
         strcpy(minibuffer, buffer+8);
        
         if((whichres&MAT_RES_1) && !strncmp(minibuffer, res1, 3))         
         {
            inResidue1 = TRUE;
         }
         else
         {
            inResidue1 = FALSE;
         }

         if((whichres&MAT_RES_2) && !strncmp(minibuffer, res2, 3))
         {
            inResidue2 = TRUE;       
         }
         else
         {
            inResidue2 = FALSE;
         }
      }

      if(inResidue1)
      {
         if(type&MAT_READ_DONOR1)
         { 
            if(!strncmp(buffer, "donate", 6))
            {
               sscanf(buffer, "%s %d %d %d %d", junk, &x, &y, &z, &count);
               gDonate[x][y][z] = count; 
            }
            if(!strncmp(buffer, "partnertodonate", 15))
            {
               sscanf(buffer, "%s %d %d %d %d", junk, &x, &y, &z, &count);          
               gPartnertoDonate[x][y][z] = count;
            }
         }
         if(type&MAT_READ_ACCEPTOR1)
         {
            if(!strncmp(buffer, "accept", 6))
            {
               sscanf(buffer, "%s %d %d %d %d", junk, &x, &y, &z, &count);
               gAccept[x][y][z] = count; 
            }
            if(!strncmp(buffer, "partnertoaccept", 15))
            {
               sscanf(buffer, "%s %d %d %d %d", junk, &x, &y, &z, &count);          
               gPartnertoAccept[x][y][z] = count;
            }
         }
         
         found_residue1 = TRUE;
      }
      
      if(inResidue2)
      {
         if(type&MAT_READ_DONOR2)
         { 
            if(!strncmp(buffer, "donate", 6))
            {
               sscanf(buffer, "%s %d %d %d %d", junk, &x, &y, &z, &count);
               gDonate[x][y][z] = count; 
            }
            if(!strncmp(buffer, "partnertodonate", 15))
            {
               sscanf(buffer, "%s %d %d %d %d", junk, &x, &y, &z, &count);          
               gPartnertoDonate[x][y][z] = count;
            }
         }
         if(type&MAT_READ_ACCEPTOR2)
         {
            if(!strncmp(buffer, "accept", 6))
            {
               sscanf(buffer, "%s %d %d %d %d", junk, &x, &y, &z, &count);
               gAccept[x][y][z] = count; 
            }
            if(!strncmp(buffer, "partnertoaccept", 15))
            {
               sscanf(buffer, "%s %d %d %d %d", junk, &x, &y, &z, &count);          
               gPartnertoAccept[x][y][z] = count;
            }
         }
 
         found_residue2 = TRUE; 
      }
   }
   
   if(whichres == MAT_RES_BOTH)
   {
      if(found_residue1 && found_residue2)
         return(TRUE);
   }
   else if(whichres == MAT_RES_1)
   {
      if(found_residue1)
         return(TRUE);
   }
   else if(whichres == MAT_RES_2)
   {
      if(found_residue2)
         return(TRUE);
   }
         
   return(FALSE);
}

/************************************************************************/
/* function that sets array elements to 0 */
void ClearArrays(int array)
{
   int i, j, k;
   
   for(i = 0; i < MAXSIZE; i++)
   {
      for(j = 0; j < MAXSIZE; j++)
      {
         for(k = 0; k < MAXSIZE; k++)
         {
            gDonate[i][j][k] = 0;
            gAccept[i][j][k] = 0;
            gPartnertoDonate[i][j][k] = 0;
            gPartnertoAccept[i][j][k] = 0;
         }
      }
   }  
}

/************************************************************************/
/* function to open matrix file. Opens default matrix file if none 
   specified on command line 
*/

FILE *OpenMatrixFile(char *matrix_file, char *def_matrix_file)
{
   FILE *fp;

   /* if filename has been specified on command line, open and return */
   if(matrix_file !=NULL && matrix_file[0])
   {      
      fp = fopen(matrix_file, "r");
   }
   /* if no filename specified, open default matrix file */
   else
   {
      fp = fopen(def_matrix_file, "r");
   }

   return(fp);
}

/************************************************************************/
/* function to parse the command line */
BOOL ParseCmdLine(int argc, char **argv, REAL *cutoff, BOOL *hbplus,
                  char *hatom1, char *hatom2,
                  char *matrix_file, char *matrix_file2, char *pdbfile,
                  char *locres1, char *locres2, char *res2, 
                  char *outputfile)
{
   argc--;
   argv++;

   *cutoff = DEFAULT_CUTOFF_VALUE;
   *hbplus = FALSE;

   matrix_file[0] = '\0';
   pdbfile[0] = outputfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'c':
            argc--;
            argv++;
            if((!argc) || !sscanf(argv[0], "%lf", cutoff))
               return(FALSE);
            break;
         case 'p':
            *hbplus = TRUE;
            argc--;
            argv++;
            if((!argc) || !sscanf(argv[0], "%s", hatom1))
               return(FALSE);
            argc--;
            argv++;
            if((!argc) || !sscanf(argv[0], "%s", hatom2))
               return(FALSE);
            break;
         case 'm':
            argc--;
            argv++;
            strcpy(matrix_file, argv[0]);
            break;
         case 'n':
            argc--;
            argv++;
            strcpy(matrix_file2, argv[0]);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* check there are 5 or 6 arguments remaining */
         if((argc > 5) || (argc < 4))
            return(FALSE);
         
         strcpy(pdbfile, argv[0]);
         argc--;
         argv++;
         strcpy(locres1, argv[0]);
         UPPER(locres1);
         argc--;
         argv++;
         strcpy(locres2, argv[0]);
         UPPER(locres2);
         argc--;
         argv++;
         strcpy(res2, argv[0]);
         UPPER(res2);
         argc--;
         argv++;
              
         if(argc)
         {
            strcpy(outputfile, argv[0]);
         }
         return(TRUE);
      }
      argc--;
      argv++;
   }
   return(FALSE);
}

/************************************************************************/
/* function to display a usage message */
void Usage(void)
{
   fprintf(stderr, "\nCheckHBond V2.0 (c) 2002-6, Alison Cuff, University of Reading\n\n");
   fprintf(stderr, "V2.0 changes by Andrew Martin, UCL\n\n")
   fprintf(stderr, "Usage: checkhbond [-c cutoff] [-p hatom1 hatom2][-m matrix file]\n");
   fprintf(stderr, "   pdbfile residue1 residue2 nameres2 [output file]\n\n");
   fprintf(stderr, "  -c [cutoff]: cutoff distance between hydrogen-capable atoms(default: 0.5A)\n");
   fprintf(stderr, "  -p: Parse HBplus data.\n");
   fprintf(stderr, "  Hydrogen donating atom (hatom1) and hydrogen accepting atom (hatom2) required \n");
   fprintf(stderr, "  -m [matrix file]: matrix file (if not using default file\n");
   fprintf(stderr, "  pdbfile:  pdb file of protein structure\n");
   fprintf(stderr, "  residue1: First residue (chain, residue number, insert)\n");
   fprintf(stderr, "  residue2: Second residue (chain, residue number, insert)\n");
   fprintf(stderr, "  nameres2:  Name of amino acid to test at residue 2\n");
   fprintf(stderr, "      (needs native residue 1 if creating 'pseudo-energy distribution'\n");
   fprintf(stderr, "  [output file]: for saving hydrogen-capable atoms (in pdb format)\n");
   fprintf(stderr, "      I/O is through stdout if file not specified\n\n");   
   fprintf(stderr, "Determines the validity of a hydrogen bond in a protein \n");
   fprintf(stderr, "Assesses if hydrogen bond is maintained if one amino acid is substituted\n");
   fprintf(stderr, "for another and calculates a pseudoenergy for any that are maintained\n");
   fprintf(stderr, "With -p switch, creates 'reference' pseudoenergy values for real hydrogen bonds\n");
   fprintf(stderr, "against which the pseudoenergy of hydrogen bonds maintained after mutation can be\n");
   fprintf(stderr, "assessed for their quality\n\n");
}

/************************************************************************/
/* Function that opens input file for reading and output file for writing
   to and appending 
*/
BOOL Open_Std_Files(char *infile, char *outfile, FILE **in, FILE **out)
{
   if(infile!=NULL && infile[0] && strcmp(infile,"-"))
   {
      if((*in = fopen(infile,"r"))==NULL)
      {
         fprintf(stderr,"Unable to open input file: %s\n",infile);
         return(FALSE);
      }
   }
      
   if(outfile!=NULL && outfile[0] && strcmp(outfile,"-"))
   {
      if((*out = fopen(outfile,"a+"))==NULL)
      {
         fprintf(stderr,"Unable to open output file: %s\n",outfile);
         return(FALSE);
      }
   }
   
   return(TRUE);
}

/************************************************************************/
/* Function that recognises any residues not capable of hydrogen bonding
 */
BOOL IsHBondCapable(char *residue)
{   
   int i;
   char *NotHBondRes[9];
   BOOL flag = FALSE;
   
   NotHBondRes[0] = "MET";
   NotHBondRes[1] = "CYS";
   NotHBondRes[2] = "ILE";
   NotHBondRes[3] = "VAL";
   NotHBondRes[4] = "PHE";
   NotHBondRes[5] = "GLY";
   NotHBondRes[6] = "ALA";
   NotHBondRes[7] = "PRO";
   NotHBondRes[8] = "LEU";

   for(i = 0; i < 9; i++)
   {
      if(!strcmp(residue, NotHBondRes[i]))
      {
         flag = TRUE;
         break;
      } 
      else
      {
         flag = FALSE;
      }
   }
   
   return(flag);
}

/************************************************************************/
PDB *GetResidues(PDB *pdb, char *chain1, int resnum1, char *insert1, 
                 char *chain2, int resnum2, char *insert2, int *errorcode)
{
   PDB *keep=NULL,
       *p, *q,
       *prev1, *prev2;
   
   /* First find the residue preceeding res1 */
   for(p=pdb, prev1=NULL; p!=NULL; )
   {
      if((p->resnum == resnum1) &&
         (p->chain[0]  == chain1[0]) &&
         (p->insert[0] == insert1[0]))
         break;
      prev1 = p;
      p = FindNextResidue(p);
   }
   if((prev1==NULL) || !ResiduesBonded(prev1, p))
   {
      *errorcode = ERR_NOPREVRES1;
      return(NULL);
   }

   /* First find the residue preceeding res2 */
   for(p=pdb, prev2=NULL; p!=NULL; )
   {
      if((p->resnum == resnum2) &&
         (p->chain[0]  == chain2[0]) &&
         (p->insert[0] == insert2[0]))
         break;
      prev2 = p;
      p = FindNextResidue(p);
   }
   if((prev2==NULL) || !ResiduesBonded(prev2, p))
   {
      *errorcode = ERR_NOPREVRES2;
      return(NULL);
   }


   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(((p->resnum == resnum1) &&
          (p->chain[0]  == chain1[0]) &&
          (p->insert[0] == insert1[0])) ||
         ((p->resnum == resnum2) &&
          (p->chain[0]  == chain2[0]) &&
          (p->insert[0] == insert2[0])) ||
         ((p->resnum == prev1->resnum) &&
          (p->chain[0]  == prev1->chain[0]) &&
          (p->insert[0] == prev1->insert[0])) ||
         ((p->resnum == prev2->resnum) &&
          (p->chain[0]  == prev2->chain[0]) &&
          (p->insert[0] == prev2->insert[0])))
      {
         if(keep == NULL)
         {
            INIT(keep, PDB);
            q = keep;
         }
         else
         {
            ALLOCNEXT(q, PDB);
         }
         if(q==NULL)
         {
            *errorcode = ERR_NOMEM;
            return(NULL);
         }
         CopyPDB(q, p);
      }
   }
   FREELIST(pdb, PDB);

   return(keep);
}


/************************************************************************/
void FindRes1Type(PDB *pdb, char *chain, int resnum, char *insert, 
                  char *res)
{
   PDB *p;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->resnum == resnum) &&
         (p->chain[0] == chain[0]) &&
         (p->insert[0] == insert[0]))
      {
         strncpy(res, p->resnam, 3);
         res[3] = '\0';
         return;
      }
   }
}

/************************************************************************/
/* Tests if two residues are bonded through C...N
   The residues must be specified in the correct order!
*/
BOOL ResiduesBonded(PDB *pdb1, PDB *pdb2)
{
   PDB *c, *n, *next1, *next2;
   
   next1 = FindNextResidue(pdb1);
   next2 = FindNextResidue(pdb2);
   for(c=pdb1; c!=next1; NEXT(c))
   {
      if(!strncmp(c->atnam, "C   ", 4))
         break;
   }
   for(n=pdb2; n!=next2; NEXT(n))
   {
      if(!strncmp(n->atnam, "N   ", 4))
         break;
   }
   if((c==next1) || (n==next2))
      return(FALSE);

   if(DISTSQ(c,n) < CNBONDSQ)
      return(TRUE);

   return(FALSE);
}

/************************************************************************/
/* function that calculates the vector from the N of res1 to CA of res2 
 */
void CalculateNToCaVector(PDB *res1_start, PDB *res1_stop,
                          PDB *res2_start, PDB *res2_stop, 
                          VEC3F *NtoCAVector)
{
   VEC3F res1_n, 
      res2_calpha;

   /*find co-ordinates of N  atom */   
   if(!FindAtom(res1_start, res1_stop, "N   ", &res1_n))
   {
      fprintf(stderr, "Error: Can't find N atom of key residue\n");
   }
   
   /*find co-ordinates of CA atom */   
   if(!FindAtom(res2_start, res2_stop, "CA  ", &res2_calpha))
   {
      fprintf(stderr, "Error: Can't find c-alpha atom of partner residue\n");
   }
   
   NtoCAVector->x =  res2_calpha.x - res1_n.x;
   NtoCAVector->y =  res2_calpha.y - res1_n.y;
   NtoCAVector->z =  res2_calpha.z - res1_n.z;
}


/************************************************************************/
/* function that calculates the vector from the C of res1 to CA of res2 
 */
void CalculateCToCaVector(PDB *res1_start, PDB *res1_stop,
                          PDB *res2_start, PDB *res2_stop, 
                          VEC3F *CtoCAVector)
{
   VEC3F res1_n, 
      res2_calpha;

   /*find co-ordinates of N atom */   
   if(!FindAtom(res1_start, res1_stop, "C   ", &res1_n))
   {
      fprintf(stderr, "Error: Can't find C atoms of key residue\n");
   }
   
   /*find co-ordinates of CA atom */   
   if(!FindAtom(res2_start, res2_stop, "CA  ", &res2_calpha))
   {
      fprintf(stderr, "Error: Can't find c-alpha atoms of partner residue\n");
   }
   
   CtoCAVector->x =  res2_calpha.x - res1_n.x;
   CtoCAVector->y =  res2_calpha.y - res1_n.y;
   CtoCAVector->z =  res2_calpha.z - res1_n.z;
}


/************************************************************************/
BOOL AnalyzeMCDonorPair(int resnum1, int resnum2, PDB *pdb, 
                        FILE *matrix, FILE *matrix2,
                        char *chain1, char *chain2, 
                        char *insert1, char *insert2, BOOL *hbplus,
                        char *hatom1, char *hatom2, REAL cutoff, 
                        char *res1, char *res2,
                        FILE *OUT)
{
   VEC3F NtoCAVector;

   PDB *res1_start, *res1_stop, *res2_start, *res2_stop, *prevres1;
   
   /* Find the key residue which is the mainchain donor and the previous residue */
   for(res1_start =  pdb, prevres1 = NULL;
       res1_start != NULL;
       res1_start =  res1_stop)
   {
      res1_stop = FindNextResidue(res1_start);
                        
      /* if *key* residue of choice */
      if((resnum1 == res1_start->resnum)
         && (chain1[0]==res1_start->chain[0])
         && (insert1[0] == res1_start->insert[0]))
      {
         break;
      }
      prevres1 = res1_start;
   }
   if(res1_start==NULL)
   {
      fprintf(stderr, "Error: Can't find residue 1  in PDB file\n");
      return(FALSE);
   }
   
   /*Find the partner residue which is the sidechain acceptor */
   for(res2_start = pdb; res2_start !=NULL;
       res2_start=res2_stop)
   {
      res2_stop = FindNextResidue(res2_start);
      
      /* if *partner* residue of choice */
      if((resnum2 == res2_start->resnum)
         && (chain1[0]==res1_start->chain[0])
         && (insert1[0] == res1_start->insert[0]))
      {
         break;
      }
   }
   if(res2_start==NULL)
   {
      fprintf(stderr, "Error: Can't find residue 2 in PDB file\n");
      return(FALSE);
   }
   
   if(!ReadInMatrices(res1, res2, matrix, MAT_READ_ACCEPTOR2, MAT_RES_2))
   {
      /* ACRM 08.09.05 - error message */
      fprintf(stderr, "Unable to find residue in the sc/sc matrix file: %s\n", res2);
      return(FALSE);
   }
   
   if(!ReadInMatrices(res1, res2, matrix2, MAT_READ_DONOR1, MAT_RES_1))
   {
      /* ACRM 08.09.05 - error message */
      fprintf(stderr, "Unable to find residue in the mc donor matrix file: %s\n", res1);
      return(FALSE);
   }
   
   OrientatePDB(pdb, res2_start, res2_stop);
#ifndef NOCULL
   CullArrays(pdb, res1_start, res2_start, gPartnertoAccept,
              gAccept);
#endif
   OrientateN_PDB(pdb, prevres1, res1_start, res1_stop);
#ifndef NOCULL
   CullArrays(pdb, res1_start, res2_start,
              gDonate, gPartnertoDonate);
#endif   
   CalculateNToCaVector(res1_start, res1_stop, res2_start,
                         res2_stop, &NtoCAVector);
   
   CreateRotationMatrix(pdb, res1_start, res1_stop, res2_start,
                        res2_stop, NtoCAVector, ATOMS_CNCA, ATOMS_NCACB, prevres1);

   if(*hbplus == TRUE)    /* ACRM 23.07.04 Dereference pointer! */
   {
      OrientatePDB(pdb, res2_start, res2_stop);
      if(!CalculateHBondEnergy(pdb, res1_start, res1_stop,
                               res2_start,
                               res2_stop, NtoCAVector, cutoff,
                               hatom1, hatom2, OUT))
      {
         fprintf(stderr, "No hydrogen bonds\n");
      }
   }
   else
   {
      if(!CheckValidHBond(NtoCAVector, cutoff, OUT,
                          gDonate, gPartnertoAccept))
      {
         if(!CheckValidHBond(NtoCAVector, cutoff, OUT,
                             gAccept, gPartnertoDonate))
         {
            return(FALSE);
         }
      }
   }
   return(TRUE);
}

/************************************************************************/
BOOL AnalyzeMCAcceptorPair(int resnum1, int resnum2, PDB *pdb, 
                           FILE *matrix, FILE *matrix2,
                           char *chain1, char *chain2, 
                           char *insert1, char *insert2, BOOL *hbplus,
                           char *hatom1, char *hatom2, REAL cutoff, 
                           char *res1, char *res2,
                           FILE *OUT)
{
   VEC3F CtoCAVector;

   PDB *res1_start, *res1_stop, *res2_start, *res2_stop;
   
   /* Find the key residue which is the mainchain acceptor */
   for(res1_start =  pdb;
       res1_start != NULL;
       res1_start =  res1_stop)
   {
      res1_stop = FindNextResidue(res1_start);
                        
      /* if *key* residue of choice */
      if((resnum1 == res1_start->resnum)
         && (chain1[0]==res1_start->chain[0])
         && (insert1[0] == res1_start->insert[0]))
      {
         break;
      }
   }
   if(res1_start==NULL)
   {
      fprintf(stderr, "Error: Can't find residue 1  in PDB file\n");
      return(FALSE);
   }
   
   /*Find the partner residue which is the sidechain acceptor */
   for(res2_start = pdb; res2_start !=NULL;
       res2_start=res2_stop)
   {
      res2_stop = FindNextResidue(res2_start);
      
      /* if *partner* residue of choice */
      if((resnum2 == res2_start->resnum)
         && (chain1[0]==res1_start->chain[0])
         && (insert1[0] == res1_start->insert[0]))
      {
         break;
      }
   }
   if(res2_start==NULL)
   {
      fprintf(stderr, "Error: Can't find residue 2 in PDB file\n");
      return(FALSE);
   }
   
   if(!ReadInMatrices(res1, res2, matrix, MAT_READ_DONOR2, MAT_RES_2))
   {
      /* ACRM 08.09.05 - error message */
      fprintf(stderr, "Unable to find residue in the sc/sc matrix file: %s\n", res2);
      return(FALSE);
   }
   
   if(!ReadInMatrices(res1, res2, matrix2, MAT_READ_ACCEPTOR1, MAT_RES_1))
   {
      /* ACRM 08.09.05 - error message */
      fprintf(stderr, "Unable to find residue in the mc acceptor matrix file: %s\n", res1);
      return(FALSE);
   }
   
   OrientatePDB(pdb, res2_start, res2_stop);
#ifndef NOCULL
   CullArrays(pdb, res1_start, res2_start, gDonate,
              gPartnertoDonate);
#endif
   OrientateCO_PDB(pdb, res1_start, res1_stop);
#ifndef NOCULL
   CullArrays(pdb, res1_start, res2_start,
              gAccept, gPartnertoAccept);
#endif   
   CalculateCToCaVector(res1_start, res1_stop, res2_start,
                         res2_stop, &CtoCAVector);
   
   CreateRotationMatrix(pdb, res1_start, res1_stop, res2_start,
                        res2_stop, CtoCAVector, ATOMS_CACO, ATOMS_NCACB, NULL);

   if(*hbplus == TRUE)    /* ACRM 23.07.04 Dereference pointer! */
   {
      OrientatePDB(pdb, res2_start, res2_stop);
      if(!CalculateHBondEnergy(pdb, res1_start, res1_stop,
                               res2_start,
                               res2_stop, CtoCAVector, cutoff,
                               hatom1, hatom2, OUT))
      {
         fprintf(stderr, "No hydrogen bonds\n");
      }
   }
   else
   {
      if(!CheckValidHBond(CtoCAVector, cutoff, OUT,
                          gDonate, gPartnertoAccept))
      {
         if(!CheckValidHBond(CtoCAVector, cutoff, OUT,
                             gAccept, gPartnertoDonate))
         {
            return(FALSE);
         }
      }
   }
   return(TRUE);
}

/************************************************************************/
/* Takes a PDB linked list of a single residue (pdb) and removes the
   C atom. Then takes a copy of the C atom from another residue (prev)
   and prepends it onto the first linked list
*/
PDB *SwapCarbon(PDB *pdb, PDB *prevres)
{
   PDB *p, *c, *newc, *stop, *prev;
   
   /* First we drop the existing carbon */
   for(p=pdb, prev=NULL; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam, "C   ", 4))
      {
         if(prev==NULL)
         {
            pdb = p->next;
            free(p);
            break;
         }
         else
         {
            prev->next = p->next;
            free(p);
            break;
         }
      }
      prev = p;
   }

   stop = FindNextResidue(prevres);
   if((c = FindAtom(prevres, stop, "C   ", NULL))==NULL)
   {
      return(NULL);
   }
   INIT(newc, PDB);
   if(newc==NULL)
   {
      return(NULL);
   }
   CopyPDB(newc,c);
   newc->next = pdb;
   
   return(newc);
}
