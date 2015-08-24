/********************************************************

Description:

Checkhbond program

If a hydrogen-capable residue is
substituted for another will the hydrogen-bond between the replacement residue 
and the partner residue of the original be maintained

Steps:

Orientate protein so that res 1 has its CA at the origin.
Cull matrices to ensure that there are no steric clashes

vector is calculated from CA of res 1 to CA of res 2
res 2 is moved with CA at the origin
matfit is then used to obtain a rotation matrix (with res 2 as the 'mobile' residue)

Run though grid for res 2 for all non-zero positions
for each position:
    convert to real x, y, z co-ordinates
    multiply the resulting vector by the rotation matrix -- res 2 grid points should be 
    at same orientation now as res 1 grid points
    move the vector by the CatoCa vector so that both grids are with CA at the origin
    Convert real co-ordinates back to grid locations and match with locations in other grid
*/

/*******************************************************/
 
/*includes */

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

/*************************************/

/* defines and macros */
/* default matrix file */
/*#define MATRIXFILE "/home/alison/hydrogen_bonding/matrices_05new.txt"*/
/*#define MATRIXFILE "/home/alison/hydrogen_bonding/hydrogen_bonding_new/testmatrices.txt"*/
#define MATRIXFILE "/home/alison/hydrogen_bonding/hydrogen_bonding_new/matrices2.txt"
#define MAXBUFF 300
/* radius of atom */
#define RAD 25
/* default cutoff dostance between atoms */
#define DEFAULT_CUTOFF_VALUE 0.5
/* calculates distance in angstroms between atoms */
#define DISTSQ_HATOMS(a,b) (a.x - b.x) * (a.x - b.x) + \
                       (a.y - b.y) * (a.y - b.y) + \
                       (a.z - b.z) * (a.z - b.z)

/* for debugging purposes */
/*#define DEBUG */

/*************************************/

/* global variables */

/* matrices that store all acceptor and donor atoms */
int gDonate[MAXSIZE][MAXSIZE][MAXSIZE];
int gAccept[MAXSIZE][MAXSIZE][MAXSIZE];
/* matrix storing partner atoms to hydrogen accepting atoms */
int gPartnertoAccept[MAXSIZE][MAXSIZE][MAXSIZE];
int gPartnertoDonate[MAXSIZE][MAXSIZE][MAXSIZE];
/* rotation matrix */
REAL gRotation_matrix[3][3];

/*************************************/

/* prototypes */

int main (int argc, char *argv[]);
void Usage(void);
BOOL ReadInMatrices(char *res1, char *res2, FILE *matrix);
void CullArrays(PDB *pdb, PDB *res1, PDB *res2,
                int donate_array[MAXSIZE][MAXSIZE][MAXSIZE],
                int accept_array[MAXSIZE][MAXSIZE][MAXSIZE]);
void CalculateCaToCaVector(PDB *res1_start, PDB *res1_stop,
                           PDB *res2_start, PDB *res2_stop,
                           VEC3F *CAtoCAVector);
BOOL FindCAAtom(PDB *start, PDB *stop, VEC3F *c_alpha);
BOOL CreateRotationMatrix(PDB *pdb, PDB *res1_start,
                          PDB *res1_stop, PDB *res2_start,
                          PDB *res2_stop, VEC3F CAtoCAVector);
BOOL CheckValidHBond(VEC3F CAtoCAVector, REAL cutoff,
                     FILE *out,
                     int keyarray[MAXSIZE][MAXSIZE][MAXSIZE],
                     int partnerarray[MAXSIZE][MAXSIZE][MAXSIZE]);
void ClearArrays();
BOOL ParseCmdLine(int argc, char **argv, REAL *cutoff,
                  BOOL *hbplus, char *hatom1,
                  char *hatom2, char *matrix_file, char *pdbfile,
                  char *locres1, char *locres2, char *res1,char *res2, char *outputfile);
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
FILE *OpenMatrixFile(char *matrix_file);
BOOL Open_Std_Files(char *infile, char *outfile, FILE **in, FILE **out);
BOOL IsHBondCapable(char *residue);

/*************************************/

int main (int argc, char *argv[])
{
   FILE *PDBFILE = stdin,
      *OUT = stdout,
      *matrix = stdout;
   char locres1[6], locres2[6],res1[4], res2[6],
      pdbfile[MAXBUFF], outputfile[MAXBUFF];
   char chain1[6], insert1[6], chain2[6], insert2[6], hatom1[6], hatom2[6];
   char matrix_file[MAXBUFF];
   PDB *pdb, *res1_start, *res1_stop, *res2_start, *res2_stop;
   int natoms, resnum1, resnum2;
   REAL cutoff;
   VEC3F CAtoCAVector;
   BOOL hbplus = FALSE;
   
   /* set array elements to 0 */
   ClearArrays();
   
   if(ParseCmdLine(argc, argv,  &cutoff, &hbplus,  hatom1, hatom2,
                   matrix_file, pdbfile, locres1, locres2, res1, res2, outputfile))
   {
      if((matrix = OpenMatrixFile(matrix_file)))
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
                     /*step though each residue in list
                       is it *key* residue? */
                     for(res1_start = pdb; res1_start !=NULL;
                         res1_start=res1_stop)
                     {
                        
                        res1_stop = FindNextResidue(res1_start);
                        
                        /* if *key* residue of choice */
                        if((resnum1 == res1_start->resnum)
                           && (chain1[0]==res1_start->chain[0])
                           && (insert1[0] == res1_start->insert[0]))
                        {
                           strcpy(res1, res1_start->resnam);
                           break;
                        }
                        else
                        {
                           if(res1_start==NULL)
                           {
                              fprintf(stderr, "Error: Can't find residue 1  in PDB file\n");
                           }
                           
                        }
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
                        else
                        {
                           if(res2_start==NULL)
                           {
                              fprintf(stderr, "Error: Can't find residue 2 in PDB file\n");
                           }
                        }
                     }
                     
                     if(IsHBondCapable(res1) != 0)
                     {
                        fprintf(OUT, "No hydrogen bonds\n");
                        return(1);
                     }
                     
                     if(IsHBondCapable(res2) !=0)
                     {
                        fprintf(OUT, "No hydrogen bonds\n");
                        return(1);
                     }
                     
                     if(!ReadInMatrices(res1, res2, matrix))
                     {
                        fprintf(OUT, "No hydrogen bonds\n");
                        return(1);
                     }
                     
                     OrientatePDB(pdb, res2_start, res2_stop);
                     CullArrays(pdb, res1_start, res2_start, gPartnertoDonate,
                                gPartnertoAccept);
                     
                     OrientatePDB(pdb, res1_start, res1_stop);
                     CullArrays(pdb, res1_start, res2_start,
                                gDonate, gAccept);
                     
                     CalculateCaToCaVector(res1_start, res1_stop, res2_start,
                                           res2_stop, &CAtoCAVector);
                     CreateRotationMatrix(pdb, res1_start, res1_stop, res2_start,
                                          res2_stop, CAtoCAVector );
                     
                     if(hbplus)
                     {
                        
                        OrientatePDB(pdb, res2_start, res2_stop);
                        if(!CalculateHBondEnergy(pdb, res1_start, res1_stop,
                                                 res2_start,
                                                 res2_stop, CAtoCAVector, cutoff,
                                                 hatom1, hatom2, OUT))
                        {
                           fprintf(OUT, "No hydrogen bonds\n");
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
                              fprintf(OUT, "No hydrogen bonds\n");
                           }
                           
                        }
                        
                     }
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
/*************************************/

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
      fprintf(out, "Pseudoenergy of best quality hydrogen bond: %f\n",
              pseudoenergy);
      printf("Pseudoenergy of best quality hydrogen bond: %f\n",
             pseudoenergy);
      return(TRUE);
      
   }
   else
   {
      fprintf(out, "No hydrogen bonds\n");
      return(FALSE);
   }
   
      
}
/*************************************/
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

/*************************************/
/* function that calculates the vector from the CA of res1
 to CA of res2 */
void CalculateCaToCaVector(PDB *res1_start, PDB *res1_stop,
 PDB *res2_start, PDB *res2_stop, VEC3F *CAtoCAVector)
{
   VEC3F res1_calpha, 
      res2_calpha;

   /*find co-ordinates of CA atom */   
   if(!FindCAAtom(res1_start, res1_stop, &res1_calpha))
   {
      fprintf(stderr, "Error: Can't find c-alpha atoms of key residue\n");
   }
   
   if(!FindCAAtom(res2_start, res2_stop, &res2_calpha))
   {
      fprintf(stderr, "Error: Can't find c-alpha atoms of partner residue\n");
   }
   
   CAtoCAVector->x =  res2_calpha.x - res1_calpha.x;
   CAtoCAVector->y =  res2_calpha.y - res1_calpha.y;
   CAtoCAVector->z =  res2_calpha.z - res1_calpha.z;

       
}

/*************************************/
/* function that creates a rotation matrix (global) to fit res 1
   (*key* residue) onto res 2 (*partner* residue). A weight of 1.0
   is added to atoms CA and N and 0.1 to atom C. 
   Includes option to print out residues at set stages for debugging
   purposes.
*/
BOOL CreateRotationMatrix(PDB *pdb, PDB *res1_start, PDB *res1_stop,
 PDB *res2_start, PDB *res2_stop,VEC3F CAtoCAVector)
{
   
   VEC3F tempv;   
   
   PDB *keyres1_pdb = NULL, 
      *partnerres2_pdb = NULL;
   

   

   char *sel[4];


   
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
   
 
   
  	
    
   SELECT(sel[0], "N   ");
   SELECT(sel[1], "CA  ");
   SELECT(sel[2], "C   ");
 
 
   if(sel[0] ==NULL)
   {
      return(FALSE);
   }
   /* create linked list containing just C, N and CA atoms of residue 1 */
   if((keyres1_pdb = SelectAtomsResidue(res1_start, res1_stop,
                                        3, sel, &natoms)) == NULL)
   {
      return(FALSE);
      
   }   

   /* create linked list containing just C, N and CA atoms of residue 2 */
   if((partnerres2_pdb = SelectAtomsResidue(res2_start, res2_stop,
                                            3, sel, &natoms)) == NULL)
   {
      return(FALSE);
      
   }


#ifdef DEBUG
   printf("original\n");
   WritePDB(stdout, keyres1_pdb);
   WritePDB(stdout, partnerres2_pdb);
#endif
   
   /*moving partnerres2_pdb  CA atom to the origin */
   tempv.x = -CAtoCAVector.x;
   tempv.y = -CAtoCAVector.y;
   tempv.z = -CAtoCAVector.z;
   
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
      if(partnerres2_coor) free(partnerres2_coor);
      if(keyres1_coor) free(keyres1_coor);
      return(FALSE);
   }
   
   resnum = partnerres2_pdb->resnum;
   insert = partnerres2_pdb->insert[0];
   start = partnerres2_pdb;
   natom = 0;
   i = 0;

   for(p=partnerres2_pdb; p!=NULL; NEXT(p))
   {
      if(p->resnum !=resnum || p->insert[0] !=insert)
      {
         for(q=start; q!=p; NEXT(q))
         {
            if(!strncmp(q->atnam, "C   ", 4))
               weight[i] = (REAL)(0.1);
            if(!strncmp(q->atnam, "CA  ", 4))
               weight[i] = (REAL)(1.0);
            if(!strncmp(q->atnam, "N   ", 4))
               weight[i] = (REAL)(1.0);
            i++;
            
         }
         natom = 0;
         start = p;
         
      }
      natom++;
   }
   
   for(q = start; q!=NULL; NEXT(q))
   {
      if(!strncmp(q->atnam, "C   ", 4))
         weight[i] = (REAL)(0.1);
      if(!strncmp(q->atnam, "CA  ", 4))
         weight[i] = (REAL)1.0;
      if(!strncmp(q->atnam, "N   ", 4))
         weight[i] = (REAL)1.0;
      i++;
   }
   
   /* create rotation matrix. */
   if(!matfit(partnerres2_coor, keyres1_coor, gRotation_matrix,
              NumCo_ord2, weight, FALSE))
   {       
      return(FALSE);
   }

#ifdef DEBUG   
   ApplyMatrixPDB(keyres1_pdb, gRotation_matrix);
   
   printf("rotate\n");
   
   WritePDB(stdout, keyres1_pdb);
   WritePDB(stdout, partnerres2_pdb);
#endif
   
   FREELIST(partnerres2_pdb, PDB);
   FREELIST(keyres1_pdb, PDB);
   free(weight);
   free(sel[0]);
   free(sel[1]);
   free(sel[2]);
   
   return(TRUE);
   
}

/*************************************/
/* function that assesses whether a hydrogen-bond between the *key* residue
   and *partner* residue is valid. It compares the location of the *key*
   residue hydrogen-capable atoms  with that the *partner* hydrogen atoms
   (i.e. where the *partner* residue likes its partner hydrogen atoms to be).
   Fits *partner* residue matrix on top of *key* residue matrix and calculates
   the difference between their respective hydrogen-capable atoms.
   A valid hydrogen bond is said to exist if the distance between the two atoms
   are within a certain cutoff distance */
BOOL CheckValidHBond(VEC3F CAtoCAVector, REAL cutoff, FILE *out,
                     int keyarray[MAXSIZE][MAXSIZE][MAXSIZE],
                     int partnerarray[MAXSIZE][MAXSIZE][MAXSIZE])
{
   
   int x, y, z, totalcount1 = 0, totalcount2 = 0, atnum = 0, resnum = 0;

   FILE *OUT = out;
      
   REAL pseudoenergy, final_penergy  = 9999.9999;

   VEC3F partner_coord;
   

   /* first ... compare *key* residue hydrogen-donor atoms
      with *partner* residue *partner-to-accept* donor atoms */ 
   

   totalcount1 = CalculateTotalCounts(keyarray);
   totalcount2 = CalculateTotalCounts(partnerarray);
      
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
      fprintf(out, "Pseudoenergy of best quality hydrogen bond: %f\n", 
              final_penergy);
      return(TRUE);
      
   }
   else
      return(FALSE);
      
}

/****************************************/

REAL DoCheckHBond(int x, int y, int z, VEC3F CAtoCAVector, VEC3F *partner_coord, int totalcount1,
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

   int atnum = 0, resnum = 0;
         
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
   fprintf(stdout, "ATOM  %5d  C   THR  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
           atnum++, resnum++, rotated_coord.x, rotated_coord.y,  rotated_coord.z, 1.00, 2.00);
#endif
            
 
         }
         
         else
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
   fprintf(stdout, "ATOM  %5d  C   THR  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
           atnum++, resnum++, rotated_coord.x, rotated_coord.y,  rotated_coord.z, 1.00, 3.00);
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
        
/***************************************/
REAL CalcEnergy(int count1, int totalcount1, int count2, int totalcount2)
{
   REAL potential1, potential2, log_potential;
   
   potential1 = ((REAL)count1/(REAL)totalcount1);

   
   potential2 = ((REAL)count2/(REAL)totalcount2);  
 
             
   log_potential = -log(potential1)-log(potential2);
   
   return(log_potential);
   

}
/***************************************/

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

/**************************************/
/* function that calculates distance (in angstroms)
   between two grid points */
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

/*************************************/
/* function that finds location of CA atom of a residue */
BOOL FindCAAtom(PDB *res1_start, PDB *res1_next, VEC3F *c_alpha)
{
   PDB *p;

   BOOL found = FALSE;
        
   for(p = res1_start; p !=res1_next; NEXT(p))
   {
      if(!strncmp(p->atnam, "CA", 2))
      {
         c_alpha->x = p->x;
         c_alpha->y = p->y;
         c_alpha->z = p->z;
         found = TRUE;
         break;
         
      }
      
   }
   return(found);
  
}
   
/*************************************/
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
            
/*************************************/
/* function that populates matrices from text file 
   (created in hydrogen_matrices.c program */
BOOL ReadInMatrices(char *res1, char *res2, FILE *matrix)
{
   char buffer[MAXBUFF], junk[15], minibuffer[6];
   int x, y, z, count;
   BOOL inResidue1, inResidue2,
        found_residue1 = FALSE, 
        found_residue2 = FALSE;
    
   while(fgets(buffer, MAXBUFF, matrix))
   {     
      TERMINATE(buffer);
      if(!strncmp(buffer, "residue", 7))
      {
         strcpy(minibuffer, buffer+8);
        
         if(!strncmp(minibuffer, res1, 3))         
         {
            inResidue1 = TRUE;
            
         }
         else
         {
            inResidue1 = FALSE;
         }
         if(!strncmp(minibuffer, res2, 3))
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
         /* storing location of where, if res 1 is a hydrogen acceptor,
            the *partner* donor atoms like to go */
         if(!strncmp(buffer, "donate", 6))
         {
            sscanf(buffer, "%s %d %d %d %d", junk, &x, &y, &z, &count);
            gDonate[x][y][z] = count; 
         }
         else if(!strncmp(buffer, "accept", 6))
         {
            sscanf(buffer, "%s %d %d %d %d", junk, &x, &y, &z, &count);
            gAccept[x][y][z] = count; 
         }
         found_residue1 = TRUE;
      }         
       if(inResidue2)
      {
         if(!strncmp(buffer, "partnertoaccept", 15))
         {
            /* storing res 2 donor atoms */
            sscanf(buffer, "%s %d %d %d %d", junk, &x, &y, &z, &count);          
            gPartnertoAccept[x][y][z] = count;
         }
         else if(!strncmp(buffer, "partnertodonate", 15))
         {
            /* storing res 2 donor atoms */
            sscanf(buffer, "%s %d %d %d %d", junk, &x, &y, &z, &count);          
            gPartnertoDonate[x][y][z] = count;
         }
         found_residue2 = TRUE; 
      }          
   }
   if(found_residue1 && found_residue2)
      return(TRUE);
   return(FALSE);
}

/********************************************************/
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
/**********************************************************/

/* function to open matrix file. Opens default matrix file
   if none specified on command line */

FILE *OpenMatrixFile(char *matrix_file)
{

   FILE *fp;

   char default_file[160];
   
   
   /* if filename has been specified on command line, open and return */
   if(matrix_file !=NULL && matrix_file[0])
   {      
      fp = fopen(matrix_file, "r");
        
   }
   /* if no filename specified, open default matrix file */
   else
   {
      strcpy(default_file, MATRIXFILE);     
      fp = fopen(default_file, "r");
     
   }
   return(fp);
   
      
}

/**********************************************************/

/* function to parse the command line */

BOOL ParseCmdLine(int argc, char **argv, REAL *cutoff, BOOL *hbplus,
                  char *hatom1, char *hatom2,char *matrix_file,char *pdbfile,
                  char *locres1, char *locres2, char *res1, char *res2, char *outputfile)
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
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* check there are 5 or 6 arguments remaining */
         if((argc > 6) || (argc < 5))
            return(FALSE);
         
         strcpy(pdbfile, argv[0]);
         argc--;
         argv++;
         strcpy(locres1, argv[0]);
         argc--;
         argv++;
         strcpy(locres2, argv[0]);
         argc--;
         argv++;
         strcpy(res1, argv[0]);
         argc--;
         argv++;
         strcpy(res2, argv[0]);
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
   return(TRUE);
}

/*************************************/

/* function to display a usage message */

void Usage(void)
{
   
   fprintf(stderr, "\nCheckHBond v1.0 (c) 2002, Alison Cuff,University of Reading\n\n");
   fprintf(stderr, "USAGE: checkhbond [-c cutoff] [-p hatom1 hatom2][-m matrix file]\n");
   fprintf(stderr, "   pdbfile residue1 residue2 nameres1 nameres2/sub [output file]\n\n");
   fprintf(stderr, "  -c [cutoff]: cutoff distance between hydrogen-capable atoms(default: 0.5A)\n");
   fprintf(stderr, "  -p: Parse HBplus data.\n");
   fprintf(stderr, "  Hydrogen donating atom (hatom1) and hydrogen accepting atom (hatom2) required \n");
   fprintf(stderr, "  -m [matrix file]: matrix file (if not using default file\n");
   fprintf(stderr, "  pdbfile:  pdb file of protein structure\n");
   fprintf(stderr, "  residue1:  First residue (chain, residue number, insert)\n");
   fprintf(stderr, "  residue2: Second residue (chain, residue number, insert)\n");
   fprintf(stderr, "  nameres1: Name of first residue\n");
   fprintf(stderr, "  nameres2/sub:  Name of replacement amino acid (mutation of residue 2),\n");
   fprintf(stderr, "      or native residue 1 if creating 'pseudo-energy distribution'\n");
   fprintf(stderr, "  [output file]: for saving hydrogen-capable atoms (in pdb format)\n");
   fprintf(stderr, "      I/O is through stdout if file not specified\n\n");   
   fprintf(stderr, "Determines the validity of a hydrogen bond in a protein \n");
   fprintf(stderr, "Assesses if hydrogen bond is maintained if one amino acid is substituted\n");
   fprintf(stderr, "for another and calculates a pseudoenergy for any that are maintained\n");
   fprintf(stderr, "With -p switch, creates 'reference' pseudoenergy values for real hydrogen bonds\n");
   fprintf(stderr, "against which the pseudoenergy of hydrogen bonds maintained after mutation can be\n");
   fprintf(stderr, "assessed for their quality\n\n");
   
}

/****************************************************/
      
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

/******************************************************/


BOOL IsHBondCapable(char *residue)
{   
   int i;
   
   char *NotHBondRes[3];
   
   BOOL flag = FALSE;
   
   NotHBondRes[0] = "MET";
   NotHBondRes[1] = "CYS";
   NotHBondRes[2] = "HIS";
   /* NotHBondRes[3] = "ILE";
   NotHBondRes[4] = "VAL";
   NotHBondRes[5] = "PHE";
   NotHBondRes[6] = "GLY";
   NotHBondRes[7] = "ALA";
   NotHBondRes[8] = "PRO";
   NotHBondRes[9] = "LEU";*/
   

   for(i = 0; i < 3; i++)
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

