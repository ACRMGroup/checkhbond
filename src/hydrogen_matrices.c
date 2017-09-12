/*************************************************************************

   Program:    hydrogen_matrices
   File:       hydrogen_matrices.c
   
   Version:    V2.1
   Date:       12.09.17
   Function:   Generate matrices of hydrogen bond information for use
               by checkhbond
   
   Copyright:  (c) University of Reading / Alison L. Cuff 2002-2017
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
creates four matrices ....

gAccept: stores all acceptor atoms of *key* residue
gDonate: stores all donor atoms of *key* residue
gPartnertoAccept: stores location of *partner* residue donor atoms 
                  i.e. those atoms that hydrogen bond with 
                  acceptor atoms in gAccept
gPartnertoDonate: stores location of *partner* residue acceptor atoms 
                  i.e. those atoms that hydrogen bond with 
                  acceptor atoms in gDonate
Input: cath list of protein structures
Output: text file printing matrices (see above) for hydrogen-capable
        residues in turn:

        for example:
        residue ASN
        gDonate (stores location of ASN hydrogen donor atoms)
        gAccept (stores location of ASN hydrogen acceptor atoms)
        gPartnertoDonor (stores location of *partner* acceptor atoms
                         that hydrogen-bond with ASN donor atoms)
        gPartnertoAccept (stores location of *partner* donor atoms 
                          that hydrogen-bond with ASN acceptor atoms)

output file to be used in program checkhbond.c

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  17.06.03 Original version as used for thesis
   V1.1  19.08.05 Various bug fixes By: ACRM
   V2.0  24.01.05 Modified to allow mc/sc matrices to be generated
   V2.1  12.09.17 Updated for new Bioplib and some cleanup

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include <math.h>
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/matrix.h"
#include "bioplib/hbond.h"
#include "orientate.h"
#include "hbondmat2.h"
#include "cavallo_userfunc.h"

/************************************************************************/
/* Defines and macros
*/
#define NOISY

struct hbond_data
{
   char *residue;
   char *hatom;
   BOOL select;
   BOOL accept;
   BOOL donate;
   char *s_donh1;
   char *s_donh2;
   char *s_donh3;
   char *s_ante;
   struct hbond_data *next;
   
};

typedef struct hbond_data HBOND;

struct pdb_names
{
   char filename[8];
   struct pdb_names *next;
};

typedef struct pdb_names NAMES;

/************************************************************************/
/* Globals
*/
/* matrices that store all acceptor and donor atoms */
static int gAccept[MAXSIZE][MAXSIZE][MAXSIZE];
static int gDonate[MAXSIZE][MAXSIZE][MAXSIZE];

/* matrix storing partner atoms to hydrogen accepting atoms */
static int gPartnertoAccept[MAXSIZE][MAXSIZE][MAXSIZE];
/* matrix storing partner atoms to hydrogen donating heavy atoms */
static int gPartnertoDonate[MAXSIZE][MAXSIZE][MAXSIZE];

/************************************************************************/
/* Prototypes
*/
int main (int argc, char *argv[]);
HBOND *InitializeHbondTypes(void);
BOOL ParseCmdLine(int argc, char **argv, char *inputfile, char *outputfile);
void Usage(void);
NAMES *InitializeDomainList(FILE *fp);
char *FindStructureLocation(NAMES *names, BOOL *tempflag);
BOOL CalcAndStoreHBondData(HBOND *hb, NAMES *names, FILE *out);
BOOL isDonor(PDB *d, HBOND *hb, char *h_name1, char *h_name2, char *h_name3);
BOOL isAcceptor(PDB *a, HBOND *hb, char *p_name);
PDB *FindAtomInRange(PDB *start, PDB *stop, char *name);
void FindHAtoms(PDB *resA, PDB *stopA, PDB *resB, PDB *stopB, HBOND *hb);
void ClearArrays(void);
void StoreHBondingPosition(PDB *start, PDB *next, HBOND *hb);
void PrintMatrix(HBOND *h, FILE *out);
void StorePartnertoDonatePosition(PDB *a);
void StorePartnertoAcceptPosition(PDB *d);
void FindMCDonorHAtoms(PDB *resA, PDB *stopA, PDB *resB, PDB *stopB, HBOND *hb);
void StoreHBondingNPosition(PDB *start, PDB *stop, HBOND *hb);
void FindMCAcceptorAtoms(PDB *resA, PDB *stopA, PDB *resB, PDB *stopB, HBOND *hb);
void StoreHBondingCOPosition(PDB *start, PDB *stop, HBOND *hb);
void FindHAtomsSCMC(PDB *resA, PDB *stopA, PDB *resB, PDB *stopB, HBOND *hb);





/************************************************************************/
int main (int argc, char *argv[])
{
   FILE *in = stdin, *out = stdout;
   HBOND *hb;
   char inputfile[160], outputfile[160];
   NAMES *names;
   
   inputfile[0] = outputfile[0] = '\0';
   
   if(ParseCmdLine(argc, argv, inputfile, outputfile))
   {
      if(blOpenStdFiles(inputfile, outputfile, &in, &out))
      {
         if((hb = InitializeHbondTypes()) !=NULL)
         {
            if((names = InitializeDomainList(in)) !=NULL)
            {
               if(!CalcAndStoreHBondData(hb, names, out))
               {
                  return(1);
               }
            }
            else
            {
               printf("ERROR: Can't read protein domains linked list\n");
               printf("       (perhaps no high resolution structures?)\n");
               return(1);
            }
         }
         else
         {
            printf("ERROR: Internal error in building hydrogen atoms linked list\n");
            return(1);
         }
      }
      else
      {
         printf("ERROR: Unable to open specified input or output files\n");
         return(1);
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}
   
               
/************************************************************************/
BOOL CalcAndStoreHBondData(HBOND *hb, NAMES *names, FILE *out)
{
   NAMES *n;
   HBOND *h;
   int natoms, natoms2, nHatoms;
   PDB *pdb, *pdb2, *start, *next, *nextres, *stop, *prev;
   FILE *fp1 = NULL,
        *fp2 = NULL;
   char *location;
   BOOL noenv, tempflag;
 
   if((fp1 = blOpenFile(PGPFILE, "DATADIR", "r", &noenv)) == NULL)
   {
      fprintf(stderr, "ERROR: Can't open pgp file\n");
      if(noenv)
      {
         fprintf(stderr, "       DATADIR environment variable not set\n");
      }
      
      return(FALSE);
   }
   
   for(h = hb; h!=NULL; NEXT(h))
   {
#if defined(MCDONOR) || defined(MCACCEPTOR)
      if(h->select)
#else
      if(h->select && (h->accept || h->donate))
#endif
      {
         fprintf(stderr,"INFO: Processing residue type %s\n",h->residue);

         ClearArrays();
         
         for(n=names; n !=NULL; NEXT(n))
         {
            if((location = FindStructureLocation(n, &tempflag)) == NULL)
            {
               /* 19.08.05 ACRM: Corrected from 'location' to 'n->filename' 
                  Also, FindStructureLocation() will generate a warning, so
                  this is a continuation.
                */
               fprintf(stderr,"         File not processed for: %s\n", n->filename);
            }
            else
            {
#ifdef NOISY
               fprintf(stderr,"INFO: Processing file %s\n", location);
#endif
               /* open protein domain file */
               if((fp2 = fopen(location, "r")) !=NULL)
               { 
                  /* create linked list of pdb file */  
                  if((pdb = blReadPDBAtoms(fp2, &natoms)) !=NULL)
                  {   
                     /* strip any hydrogens present in protein domain file */
                     if((pdb2 = blStripHPDBAsCopy(pdb, &natoms2)) !=NULL)
                     {
                        FREELIST(pdb, PDB);
                        pdb  = pdb2;
                        prev = NULL;
                        
                        if((nHatoms = blHAddPDB(fp1, pdb)) !=0)                    
                        {
                           for(start=pdb; start!=NULL; prev=start, start=next)
                           {
                              next = blFindNextResidue(start);

#if defined(MCDONOR)
                              /* If not proline and not first residue */
                              if(!strncmp(start->record_type, "ATOM  ", 6) &&
                                 !strncmp(start->resnam, h->residue, 3) &&
                                 strncmp(start->resnam, "PRO", 3) && 
                                 (prev != NULL) &&
                                 ResiduesBonded(prev, start))
                              {
                                 /* return true if backbone atoms cannot be found
                                    .. stops program
                                    progressing to next stage 

                                    This key residue is the one NOT being mutated.
                                    So if we are doing sc/mc HBonds, this is always
                                    the m/c residue
                                 */
                                 if((OrientateN_PDB(pdb, prev, start, next) != FALSE))
                                 {
                                    StoreHBondingNPosition(start, next, hb);
                                    
                                    for(nextres=pdb; nextres!=NULL; nextres=stop)
                                    {
                                       stop = blFindNextResidue(nextres);
                                       
                                       if((nextres !=start))
                                       {
                                          /*
                                            printf("%s %s\n", start->resnam, 
                                            nextres->resnam);
                                          */
                                          
                                          if((blIsMCDonorHBonded(start,nextres,HBOND_SIDE2))
                                             !=0)
                                          {
                                             FindMCDonorHAtoms(start, next, nextres, 
                                                               stop, hb);
                                          }
                                       }
                                    }
                                 }
                                 else /* If we couldn't orientate the residue */
                                 {
                                    printf("WARNING: backbone atoms can't be found for %c%d%c\
 (PDB file: %s)\n",
                                           start->chain[0], start->resnum, start->insert[0],
                                           location);
                                 }
                              }  /* If the residue type matches */
#elif defined(MCACCEPTOR)
                              if(!strncmp(start->record_type, "ATOM  ", 6) &&
                                 !strncmp(start->resnam, h->residue, 3))
                              {
                                 /* return true if backbone atoms cannot be found
                                    .. stops program
                                    progressing to next stage 
                                 */
                                 if((OrientateCO_PDB(pdb, start, next) != FALSE))
                                 {
                                    StoreHBondingCOPosition(start, next, hb);
                                    
                                    for(nextres=pdb; nextres!=NULL; nextres=stop)
                                    {
                                       stop = blFindNextResidue(nextres);
                                       
                                       if((nextres !=start))
                                       {
                                          /*
                                            printf("%s %s\n", start->resnam, 
                                            nextres->resnam);
                                          */
                                          
                                          if((blIsMCAcceptorHBonded(start,nextres,HBOND_SIDE2))
                                             !=0)
                                          {
                                             FindMCAcceptorAtoms(start, next, nextres, 
                                                                 stop, hb);
                                          }
                                       }
                                    }
                                 }
                                 else /* If we couldn't orientate the residue */
                                 {
                                    printf("WARNING: backbone atoms can't be found for %c%d%c\
 (PDB file: %s)\n",
                                           start->chain[0], start->resnum, start->insert[0],
                                           location);
                                 }
                              }
#else /* SCSC */
                              if((!strncmp(start->resnam, h->residue, 3)))
                              {
                                 /* return true if backbone atoms cannot be found
                                    .. stops program
                                    progressing to next stage 
                                 */
                                 if((OrientatePDB(pdb, start, next) != FALSE))
                                 {
                                    StoreHBondingPosition(start, next, hb);
                                    
                                    for(nextres=pdb; nextres!=NULL; nextres=stop)
                                    {
                                       stop = blFindNextResidue(nextres);
                                       
                                       if((nextres !=start))
                                       {
                                          /*
                                            printf("%s %s\n", start->resnam, 
                                            nextres->resnam);
                                          */
                                          
#ifdef SCMC
                                          if((blIsHBonded(start,nextres,HBOND_SIDECHAIN))
                                             !=0)
                                          {
                                             FindHAtomsSCMC(start, next, nextres, 
                                                            stop, hb);
                                          }
#else
                                          if((blIsHBonded(start,nextres,HBOND_SS))
                                             !=0)
                                          {
                                             FindHAtoms(start, next, nextres, 
                                                        stop, hb);
                                          }
#endif
                                       }
                                    }
                                 }
                                 else /* If we couldn't orientate the residue */
                                 {
                                    printf("WARNING: backbone atoms can't be found for %c%d%c\
 (PDB file: %s)\n",
                                           start->chain[0], start->resnum, start->insert[0],
                                           location);
                                 }
                              }  /* If the residue type matches */
#endif
                           }  /* For each residue in the PDB */
                        }  /* If we added the hydrogens */
                     }  /* If we stripped the hydrogens */
                  }
                  else /* if we didn't read the PDB file */
                  {
                     printf("WARNING: Can't read atom list from PDB file %s\n",
                            location);
                  }
                  if(pdb != NULL) FREELIST(pdb, PDB); 
               }
            }
            if(fp2 !=NULL) fclose(fp2);
            /* ACRM 19.08.05 - Moved this in here. Previously only one
               temporary file was created and re-used each time. After
               Antonio's changes, a new filename was created for each
               temp file being processed, so they weren't getting
               deleted. Similarly, the location variable needs to
               be freed on each cycle.
            */
            /* remove protein domain files created by getchain */
            if(tempflag)
            {
               if(remove(location))
               {
                  fprintf(stderr, "WARNING: Unable to remove temporary protein \
domain file %s\n", location);
               } 
            }
            free(location);
         }  /* foreach name, n */
         PrintMatrix(h, out);
      }  /* if(h->select) */
   }  /* For each hydrogen bond type, h */
   
   fclose(fp1);
   return(TRUE);
}

/************************************************************************/
void FindHAtoms(PDB *resA, PDB *stopA, PDB *resB, PDB *stopB, HBOND *hb)
{
   PDB *d, *a, *h1, *h2, *h3, *p;
   
   char h_name1[6], h_name2[6], h_name3[6], p_name[6];

   for(d=resA; d!=stopA; NEXT(d))
   {
      /* if residue A heavy atom is a  donor */
      if(isDonor(d, hb, h_name1, h_name2, h_name3))
      {
         for(a=resB; a!=stopB; NEXT(a))
         {
            /* and  residue B is a acceptor */
            if(isAcceptor(a, hb, p_name))
            {
               h1 = FindAtomInRange(resA, stopA, h_name1);
               h2 = FindAtomInRange(resA, stopA, h_name2);
               h3 = FindAtomInRange(resA, stopA, h_name3);
              
               p = FindAtomInRange(resB, stopB, p_name);
               
               /* if the hydrogen bond is valid */
               if(blValidHBond(h1, d, a, p))
               {
                  /* store position of partner acceptor atom */          
#ifdef NOISY
                  fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                          d->resnam, d->resnum, d->atnam, 
                          a->resnam, a->resnum,  a->atnam, 
                          ((h1==NULL?"???":h1->atnam)), p->atnam); 
#endif
                  StorePartnertoDonatePosition(a);
               }
               else if(h2!=NULL)
               {
                  if(blValidHBond(h2, d, a, p))
                  {                     
#ifdef NOISY
                     fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                            d->resnam, d->resnum, d->atnam, 
                            a->resnam, a->resnum, a->atnam, 
                            h2->atnam, p->atnam);
#endif
                     StorePartnertoDonatePosition(a);
                  }
               }
               else if(h3!=NULL)
               {
                  if(blValidHBond(h3, d, a, p))
                  {
#ifdef NOISY                     
                     fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                            d->resnam, d->resnum, d->atnam, 
                            a->resnam, a->resnum, a->atnam, 
                            h3->atnam, p->atnam);
#endif
                     StorePartnertoDonatePosition(a);
                  }
               }
            }
         }
      }
   }

   /* if residue A is an acceptor */
   for(a=resA; a!=stopA; NEXT(a))
   {
      if(isAcceptor(a, hb, p_name))    
      {
         /* and residue B is a donor */

         for(d=resB; d!=stopB; NEXT(d))
         {
            if(isDonor(d, hb, h_name1, h_name2, h_name3))
            {
               h1 = FindAtomInRange(resB, stopB, h_name1);
               h2 = FindAtomInRange(resB, stopB, h_name2);
               h3 = FindAtomInRange(resB, stopB, h_name3);
               
               p = FindAtomInRange(resA, stopA, p_name);
               
               /* and the hydrogen bond is valid */

               if(blValidHBond(h1, d, a, p))
               {
                  /* store position of partners hydrogen donating heavy atom */
#ifdef NOISY
                  fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                          d->resnam, d->resnum, d->atnam, 
                          a->resnam, a->resnum, a->atnam, 
                          ((h1==NULL?"???":h1->atnam)), p->atnam);
#endif
                  StorePartnertoAcceptPosition(d);
               }
               else if(h2!=NULL)
               {
                  if(blValidHBond(h2, d, a, p))
                  {
#ifdef NOISY
                     fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                             d->resnam, d->resnum, d->atnam, 
                             a->resnam, a->resnum, a->atnam, 
                             h2->atnam, p->atnam);
#endif
                     StorePartnertoAcceptPosition(d);
                  }
               }
               else if(h3!=NULL)
               {
                  if(blValidHBond(h3, d, a, p))
                  {
#ifdef NOISY
                     fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                            d->resnam, d->resnum, d->atnam, 
                            a->resnam, a->resnum, a->atnam, 
                            h3->atnam, p->atnam);
#endif
                     StorePartnertoAcceptPosition(d);
                  }
               }
            }
         }
      }
   }
}

/************************************************************************/
BOOL isDonor(PDB *d, HBOND *hb, char *h_name1, char *h_name2, char *h_name3)
{
   HBOND *h;
   
   for(h=hb; h!=NULL; NEXT(h))
   {
      if(h->donate && !strncmp(h->residue, d->resnam, 3) &&
         !strncmp(h->hatom, d->atnam, 3))
      {
         if(h->s_donh1 !=NULL)
            strcpy(h_name1, h->s_donh1);
         else
            h_name1[0] = '\0';
         if(h->s_donh2 !=NULL)
            strcpy(h_name2, h->s_donh2);
         else
            h_name2[0] = '\0';
         if(h->s_donh3 !=NULL)
            strcpy(h_name3, h->s_donh3);
         else
            h_name3[0] = '\0';

         return(TRUE);
      }
   }
   h_name1[0] = '\0';
   return(FALSE);
}

/************************************************************************/
BOOL isAcceptor(PDB *a, HBOND *hb, char *p_name)
{
   HBOND *h;
   
   for(h=hb; h!=NULL; NEXT(h))
   {
      if(h->accept && !strncmp(h->residue, a->resnam, 3) &&
         !strncmp(h->hatom, a->atnam, 3))
      {
         if(h->s_ante !=NULL)
            strcpy(p_name, h->s_ante);
         else
            p_name[0] = '\0';

         return(TRUE);
      }
   }
   p_name[0] = '\0';
   return(FALSE);
}

/************************************************************************/
PDB *FindAtomInRange(PDB *start, PDB *stop, char *name)
{
   PDB *p;
   
   if(!name[0])
      return(NULL);
   
   for(p=start; p!=stop; NEXT(p))
   {
      if(!strncmp(p->atnam, name, 4))
      {
         return(p);
      }
   }
   return(NULL);
}

/************************************************************************/
void StoreHBondingPosition(PDB *start, PDB *stop, HBOND *hb)
{
   PDB *p;
   
   HBOND *h;
     
   for(h=hb; h!=NULL; NEXT(h))
   {   
      if(!strncmp(start->resnam, h->residue, 3))
      {
         for(p = start; p!=stop; NEXT(p))
         {
            /* search for hydrogen-capable /residue/atoms listed in 
               *hydrogen residue* linked list 
             */
            if(!strncmp(p->atnam, h->hatom, 3))
            {
               /* if hydrogen-capable atom found, is it a donor */ 
               if(h->donate)
               {
                  /* store location */
                  gDonate[(int)(p->x/DIV) + OFFSET] 
                     [(int)(p->y/DIV) + OFFSET]
                     [(int)(p->z/DIV) + OFFSET]++;
               }

               /* or acceptor */
               if(h->accept)
               {
                  /* store location */
                  gAccept[(int)(p->x/DIV) + OFFSET]  
                     [(int)(p->y/DIV) + OFFSET]
                     [(int)(p->z/DIV) + OFFSET]++;
               }
               break;
            }
         }
      }
   }

   /* ACRM 18.08.05 removed this! Luckily it did no harm as
      'h' would always be NULL when we got this far
   */
/*   FREELIST(h, HBOND); */
}

/************************************************************************/
void StorePartnertoAcceptPosition(PDB *d)
{
   gPartnertoAccept[(int)(d->x/DIV) + OFFSET]
      [(int)(d->y/DIV) + OFFSET]
      [(int)(d->z/DIV) + OFFSET]++;   
}

/************************************************************************/
void StorePartnertoDonatePosition(PDB *a)
{
   gPartnertoDonate[(int)(a->x/DIV) + OFFSET]
      [(int)(a->y/DIV) + OFFSET]
      [(int)(a->z/DIV) + OFFSET]++;   
}

/************************************************************************/
void PrintMatrix(HBOND *h, FILE *out)
{
   int x, y, z;

   if(h->select)
   {
      fprintf(out, "residue %s\n", h->residue);
   }
      
   for(x = 0; x < MAXSIZE; x++)
   {
      for(y = 0; y < MAXSIZE; y++)
      {
         for(z = 0; z < MAXSIZE; z++)
         {
            if(gDonate[x][y][z] > 0)
            {
               fprintf(out, "donate\t%8d\t%8d\t%8d\t%6d\n",  x, y, z,  
                       gDonate[x][y][z]);

               /*
                 GRID_2_COORD(x, grid.x);
                 GRID_2_COORD(y, grid.y);
                 GRID_2_COORD(z, grid.z);
               
                 fprintf(out, "ATOM  %5d  CA  THR  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                 atnum++, resnum++,  grid.x, grid.y, grid.z, 1.00, 1.00);
               */
            }
         }
      }
   }
   /* printing out matrix containing partner atoms to hydrogen-donating 
      heavy atoms 
   */     
   for(x = 0; x < MAXSIZE; x++)
   {
      for(y = 0; y < MAXSIZE; y++)
      {
         for(z = 0; z < MAXSIZE; z++)
         {
            if(gPartnertoDonate[x][y][z] > 0)
            {
               fprintf(out, "partnertodonate\t%8d\t%8d\t%8d\t%6d\n", 
                       x, y, z,  gPartnertoDonate[x][y][z]);

               /* 
                  GRID_2_COORD(x, grid.x);
                  GRID_2_COORD(y, grid.y);
                  GRID_2_COORD(z, grid.z);
              
                  fprintf(out, "ATOM  %5d  CA  THR  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                  atnum++, resnum++,  grid.x, grid.y, grid.z, 1.00, 2.00);
               */
            }
         }
      }
   }
         
   for(x = 0; x < MAXSIZE; x++)
   {
      for(y = 0; y < MAXSIZE; y++)
      {
         for(z = 0; z < MAXSIZE; z++)
         {
            if(gAccept[x][y][z] > 0)
            {
               fprintf(out, "accept\t%8d\t%8d\t%8d\t%6d\n",   
                       x, y, z,  gAccept[x][y][z]);

               /*GRID_2_COORD(x, grid.x);
               GRID_2_COORD(y, grid.y);
               GRID_2_COORD(z, grid.z);
               
               fprintf(out, "ATOM  %5d  CA  THR  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
               atnum++, resnum++,  grid.x, grid.y, grid.z, 1.00, 3.00);*/
               
            }
         }
      }
   }
   
   /* printing out matrix containing partner atoms to hydrogen 
      accepting atoms 
   */
   for(x = 0; x < MAXSIZE; x++)
   {
      for(y = 0; y < MAXSIZE; y++)
      {
         for(z = 0; z < MAXSIZE; z++)
         {
            if(gPartnertoAccept[x][y][z] > 0)
            {
              fprintf(out, "partnertoaccept\t%8d\t%8d\t%8d\t%6d\n", 
                      x, y, z,  gPartnertoAccept[x][y][z]);

               /*GRID_2_COORD(x, grid.x);
               GRID_2_COORD(y, grid.y);
               GRID_2_COORD(z, grid.z);
               
               fprintf(out, "ATOM  %5d  CA  THR  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
               atnum++, resnum++,  grid.x, grid.y, grid.z, 1.00, 4.00);*/
            }
         }
      }
   }
}

/************************************************************************/
void ClearArrays()
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
/* function to create linked list of protein domain files */
NAMES *InitializeDomainList(FILE *fp)
{
   NAMES *names = NULL;
   
   NAMES *n;
   FILE *fn = fp;
   char *p;
   char buffer[MAXBUFF];
   float resolution;
   
   while(fgets(buffer, MAXBUFF, fn))
   {
      if(buffer[0] == '#')
         continue;
      
      TERMINATE(buffer);
      /* ACRM 18.08.05 Changed to use sscanf as the field widths
         in the CATH file have changed over time... (Also was
         really 54 not 53!)
      */
/*
      strcpy(resol, buffer+53);
      resolution = atof(resol);
*/
      sscanf(buffer, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %f", &resolution);
    
      if(resolution <= RESOL)
      {
         if((p = strchr(buffer, ' '))!=NULL)
            *p='\0';
         
         if(names == NULL)
         {
            INIT(names, NAMES);
            n = names;
         }
         else
         {
            ALLOCNEXT(n, NAMES);
         }
         if(n==NULL)
         {
            return(NULL);
         }
         strncpy(n->filename, buffer, 7);
      }
      else
      {
         buffer[7] = '\0';
         fprintf(stderr,"INFO: %s Discarded as resolution is %.2f\n",
                 buffer, resolution);
      }
   }

   return(names);
}

/************************************************************************/
char *FindStructureLocation (NAMES *names, BOOL *tempflag)
{
   char *filename=(char *)malloc((size_t)(sizeof(char)*(6+1)));
   char *result,*tbuffer,*tbuffer1,chainlabel;
   *tempflag = 0;
   
   if(!filename)
   {
      fprintf(stderr, "WARNING: Unable to allocate memory for name of \
structure file\n");
      return(NULL);
      
   }
   
   if(strlen(names->filename) !=6) 
   {
      fprintf(stderr,"WARNING: Name from CATH file does not contain 6 char: %s\n",
              names->filename);
      return(NULL);
   } 

   filename=strcpy(filename,names->filename);
   if(filename[5]!='0')
   {
      /*dompdb */
      result=multiappend("%s%s",DOMAINLOC,filename);		  
      free(filename);
      return(result);
   }
   else
   {
      if(filename[4]!='0')
      {
         /* getchain */
         chainlabel=filename[4];
         filename[4]='\0';
         tbuffer=multiappend("getchain %c %s%s%s",
                             chainlabel,PDBLOC,filename,PDBEXT);
         tbuffer1=multiappend("%s/%s",
                              tmpdir(),filename);
         doIt(tbuffer,tbuffer1);
         free(tbuffer);
         free(filename);
         *tempflag = 1;
         return(tbuffer1);
      }
      else
      {
         /* pdb */
         filename[4]='\0';
         result=multiappend("%s%s%s",PDBLOC,filename,PDBEXT);		  
         free(filename);
         return(result);
      }
   }	  

   fprintf(stderr,"WARNING: Internal error in finding file\n");
   return(NULL);
}

/************************************************************************/
/* function that creates a linked list of all the hydrogen-capable 
   atoms present in each amino acid 
*/
HBOND *InitializeHbondTypes()
{
   /* ACRM 18.08.05 Fixed donh1 in h17 and h18 to add the missing trailing space
      Tidied the columns so we can see what is going on!
      
      Also changed the following entries since we don't know precisely where the
      hydrogens go:

   static HBOND h2  = {"ARG", "NH1 ", 0,  0,  1,  "HH11", "HH12", NULL,   NULL,  NULL};
   static HBOND h3  = {"ARG", "NH2 ", 0,  0,  1,  "HH21", "HH22", NULL,   NULL,  NULL};
   static HBOND h5  = {"ASN", "ND2 ", 1,  0,  1,  "HD21", "HD22", NULL,   NULL,  NULL};
   static HBOND h11 = {"GLN", "NE2 ", 1,  0,  1,  "HE21", "HE22", NULL,   NULL,  NULL};
   static HBOND h13 = {"LYS", "NZ  ", 1,  0,  1,  "HZ1 ", "HZ2 ", "HZ3 ", NULL,  NULL};
   */

   /*                  res    hatom   sel acc don donh1   donh2   donh3   ante   next*/
   static HBOND h1  = {"ARG", "NE  ", 1,  0,  1,  "HE  ", NULL,   NULL,   NULL,  NULL};
   static HBOND h2  = {"ARG", "NH1 ", 0,  0,  1,  NULL,   NULL,   NULL,   NULL,  NULL};
   static HBOND h3  = {"ARG", "NH2 ", 0,  0,  1,  NULL,   NULL,   NULL,   NULL,  NULL};
   static HBOND h4  = {"THR", "OG1 ", 1,  1,  1,  NULL,   NULL,   NULL,   "CB  ",NULL};
   static HBOND h5  = {"ASN", "ND2 ", 1,  0,  1,  NULL,   NULL,   NULL,   "CG  ",NULL}; /* Acceptor? */
   static HBOND h6  = {"ASN", "OD1 ", 0,  1,  0,  NULL,   NULL,   NULL,   "CG  ",NULL};
   static HBOND h7  = {"ASP", "OD1 ", 1,  1,  0,  NULL,   NULL,   NULL,   "CG  ",NULL};
   static HBOND h8  = {"ASP", "OD2 ", 0,  1,  0,  NULL,   NULL,   NULL,   "CG  ",NULL};
   static HBOND h9  = {"GLU", "OE1 ", 1,  1,  0,  NULL,   NULL,   NULL,   "CD  ",NULL};
   static HBOND h10 = {"GLU", "OE2 ", 0,  1,  0,  NULL,   NULL,   NULL,   "CD  ",NULL};
   static HBOND h11 = {"GLN", "NE2 ", 1,  0,  1,  NULL,   NULL,   NULL,   "CD  ",NULL}; /* Acceptor? */
   static HBOND h12 = {"GLN", "OE1 ", 0,  1,  0,  NULL,   NULL,   NULL,   "CD  ",NULL};
   static HBOND h13 = {"LYS", "NZ  ", 1,  0,  1,  NULL,   NULL,   NULL,   NULL,  NULL};
   static HBOND h14 = {"SER", "OG  ", 1,  1,  1,  NULL,   NULL,   NULL,   "CB  ",NULL};
   static HBOND h15 = {"TRP", "NE1 ", 1,  0,  1,  "HE1 ", NULL,   NULL,   NULL,  NULL}; /* Acceptor? */
   static HBOND h16 = {"TYR", "OH  ", 1,  1,  1,  NULL,   NULL,   NULL,   "CZ  ",NULL};
   static HBOND h17 = {"HIS", "ND1 ", 1,  1,  1,  "HD1 ", NULL,   NULL,   "CG  ",NULL};
   static HBOND h18 = {"HIS", "NE1 ", 1,  1,  1,  "HE1 ", NULL,   NULL,   "CD2 ",NULL};
   /* Backbone only HBonds */
   static HBOND h19 = {"ALA", "    ", 1,  0,  0,  NULL,   NULL,   NULL,   NULL,  NULL};
   static HBOND h20 = {"CYS", "    ", 1,  0,  0,  NULL,   NULL,   NULL,   NULL,  NULL};
   static HBOND h21 = {"PHE", "    ", 1,  0,  0,  NULL,   NULL,   NULL,   NULL,  NULL};
   static HBOND h22 = {"GLY", "    ", 1,  0,  0,  NULL,   NULL,   NULL,   NULL,  NULL};
   static HBOND h23 = {"ILE", "    ", 1,  0,  0,  NULL,   NULL,   NULL,   NULL,  NULL};
   static HBOND h24 = {"LEU", "    ", 1,  0,  0,  NULL,   NULL,   NULL,   NULL,  NULL};
   static HBOND h25 = {"MET", "    ", 1,  0,  0,  NULL,   NULL,   NULL,   NULL,  NULL};
   static HBOND h26 = {"VAL", "    ", 1,  0,  0,  NULL,   NULL,   NULL,   NULL,  NULL};
   static HBOND h27 = {"PRO", "    ", 1,  0,  0,  NULL,   NULL,   NULL,   NULL,  NULL};



   HBOND *first;
   
   first = &h1; 
   h1.next = &h2;
   h2.next = &h3;
   h3.next = &h4;
   h4.next = &h5;
   h5.next = &h6;
   h6.next = &h7;
   h7.next = &h8;
   h8.next = &h9;
   h9.next = &h10;
   h10.next = &h11;
   h11.next = &h12;
   h12.next = &h13;
   h13.next = &h14;
   h14.next = &h15;
   h15.next = &h16;
   h16.next = &h17;
   h17.next = &h18;
   h18.next = &h19;
   h19.next = &h20;
   h20.next = &h21;
   h21.next = &h22;
   h22.next = &h23;
   h23.next = &h24;
   h24.next = &h25;
   h25.next = &h26;
   h26.next = &h27;
   h27.next = NULL;
   
   return(first);
}

/************************************************************************/
/* function to display a usage message */
void Usage(void)
{
   fprintf(stderr, "\nHydrogen Matrices V2.0 (c) 2002-6, Alison Cuff, University of Reading\n");
   fprintf(stderr, "V1.1/2.0 modifications, Andrew C.R. Martin, University College London\n\n");
   
   fprintf(stderr, "Usage: hydrogen_matrices [cath domain file] [output file]\n\n");
   fprintf(stderr, "  [cath domain file] non-redundant (e.g Sreps) cath domain list file\n");
   fprintf(stderr, "  [output file] name of file to print out matrices\n");
   fprintf(stderr, "                I/O is though stdout if file not specified\n\n");   
   fprintf(stderr, "Takes a list of protein structures, and for each hydrogen-bonding residue\n");   
   fprintf(stderr, "orientates protein so that the carbon alpha atom of the residue is at the\n");
   fprintf(stderr, "origin, the carbon beta atom on the xy plane, and the nitrogen atom on the\n");
   fprintf(stderr, "x axis. Program then determines and prints out:\n");
   fprintf(stderr,"  (1) location of hydrogen donor atoms\n");
   fprintf(stderr,"  (2) location of hydrogen  acceptor atoms\n");
   fprintf(stderr,"  (3) location of *partner* acceptor atoms to each donor atom\n");
   fprintf(stderr,"  (4) location of *partner* donor atoms to each acceptor atom\n");
   fprintf(stderr,"for that residue\n\n");
    
}
      
/************************************************************************/
/* function to parse the command line */
BOOL ParseCmdLine(int argc, char **argv, char *inputfile, char *outputfile)
{
   argc--;
   argv++;
        
   if(argc > 2 || argc < 1)
   {
      return(FALSE);
   }
      
   strcpy(inputfile, argv[0]);
   argc--;
   argv++;      
   if(argc)
   {
      strcpy(outputfile, argv[0]);
   }
   return(TRUE);
   
   argc--;
   argv++;
      
   return(TRUE);
}

/************************************************************************/
void StoreHBondingCOPosition(PDB *start, PDB *stop, HBOND *hb)
{
   PDB *p;
   
   for(p = start; p!=stop; NEXT(p))
   {
      /* search for hydrogen-capable /residue/atoms listed in 
       *hydrogen residue* linked list 
       */
      if(!strncmp(p->atnam, "O  ", 3))
      {
         /* store location */
         gAccept[(int)(p->x/DIV) + OFFSET]  
            [(int)(p->y/DIV) + OFFSET]
            [(int)(p->z/DIV) + OFFSET]++;
         break;
      }
   }
}

/************************************************************************/
void FindMCAcceptorAtoms(PDB *resA, PDB *stopA, PDB *resB, PDB *stopB, HBOND *hb)
{
   PDB *d, *a, *h1, *h2, *h3, *p;
   
   char h_name1[6], h_name2[6], h_name3[6], p_name[6];

   strcpy(p_name, "C   ");

   /* Residue A *is* an acceptor (backbone oxygen) */
   a = FindAtomInRange(resA, stopA, "O   ");
   p = FindAtomInRange(resA, stopA, "C   ");

   /* If residue B is a donor */
   for(d=resB; d!=stopB; NEXT(d))
   {
      if(isDonor(d, hb, h_name1, h_name2, h_name3))
      {
         h1 = FindAtomInRange(resB, stopB, h_name1);
         h2 = FindAtomInRange(resB, stopB, h_name2);
         h3 = FindAtomInRange(resB, stopB, h_name3);
               
         /* and the hydrogen bond is valid */
         if(blValidHBond(h1, d, a, p))
         {
            /* store position of partners hydrogen donating heavy atom */
#ifdef NOISY
            fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                    d->resnam, d->resnum, d->atnam, 
                    a->resnam, a->resnum, a->atnam, 
                    ((h1==NULL?"???":h1->atnam)), p->atnam);
#endif
            StorePartnertoAcceptPosition(d);
         }
         else if(h2!=NULL)
         {
            if(blValidHBond(h2, d, a, p))
            {
#ifdef NOISY
               fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                       d->resnam, d->resnum, d->atnam, 
                       a->resnam, a->resnum, a->atnam, 
                       h2->atnam, p->atnam);
#endif
               StorePartnertoAcceptPosition(d);
            }
         }
         else if(h3!=NULL)
         {
            if(blValidHBond(h3, d, a, p))
            {
#ifdef NOISY
               fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                       d->resnam, d->resnum, d->atnam, 
                       a->resnam, a->resnum, a->atnam, 
                       h3->atnam, p->atnam);
#endif
               StorePartnertoAcceptPosition(d);
            }
         }
      }
   }
}

/************************************************************************/
void StoreHBondingNPosition(PDB *start, PDB *stop, HBOND *hb)
{
   PDB *p;
   
   /* If it's not Proline */
   if(strncmp(start->resnam, "PRO", 3))
   {
      for(p = start; p!=stop; NEXT(p))
      {
         /* search for backbone nitrogen */
         if(!strncmp(p->atnam, "N  ", 3))
         {
            /* store location */
            gDonate[(int)(p->x/DIV) + OFFSET] 
                   [(int)(p->y/DIV) + OFFSET]
                   [(int)(p->z/DIV) + OFFSET]++;
         }
      }
   }
}


/************************************************************************/
void FindMCDonorHAtoms(PDB *resA, PDB *stopA, PDB *resB, PDB *stopB, HBOND *hb)
{
   PDB *d, *a, *h1, *p;
   
   char p_name[6];

   for(d=resA; d!=stopA; NEXT(d))
   {
      if(!strcmp(d->atnam, "N   "))
      {
         h1 = FindAtomInRange(resA, stopA, "H   ");

         for(a=resB; a!=stopB; NEXT(a))
         {
            if(isAcceptor(a, hb, p_name))
            {
               p = FindAtomInRange(resB, stopB, p_name);
               
               /* if the hydrogen bond is valid */
               if(blValidHBond(h1, d, a, p))
               {
                  /* store position of partner acceptor atom */          
#ifdef NOISY
                  fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                          d->resnam, d->resnum, d->atnam, 
                          a->resnam, a->resnum,  a->atnam, 
                          ((h1==NULL?"???":h1->atnam)), p->atnam); 
#endif
                  StorePartnertoDonatePosition(a);
               }
            }
         }
      }
   }
}

/************************************************************************/
void FindHAtomsSCMC(PDB *resA, PDB *stopA, PDB *resB, PDB *stopB, 
                    HBOND *hb)
{
   PDB *d, *a, *h1, *h2, *h3, *p;
   
   char h_name1[6], h_name2[6], h_name3[6], p_name[6];

   for(d=resA; d!=stopA; NEXT(d))
   {
      /* if residue A heavy atom is a donor */
      if(isDonor(d, hb, h_name1, h_name2, h_name3))
      {
         h1 = FindAtomInRange(resA, stopA, h_name1);
         h2 = FindAtomInRange(resA, stopA, h_name2);
         h3 = FindAtomInRange(resA, stopA, h_name3);

         /* Residue B is a mainchain acceptor */
         a = FindAtomInRange(resB, stopB, "O   ");
         p = FindAtomInRange(resB, stopB, "C   ");
               
         /* if the hydrogen bond is valid */
         if(blValidHBond(h1, d, a, p))
         {
            /* store position of partner acceptor atom */          
#ifdef NOISY
            fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                    d->resnam, d->resnum, d->atnam, 
                    a->resnam, a->resnum,  a->atnam, 
                    ((h1==NULL?"???":h1->atnam)), p->atnam); 
#endif
            StorePartnertoDonatePosition(a);
         }
         else if(h2!=NULL)
         {
            if(blValidHBond(h2, d, a, p))
            {                     
#ifdef NOISY
               fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                       d->resnam, d->resnum, d->atnam, 
                       a->resnam, a->resnum, a->atnam, 
                       h2->atnam, p->atnam);
#endif
               StorePartnertoDonatePosition(a);
            }
         }
         else if(h3!=NULL)
         {
            if(blValidHBond(h3, d, a, p))
            {
#ifdef NOISY                     
               fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                       d->resnam, d->resnum, d->atnam, 
                       a->resnam, a->resnum, a->atnam, 
                       h3->atnam, p->atnam);
#endif
               StorePartnertoDonatePosition(a);
            }
         }
      }
   }

   /* if residue A is an acceptor */
   for(a=resA; a!=stopA; NEXT(a))
   {
      if(isAcceptor(a, hb, p_name))    
      {
         p = FindAtomInRange(resA, stopA, p_name);

         /* Residue B is a mainchain donor */
         h1 = FindAtomInRange(resB, stopB, "H   ");
         
         /* If the hydrogen bond is valid */

         if(blValidHBond(h1, d, a, p))
         {
            /* store position of partners hydrogen donating heavy atom */
#ifdef NOISY
            fprintf(stderr,"%3s %5d %4s : %3s %5d %4s %4s %4s\n", 
                    d->resnam, d->resnum, d->atnam, 
                    a->resnam, a->resnum, a->atnam, 
                    ((h1==NULL?"???":h1->atnam)), p->atnam);
#endif
            StorePartnertoAcceptPosition(d);
         }
      }
   }
}

