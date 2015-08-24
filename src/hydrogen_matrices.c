/**************hydrogen_matrices.c*******************

Creates four matrices ....

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


***************************************************/


/* Includes */

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

/********************************************/

/* Defines */

/********************************************/

/* Global variables */

/* matrices that store all acceptor and donor atoms */
static int gAccept[MAXSIZE][MAXSIZE][MAXSIZE];
static int gDonate[MAXSIZE][MAXSIZE][MAXSIZE];

/* matrix storing partner atoms to hydrogen accepting atoms */
static int gPartnertoAccept[MAXSIZE][MAXSIZE][MAXSIZE];
/* matrix storing partner atoms to hydrogen donating heavy atoms */
static int gPartnertoDonate[MAXSIZE][MAXSIZE][MAXSIZE];



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

/******************************************/

/* Prototypes */

int main (int argc, char *argv[]);
HBOND *InitializeHbondTypes(void);
BOOL ParseCmdLine(int argc, char **argv, char *inputfile, char *outputfile);
void Usage(void);
NAMES *InitializeDomainList(FILE *fp);
char *FindStructureLocation(NAMES *names);
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


/******************************************/

int main (int argc, char *argv[])
{
   FILE *in = stdin, *out = stdout;
   HBOND *hb;
   char inputfile[160], outputfile[160];
   NAMES *names;
   
   inputfile[0] = outputfile[0] = '\0';
   
   if(ParseCmdLine(argc, argv, inputfile, outputfile))
   {
      if(OpenStdFiles(inputfile, outputfile, &in, &out))
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
               printf("ERROR: can't read protein domains linked list\n");
               return(1);
               
            }
         }
         else
         {
            printf("ERROR: can't read hydrogen atoms linked list\n");
            return(1);
         }
      }
      else
      {
         Usage();
      }
      
      
      
   }
   
   else
   {
      Usage();
   }
   
   return(0);
   
}
   
               
/******************************************/

BOOL CalcAndStoreHBondData(HBOND *hb, NAMES *names, FILE *out)
{
   NAMES *n;
   HBOND *h;
   int natoms, natoms2, nHatoms, l;
   PDB *pdb, *pdb2, *start, *next, *nextres, *stop;
   FILE *fp1,*fp2;
   char *location;
   BOOL noenv; 
 
   if((fp1 = OpenFile(PGPFILE, "DATADIR", "r", &noenv)) == NULL)
   {
      fprintf(stderr, "ERROR: Can't open pgp file\n");
      return(FALSE);
   }
   
   for(h = hb; h!=NULL; NEXT(h))
   {
      if(h->select)
      {
         ClearArrays();
         
         for(n=names; n !=NULL; NEXT(n))
         {
            location = FindStructureLocation(n);            
            printf("file %s\n", location);
            
            /* open protein domain file */
            if((fp2 = fopen(location, "r")) !=NULL)
            { 
               /* create linked list of pdb file */  
               if((pdb = ReadPDB(fp2, &natoms)) !=NULL)
               {   
                  /* strip any hydrogens present in protein domain file */
                  if((pdb2 = StripHPDB(pdb, &natoms2)) !=NULL)
                  {
                     FREELIST(pdb, PDB);
                     pdb = pdb2;
                  }
                  else
                  {
                     printf("WARNING: StripHPDB failed\n");
                   
                  }
                    
                  if((nHatoms = HAddPDB(fp1, pdb)) !=0)                    
                  {
                     for(start=pdb; start!=NULL; start=next)
                     {
                        next = FindNextResidue(start);
                        
                        if((!strncmp(start->resnam, h->residue, 3)))
                        {
                           /* return true if backbone atoms cannot be found .. stops program
                              progressing to next stage */
                           if((OrientatePDB(pdb, start, next) != FALSE))
                           {
                              
                              StoreHBondingPosition(start, next, hb);
                              
                              for(nextres=pdb; nextres!=NULL; nextres=stop)
                              {                                                                 
                                 stop = FindNextResidue(nextres);
                                 
                                 if((nextres !=start))
                                 {
                                    /*printf("%s %s\n", start->resnam, nextres->resnam);*/
                                    
                                    
                                    if((IsHBonded(start, nextres, HBOND_SS)) !=0)
                                    {
                                       FindHAtoms(start, next, nextres, stop, hb);
                                       
                                    }
                                 }
                                 
                                 
                              }
                           }
                           else
                           {
                              printf("backbone atoms can't be found for %s\n", location);
                           }
                           
                        }
                        
                     }
                     
                  }
                  
               }
               
               else
               {
                  printf("ERROR: can't create PDB linked list\n");
               }
               if(pdb != NULL) FREELIST(pdb, PDB); 
            }
            if(fp2 !=NULL) fclose(fp2);
            
            /* remove protein domain files created by getchain */
            
            if(l = strlen(location) == 55)
            {
               if(remove(location))
               {
                  fprintf(stderr, "WARNING: unable to remove temporary protein domain  file\n");
                  
               } 
            }
            
            
            
         }
         
         PrintMatrix(h, out);
         
      }
      
      
   }
   
   fclose(fp1);
   return(TRUE);
   
}

/**********************************************************/

void FindHAtoms(PDB *resA, PDB *stopA, PDB *resB, PDB *stopB, HBOND *hb)
{
   PDB *d, *a, *h1, *h2, *h3, *p;
   
   char h_name1[6], h_name2[6], h_name3[6], p_name[6];

   /* if residue A heavy atom is a  donor */
       
   for(d=resA; d!=stopA; NEXT(d))
   {
      if(isDonor(d, hb, h_name1, h_name2, h_name3))
      {
         /* and  residue B is a acceptor */

         for(a=resB; a!=stopB; NEXT(a))
         {
            if(isAcceptor(a, hb, p_name))
            {
               h1 = FindAtomInRange(resA, stopA, h_name1);
               h2 = FindAtomInRange(resA, stopA, h_name2);
               h3 = FindAtomInRange(resA, stopA, h_name3);
              
               p = FindAtomInRange(resB, stopB, p_name);
               
               /* if the hydrogen bond is valid */

               if(ValidHBond(h1, d, a, p))
               {
                  /* store position of partner acceptor atom */          
                  printf("%s %d %s %s %d %s %s %s\n", d->resnam, d->resnum, d->atnam, a->resnam, a->resnum,  a->atnam, ((h1==NULL?"???":h1->atnam)), p->atnam); 
                  
                  StorePartnertoDonatePosition(a);
                                    
                  
               }
               
               
               else if(h2!=NULL)
               {
                  
                  if(ValidHBond(h2, d, a, p))
                  {                     
                     printf("%s %d  %s %s %d %s %s %s\n", d->resnam, d->resnum, d->atnam, a->resnam, a->resnum, a->atnam, h2->atnam, p->atnam);
                     StorePartnertoDonatePosition(a);
                     
                       
                  }
                  
               }
               
               
               else if(h3!=NULL)
               {
                  if(ValidHBond(h3, d, a, p))
                  {
                     
                     printf("%s %d %s %s %d %s %s %s\n", d->resnam, d->resnum, d->atnam, a->resnam, a->resnum, a->atnam, h3->atnam, p->atnam);
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

               if(ValidHBond(h1, d, a, p))
               {
                  /* store position of  partners hydrogen donating heavy atom */
                  printf("%s %d %s %s %d %s %s %s\n", d->resnam, d->resnum, d->atnam, a->resnam, a->resnum, a->atnam, ((h1==NULL?"???":h1->atnam)), p->atnam);
                  StorePartnertoAcceptPosition(d);
                  
                  
               }
               
               else if(h2!=NULL)
               {
                  
                  if(ValidHBond(h2, d, a, p))
                  {
                     printf("%s %d  %s %s %d  %s %s %s\n", d->resnam, d->resnum, d->atnam, a->resnam, a->resnum, a->atnam, h2->atnam, p->atnam);
                     StorePartnertoAcceptPosition(d);
                     
                  }
                  
               }
               
               else if(h3!=NULL)
               {
                  if(ValidHBond(h3, d, a, p))
                  {
                     printf("%s %d %s %s  %d %s %s %s\n", d->resnam, d->resnum, d->atnam, a->resnam, a->resnum, a->atnam, h3->atnam, p->atnam);
                     StorePartnertoAcceptPosition(d);
                     
                    
                    
                  }
                  
               }
               
            }
            
         }
         
      }
      
   }
   
}

/**************************************************************/

BOOL isDonor(PDB *d, HBOND *hb, char *h_name1, char *h_name2, char *h_name3)
{
   HBOND *h;
   
   for(h=hb; h!=NULL; NEXT(h))
   {
      if(h->donate && !strncmp(h->residue, d->resnam, 3) &&
         !strncmp(h->hatom, d->atnam, 3))
      {
         if(h->s_donh1 !=NULL)
         { 
            strcpy(h_name1, h->s_donh1);
         }
         
         else
         {  
            h_name1[0] = '\0';
         }
         
         if(h->s_donh2 !=NULL)
         {
            strcpy(h_name2, h->s_donh2);
         }
         
         else
         {
            h_name2[0] = '\0';
         }
         
         if(h->s_donh3 !=NULL)
         {  
            strcpy(h_name3, h->s_donh3);
         }
         
         else
         {  
            h_name3[0] = '\0';
         }
         
         return(TRUE);
         
      }
      
   }
   h_name1[0] = '\0';
   return(FALSE);
}

/***************************************************************/

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

/****************************************************************/

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

/********************************************************/

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
            /*search for hydrogen-capable /residue/atoms listed in *hydrogen residue* linked list */
            
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

   //FREELIST(h, HBOND);
   
   
}


/********************************************************/

void StorePartnertoAcceptPosition(PDB *d)
{
    
   gPartnertoAccept[(int)(d->x/DIV) + OFFSET]
      [(int)(d->y/DIV) + OFFSET]
      [(int)(d->z/DIV) + OFFSET]++;   
   
   
}

/********************************************************/

void StorePartnertoDonatePosition(PDB *a)
{

   gPartnertoDonate[(int)(a->x/DIV) + OFFSET]
      [(int)(a->y/DIV) + OFFSET]
      [(int)(a->z/DIV) + OFFSET]++;   
   
   
}

/********************************************************/

void PrintMatrix(HBOND *h, FILE *out)
{

   int x, y, z;

   VEC3F grid;
   
   int atnum = 0, resnum = 0;
   

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
               fprintf(out, "donate\t%8d\t%8d\t%8d\t%6d\n",  x, y, z,  gDonate[x][y][z]);

               /*GRID_2_COORD(x, grid.x);
               GRID_2_COORD(y, grid.y);
               GRID_2_COORD(z, grid.z);
               
               fprintf(out, "ATOM  %5d  CA  THR  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
               atnum++, resnum++,  grid.x, grid.y, grid.z, 1.00, 1.00);*/
               
            }
            
            
         }
         
      }
      
   }
   /* printing out matrix containing partner atoms to hydrogen-donating heavy atoms */     
   for(x = 0; x < MAXSIZE; x++)
   {
      for(y = 0; y < MAXSIZE; y++)
      {
         for(z = 0; z < MAXSIZE; z++)
         {
            if(gPartnertoDonate[x][y][z] > 0)
            {
               
               fprintf(out, "partnertodonate\t%8d\t%8d\t%8d\t%6d\n", x, y, z,  gPartnertoDonate[x][y][z]);

               /* GRID_2_COORD(x, grid.x);
               GRID_2_COORD(y, grid.y);
               GRID_2_COORD(z, grid.z);
              
               fprintf(out, "ATOM  %5d  CA  THR  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
               atnum++, resnum++,  grid.x, grid.y, grid.z, 1.00, 2.00);*/
               
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
               
               fprintf(out, "accept\t%8d\t%8d\t%8d\t%6d\n",   x, y, z,  gAccept[x][y][z]);

               /*GRID_2_COORD(x, grid.x);
               GRID_2_COORD(y, grid.y);
               GRID_2_COORD(z, grid.z);
               
               fprintf(out, "ATOM  %5d  CA  THR  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
               atnum++, resnum++,  grid.x, grid.y, grid.z, 1.00, 3.00);*/
               
            }
            
            
         }
         
      }
      
   }
   
   /* printing out matrix containing partner atoms to hydrogen accepting atoms */
   for(x = 0; x < MAXSIZE; x++)
   {
      for(y = 0; y < MAXSIZE; y++)
      {
         for(z = 0; z < MAXSIZE; z++)
         {
            if(gPartnertoAccept[x][y][z] > 0)
            {
               
              fprintf(out, "partnertoaccept\t%8d\t%8d\t%8d\t%6d\n", x, y, z,  gPartnertoAccept[x][y][z]);

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

/********************************************************/

void ClearArrays()
{
   memset((void*)gDonate, 0, (size_t) MAXSIZE*MAXSIZE*MAXSIZE * sizeof(int));
   memset((void*)gAccept, 0, (size_t) MAXSIZE*MAXSIZE*MAXSIZE * sizeof(int));
   memset((void*)gPartnertoDonate, 0, (size_t) MAXSIZE*MAXSIZE*MAXSIZE * sizeof(int));
   memset((void*)gPartnertoAccept, 0, (size_t) MAXSIZE*MAXSIZE*MAXSIZE* sizeof(int));
   
}

/************************************************/

/* function to create linked list of protein domain files */

NAMES *InitializeDomainList(FILE *fp)
{

    NAMES *names = NULL;
   
    NAMES *n;
   
    FILE *fn = fp;
    
    char *p;
   
    char buffer[MAXBUFF];

    char resol[6];

    float resolution;
    
    while(fgets(buffer, MAXBUFF, fn))
    {

      TERMINATE(buffer);
  
      strcpy(resol, buffer+53);
      
      resolution = atof(resol);
            
      if(resolution <=2.5)
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
    }
    
      
      return(names);
}

/**********************************************************/

char *FindStructureLocation (NAMES *names)
{
   
   char *filename=(char *)malloc((size_t)(sizeof(char)*(6+1)));
   char *result,*tbuffer,*tbuffer1,chainlabel;
   
   if(strlen(names->filename) !=6) {
      fprintf(stderr,"the name does not contain 6 char: %s\n",names->filename);
      exit(1);
   } 
   filename=strcpy(filename,names->filename);
   if (filename[5]!='0') {
      /*dompdb */
      result=multiappend("%s%s",DOMAINLOC,filename);		  
      free(filename);
      return result;
   } else {
      if (filename[4]!='0') {
         /* getchain */
         chainlabel=filename[4];
         filename[4]='\0';
         tbuffer=multiappend("getchain %c %s%s%s",\
                             chainlabel,PDBLOC,filename,PDBEXT);
         tbuffer1=multiappend("%s/%s",\
                              tmpdir(),filename);
         doIt(tbuffer,tbuffer1);
         free(tbuffer);
         free(filename);
         return tbuffer1;
      } else {
         /* pdb */
         filename[4]='\0';
         result=multiappend("%s%s%s",PDBLOC,filename,PDBEXT);		  
         free(filename);
         return result;
      }
   }	  
   fprintf(stderr,"Why are we here?\n");
   exit(1);
   
}

/******************************************/

/* function that creates a linked list of all the hydrogen-capable atoms present in each amino acid */

HBOND *InitializeHbondTypes()
{
   static HBOND h1 = {"ARG", "NE  ",1, 0, 1, "HE  ", NULL, NULL, NULL };
   static HBOND h2 = {"ARG", "NH1 ",0, 0, 1, "HH11", "HH12", NULL, NULL };
   static HBOND h3 = {"ARG", "NH2 ",0, 0, 1, "HH21", "HH22", NULL, NULL };
   static HBOND h4 = {"THR", "OG1 ",1,1, 1, NULL, NULL, NULL, "CB  " };
   static HBOND h5 = {"ASN", "ND2 ",1, 0, 1, "HD21", "HD22", NULL,  NULL };
   static HBOND h6 = {"ASN", "OD1 ",0, 1, 0, NULL, NULL, NULL, "CG  "};
   static HBOND h7 = {"ASP", "OD1 ",1, 1, 0, NULL, NULL, NULL, "CG  " };
   static HBOND h8 = {"ASP", "OD2 ",0,1, 0, NULL, NULL, NULL, "CG  " };
   static HBOND h9 = {"GLU", "OE1 ",1,1, 0, NULL, NULL, NULL, "CD  " };
   static HBOND h10 = {"GLU", "OE2 ",0,1, 0, NULL, NULL, NULL, "CD  " };
   static HBOND h11 = {"GLN", "NE2 ",1,0, 1, "HE21", "HE22", NULL, NULL};
   static HBOND h12 = {"GLN", "OE1 ",0,1, 0, NULL, NULL, NULL, "CD  "};
   static HBOND h13 = {"LYS", "NZ  ",1,0, 1, "HZ1 ","HZ2 ", "HZ3 ", NULL };
   static HBOND h14 = {"SER", "OG  ",1,1, 1, NULL, NULL, NULL, "CB  "};
   static HBOND h15 = {"TRP", "NE1 ",1,0, 1, "HE1 ", NULL, NULL, NULL };
   static HBOND h16 = {"TYR", "OH  ",1,1, 1, NULL, NULL, NULL, "CZ  " };
   
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
   h16.next = NULL;
   
   
   return(first);
   
   
}

/********************************************************/

/* function to display a usage message */

void Usage(void)
{
   fprintf(stderr, "\nHydrogen Matrices V1.0 (c) 2002, Alison Cuff, University of Reading\n\n");
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
      
/**********************************************************/

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
