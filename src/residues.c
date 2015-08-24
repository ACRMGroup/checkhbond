/*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include <math.h>
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/matrix.h"
#include "bioplib/hbond.h"
#include "bioplib/fit.h"
#include "residues.h"

/************************************************************************/

void WriteResidue(FILE *fp, PDB  *start, PDB *stop)
{
   PDB  *p;
   char  PrevChain[8];
   
   strcpy(PrevChain,start->chain);

   for(p = start ; p != stop ; NEXT(p))
   {
      if(strncmp(PrevChain,p->chain,1))
      {
         /* Chain change, insert TER card                               */
         fprintf(fp,"TER   \n");
         strcpy(PrevChain,p->chain);
      }
      WriteResidueRecord(fp,p);
   }
   fprintf(fp,"TER   \n");
}

/*************************************/

void WriteResidueRecord(FILE *fp, PDB  *pdb)
{
   fprintf(fp,"%-6s%5d %-5s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
           pdb->record_type,
           pdb->atnum,
           pdb->atnam_raw,
           pdb->resnam,
           pdb->chain,
           pdb->resnum,
           pdb->insert,
           pdb->x,
           pdb->y,
           pdb->z,
           pdb->occ,
           pdb->bval);
}

/*************************************/

void ApplyMatrixResidue(PDB  *start, PDB *stop, REAL matrix[3][3])
{
   PDB   *p;
   VEC3F incoords,
         outcoords;

   for(p=start; p!=stop; NEXT(p))
   {
      if(p->x != 9999.0 && p->y != 9999.0 && p->z != 9999.0)
      {
         incoords.x = p->x;
         incoords.y = p->y;
         incoords.z = p->z;
         MatMult3_33(incoords,matrix,&outcoords);
         p->x = outcoords.x;
         p->y = outcoords.y;
         p->z = outcoords.z;
      }
   }
}

/*************************************/

void TranslateResidue(PDB *start, PDB *stop, VEC3F tvect)
{
   PDB *p;

   for(p=start; p!=stop; NEXT(p))
   {
      if(p->x < (REAL)9999.0 && p->y < (REAL)9999.0 && p->z < (REAL)9999.0)
      {
         p->x += tvect.x;
         p->y += tvect.y;
         p->z += tvect.z;
      }
   }
}
  
/*************************************/

PDB *SelectAtomsResidue(PDB *start, PDB *stop, int nsel, char **sel, int *natom)
{
   PDB *p, *q,  *pdbout = NULL;
   int i;
   
   
   *natom = 0;
   
   /* step though the residue atom at a time */
   for(p= start; p!=stop; NEXT(p))
   {  
      for(i=0; i<nsel; i++)
      {
         
         if(!strncmp(p->atnam, sel[i], 4))
         {
            
            /* Allocate a new entry */
            if(pdbout==NULL)
            {
               INIT(pdbout, PDB);
               q = pdbout;
            }
            else
            {
               ALLOCNEXT(q, PDB);
            }
            
            /* If failed, free anything allocated and return */
            if(q==NULL)
            {
               if(pdbout != NULL) FREELIST(pdbout,PDB);
               *natom = 0;
               return(NULL);
            }
            
            /* Increment atom count*/
            (*natom)++;
            
            /* Copy the record to the output list (sets ->res1_next to NULL) */
            CopyPDB(q, p);
            
            break;
         }
         
      }
   }
      
   /* Return pointer to res1_start of output list */
   return(pdbout);
   
}









