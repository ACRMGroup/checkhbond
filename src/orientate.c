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
#include "orientate.h"

/************************************************************************/
BOOL OrientatePDB(PDB *pdb, PDB *start, PDB *next)
{
   VEC3F c_alpha,
      c_beta,
      n,
      TempVec;

   /* find co-ordinates of *key* residue */
   if(!FindNCACBAtoms(start, next, &c_alpha, &c_beta, &n))
   {
      PrintError(NULL, "Error (checkhbond): 1. Unable to find backbone atoms\n");
      return(FALSE);
   }
   
   /* move protein structure file so that c_alpha is at the origin */
   TempVec.x = -c_alpha.x;
   TempVec.y = -c_alpha.y;
   TempVec.z = -c_alpha.z;

   TranslatePDB(pdb, TempVec);
   
   /* find co-ordinates of *key* residue after the translation */
  
   if(!FindNCACBAtoms(start, next, &c_alpha, &c_beta, &n))
   {
      PrintError(NULL,"2. Unable to find backbone atoms\n");
      return(FALSE);
   }

   /* rotate n so that it is on the xz plane */
            
   RotateToXZ(pdb, &n);
   
   /* rotate n onto the x axis */

   RotateToX(pdb, &n);

   /** ACRM c_beta has moved now!!! because of the rotations in the 2 previous steps
       c-beta has to be refound!
       Add this line!
   **/
   
   if(!FindNCACBAtoms(start, next, &c_alpha, &c_beta, &n))
   {
      PrintError(NULL,"3. Unable to find backbone atoms\n");
      return(FALSE);
   }
   	 
    RotateToXY(pdb, &c_beta);

   if(!FindNCACBAtoms(start, next, &c_alpha, &c_beta, &n))
   {
      PrintError(NULL,"4. Unable to find backbone atoms\n");
      return(FALSE);
   }

#ifdef DEBUG
   printf("DEBUG: residue %s %d\n",start->resnam, start->resnum);
   printf("DEBUG: c_alpha atoms %f %f %f\n", c_alpha.x, c_alpha.y, c_alpha.z);
   printf("DEBUG: c_beta atoms %f %f %f\n", c_beta.x, c_beta.y, c_beta.z);
   printf("DEBUG: n atoms %f %f %f\n", n.x, n.y, n.z);
#endif
  
   return(TRUE); 

}

/************************************************************************/
BOOL OrientateN_PDB(PDB *pdb, PDB *prev, PDB *start, PDB *next)
{
   VEC3F c_alpha,
      c,
      n,
      TempVec;

   /* find co-ordinates of *key* residue */
   if(!FindCNCAAtoms(prev, start, next, &c, &n, &c_alpha))
   {
      PrintError(NULL,"5. Unable to find backbone atoms\n");
      return(FALSE);
   }
   
   /* move protein structure file so that N is at the origin */
   TempVec.x = -n.x;
   TempVec.y = -n.y;
   TempVec.z = -n.z;

   TranslatePDB(pdb, TempVec);
   
   /* find co-ordinates of *key* residue after the translation */
  
   if(!FindCNCAAtoms(prev, start, next, &c, &n, &c_alpha))
   {
      PrintError(NULL,"6. Unable to find backbone atoms\n");
      return(FALSE);
   }

   /* rotate c so that it is on the xz plane */
            
   RotateToXZ(pdb, &c);
   
   /* rotate c onto the x axis */

   RotateToX(pdb, &c);

   /** ACRM c_beta has moved now!!! because of the rotations in the 2 previous steps
       c-beta has to be refound!
       Add this line!
   **/
   
   if(!FindCNCAAtoms(prev, start, next, &c, &n, &c_alpha))
   {
      PrintError(NULL,"7. Unable to find backbone atoms\n");
      return(FALSE);
   }
   	 
    RotateToXY(pdb, &c_alpha);

   if(!FindCNCAAtoms(prev, start, next, &c, &n, &c_alpha))
   {
      PrintError(NULL,"8. Unable to find backbone atoms\n");
      return(FALSE);
   }

#ifdef DEBUG
   printf("DEBUG: residue %s %d\n",start->resnam, start->resnum);
   printf("DEBUG: c atoms %f %f %f\n", c.x, c.y, c.z);
   printf("DEBUG: n atoms %f %f %f\n", n.x, n.y, n.z);
   printf("DEBUG: c_alpha atoms %f %f %f\n", c_alpha.x, c_alpha.y, c_alpha.z);
#endif
  
   return(TRUE); 

}

/**************************************************************/

/* function to rotate
PDB so that the n atom is on the xz plane */

 void RotateToXZ(PDB *pdb, VEC3F *n)
{
   REAL angle;
   REAL matrix[3][3];
   VEC3F n_angle;
   
   angle = TrueAngle(n->y, n->x);
   
   CreateRotMat('z', -angle, matrix);
   
   ApplyMatrixPDB(pdb, matrix);
   MatMult3_33(*n, matrix, &n_angle);
   *n = n_angle;
   
}

/*************************************************************/

/* function to rotate PDB so that the n atom is on the x axis */

void RotateToX(PDB *pdb, VEC3F *n)
{
   REAL angle;
   REAL matrix[3][3];
   VEC3F n_angle;
   
   angle = TrueAngle(n->z, n->x);
   
   CreateRotMat('y', angle, matrix);
   
   ApplyMatrixPDB(pdb, matrix);
   MatMult3_33(*n, matrix, &n_angle);
   *n = n_angle;
}

/***********************************************************/

/* function to rotate PDB so that the c_beta atom is on the xy plane */

void RotateToXY(PDB *pdb, VEC3F *c_beta)
{
   
   REAL angle;
   REAL matrix[3][3];
   VEC3F c_beta_angle;
   
   angle = TrueAngle(c_beta->z, c_beta->y);
   
   CreateRotMat('x', -angle, matrix);
   
   ApplyMatrixPDB(pdb, matrix);
   MatMult3_33(*c_beta, matrix, &c_beta_angle);
   *c_beta = c_beta_angle;
}

/*********************************************************/

/* function that returns angle when given the opposite and adjacent
 values. Uses atan2 */

REAL TrueAngle(REAL opp, REAL adj)
{
   REAL angle;
   
   angle = atan2(opp,adj);
   
   return(angle);
}
 

/*********************************************************/
/* function to find and return co-ordinates of c_alpha, c_beta 
   and n atoms of *key* residue */
BOOL FindNCACBAtoms(PDB *res1_start, PDB *stop, VEC3F *c_alpha, VEC3F *c_beta, VEC3F *n)
{
   PDB *q;

   BOOL c_alpha_found = FALSE, 
        c_beta_found = FALSE,
        n_found = FALSE;
   
   for(q = res1_start; q != stop; NEXT(q))
   {
      if(!strncmp(q->atnam, "CA", 2))
      {
         c_alpha->x = q->x;
         c_alpha->y = q->y;
         c_alpha->z = q->z;
         c_alpha_found = TRUE;
         
         
      }
      
      if(!strncmp(q->atnam, "CB", 2))
      {
         c_beta->x = q->x;
         c_beta->y = q->y;
         c_beta->z = q->z;
         c_beta_found = TRUE;
        
      }
      
      if(!strncmp(q->atnam, "N   ", 4))
      {
         n->x = q->x;
         n->y = q->y;
         n->z = q->z;  
         n_found = TRUE;
      }
      
   }
   if(!c_alpha_found || !c_beta_found || !n_found)
   {
      return(FALSE);
   }
   return(TRUE);
   
}




/*********************************************************/

/* function to find and return co-ordinates of c, n and c_alpha
   atoms of *key* residue */

BOOL FindCNCAAtoms(PDB *res0_start, PDB *res1_start, PDB *stop, VEC3F *c, VEC3F *n, VEC3F *c_alpha)
{
   PDB *q;

   BOOL c_alpha_found = FALSE, 
        c_found = FALSE,
        n_found = FALSE;

   for(q = res0_start; q != res1_start; NEXT(q))
   {
      if(!strncmp(q->atnam, "C   ", 4))
      {
         c->x = q->x;
         c->y = q->y;
         c->z = q->z;
         c_found = TRUE;
         break;
      }
   }
   
   for(q = res1_start; q != stop; NEXT(q))
   {
      if(!strncmp(q->atnam, "N   ", 4))
      {
         n->x = q->x;
         n->y = q->y;
         n->z = q->z;  
         n_found = TRUE;
      }

      if(!strncmp(q->atnam, "CA  ", 4))
      {
         c_alpha->x = q->x;
         c_alpha->y = q->y;
         c_alpha->z = q->z;
         c_alpha_found = TRUE;
      }
   }
   if(!c_alpha_found || !c_found || !n_found)
   {
      return(FALSE);
   }
   return(TRUE);
   
}


/************************************************************************/
BOOL OrientateCO_PDB(PDB *pdb, PDB *start, PDB *next)
{
   VEC3F c_alpha,
      c,
      o,
      TempVec;

   /* find co-ordinates of *key* residue */
   if(!FindCACOAtoms(start, next, &c_alpha, &c, &o))
   {
      PrintError(NULL,"9. Unable to find backbone atoms\n");
      return(FALSE);
   }
   
   /* move protein structure file so that c is at the origin */
   TempVec.x = -c.x;
   TempVec.y = -c.y;
   TempVec.z = -c.z;

   TranslatePDB(pdb, TempVec);
   
   /* find co-ordinates of *key* residue after the translation */
  
   if(!FindCACOAtoms(start, next, &c_alpha, &c, &o))
   {
      PrintError(NULL,"10. Unable to find backbone atoms\n");
      return(FALSE);
   }

   /* rotate ca so that it is on the xz plane */
            
   RotateToXZ(pdb, &c_alpha);
   
   /* rotate ca onto the x axis */

   RotateToX(pdb, &c_alpha);
   
   if(!FindCACOAtoms(start, next, &c_alpha, &c, &o))
   {
      PrintError(NULL,"11. Unable to find backbone atoms\n");
      return(FALSE);
   }
   	 
    RotateToXY(pdb, &o);

   if(!FindCACOAtoms(start, next, &c_alpha, &c, &o))
   {
      PrintError(NULL,"12. Unable to find backbone atoms\n");
      return(FALSE);
   }

#ifdef DEBUG
   printf("DEBUG: residue %s %d\n",start->resnam, start->resnum);
   printf("DEBUG: c_alpha atoms %f %f %f\n", c_alpha.x, c_alpha.y, c_alpha.z);
   printf("DEBUG: c atoms %f %f %f\n", c.x, c.y, c.z);
   printf("DEBUG: o atoms %f %f %f\n", o.x, o.y, o.z);
#endif
  
   return(TRUE); 

}



/*********************************************************/
/* function to find and return co-ordinates of c_alpha, c
   and o atoms of *key* residue */
BOOL FindCACOAtoms(PDB *res1_start, PDB *stop, VEC3F *c_alpha, VEC3F *c, VEC3F *o)
{
   PDB *q;

   BOOL c_alpha_found = FALSE, 
        c_found = FALSE,
        o_found = FALSE;
   
   for(q = res1_start; q != stop; NEXT(q))
   {
      if(!strncmp(q->atnam, "CA", 2))
      {
         c_alpha->x = q->x;
         c_alpha->y = q->y;
         c_alpha->z = q->z;
         c_alpha_found = TRUE;
      }
      
      if(!strncmp(q->atnam, "C ", 2))
      {
         c->x = q->x;
         c->y = q->y;
         c->z = q->z;
         c_found = TRUE;
        
      }
      
      if(!strncmp(q->atnam, "O   ", 4) ||
         !strncmp(q->atnam, "O1  ", 4) ||
         !strncmp(q->atnam, "OT1 ", 4))
      {
         o->x = q->x;
         o->y = q->y;
         o->z = q->z;  
         o_found = TRUE;
      }
      
   }
   if(!c_alpha_found || !c_found || !o_found)
   {
      return(FALSE);
   }
   return(TRUE);
   
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
void PrintError(FILE *out, char *error)
{
   fprintf(stderr,"Error (checkhbond): %s", error);
   if(out != NULL)
   {
      fprintf(out,"Error (checkhbond): %s", error);
   }
}

