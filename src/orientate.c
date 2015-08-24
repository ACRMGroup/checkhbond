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
   if(!FindNCACAtoms(start, next, &c_alpha, &c_beta, &n))
   {
      printf("Unable to find backbone atoms\n");
      return(FALSE);
   }
   
   /* move protein structure file so that c_alpha is at the origin */
   TempVec.x = -c_alpha.x;
   TempVec.y = -c_alpha.y;
   TempVec.z = -c_alpha.z;

   TranslatePDB(pdb, TempVec);
   
   /* find co-ordinates of *key* residue after the translation */
  
    if(!FindNCACAtoms(start, next, &c_alpha, &c_beta, &n))
   {
      printf("Unable to find backbone atoms\n");
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
   
   if(!FindNCACAtoms(start, next, &c_alpha, &c_beta, &n))
   {
      printf("Unable to find backbone atoms");
      return(FALSE);
   }
   	 
    RotateToXY(pdb, &c_beta);

   if(!FindNCACAtoms(start, next, &c_alpha, &c_beta, &n))
   {
      printf("Unable to find backbone atoms");
      return(FALSE);
   }

   /*printf("%s %d\n",start->resnam, start->resnum);
   printf("c_alpha atoms %f %f %f\n", c_alpha.x, c_alpha.y, c_alpha.z);
   printf("c_beta atoms %f %f %f\n", c_beta.x, c_beta.y, c_beta.z);
   printf("n atoms %f %f %f\n", n.x, n.y, n.z);*/
   
  
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

BOOL FindNCACAtoms(PDB *res1_start, PDB *stop, VEC3F *c_alpha, VEC3F *c_beta, VEC3F *n)
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




