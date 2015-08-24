/************************************************************************/
/* This is a kludge to deal with non-standard PDB files like 1dxc where
   multiple occupancy atoms are listed in 2 groups rather than in atom
   pairs
*/
PDB *CleanMultipleOccupancies(PDB *pdb)
{

   PDB *p       = NULL, 
       *q       = NULL,
       *r       = NULL,
       *start   = NULL, 
       *nextres = NULL;
   BOOL GotPartial;
   

   /* Step through the PDB list one residue at a time                   */
   for(start=pdb; start!=NULL; start=nextres)
   {
      nextres = FindNextResidue(start);

      /* See if we have any partial occupancy atoms                     */
      GotPartial = FALSE;
      for(p=start; p!=nextres; NEXT(p))
      {
         if(p->altpos != ' ')
         {
            GotPartial = TRUE;
            break;
         }
      }
      
      if(GotPartial)
      {
         for(p=start; p!= nextres; NEXT(p))
         {
            /* Got first partial position                               */
            if(p->altpos == 'A')
            {
               r=p;
               for(q=p->next; q!=nextres; NEXT(q))
               {
                  if(!strncmp(p->atnam_raw, q->atnam_raw, 4))
                  {
                     r->next = q->next;
                     free(q);
                     q=r;
                     p->altpos = ' ';
                     break;
                  }
                  r=q;
               }
            }
         }
      }
   }
   return(pdb);
}


