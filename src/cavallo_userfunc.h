#ifndef __cavallo_userfuncs_xxx12
#define __cavallo_userfuncs_xxx12

#include <stdio.h>
#include <stdarg.h>

char * stemp (void) ;
char * tmpdir(void);
FILE * sfopen(char * filename, char *mode) ;
char * multiappend(char *fmt , ... ) ;

char * giveMeStructureLocation(char * name) ;
void doIt (char * command, char *fileoutput) ;

#endif
