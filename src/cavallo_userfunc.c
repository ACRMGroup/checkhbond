#include "cavallo_userfunc.h"
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/fcntl.h>
#include <unistd.h>

#define TMPDIR "/tmp/"

#define MAXBUFF 1024

#define MAKEMESSAGE \
	char *p;  \
	va_list ap;  \
	int size=100,n;  \
	if ((p=(char *)malloc( (size_t) sizeof(char)*size )) == NULL) goto WEEEE;  \
	while(1) {  \
		va_start(ap,fmt);  \
		n=vsnprintf(p,size,fmt,ap);  \
		va_end(ap);  \
		if( n > -1 && n < size)	goto WEEEE;  \
		if( n > -1 )  \
			size=n+1;   \
		else  \
			size*=2;  \
		if( (p=(char*)realloc(p,size))	== NULL) goto WEEEE;  \
	} \
WEEEE:  \
	if(p==NULL) { fprintf(stderr,"unable to allocate memory in internalMessage"); exit(4) ; } 



char * stemp (void) {
	char *buffer;		
	FILE * fp;
	int fd;
	buffer=multiappend("%s/deleteme.XXXXXX",TMPDIR);
	fd=mkstemp(buffer);
	if (fd==-1) {
		perror("error in giveMeTemporaryFile (mkstemp)");
		exit(1);
	}
	fp=fdopen(fd,"w");
	if(!fp) {
		perror("error in giveMeTemporaryFile (fdopen)");
		perror(buffer);
		exit(1);
	}
	fclose(fp);
	close(fd);
	return (buffer);
}
		
FILE *sfopen(char * filename, char *mode) {
	FILE *fp=fopen(filename,mode);
	char *buffer=(char *)malloc((size_t) ( sizeof(char) * MAXBUFF));
	sprintf(buffer,"error in sfopen %s",filename);
	if(!fp) {
		perror(buffer);
		exit(1);
	}	
	
	return fp;
}


char * multiappend(char *fmt, ...) {
	MAKEMESSAGE
	return p;	
}

void doIt (char * command,char *fileoutput) {
	char *buffer=(char *)malloc((size_t) ( sizeof(char) * MAXBUFF));
	FILE *fp;
	FILE *fpout;
	size_t bcount=1;
	
	fp=popen(command,"r");
	if(!fileoutput){
		fpout=fopen("/tmp/cavallo.out","w");
	} else {
		fpout=fopen(fileoutput,"w");
	}
	while(bcount!=0){
		bcount=fread(buffer,1,MAXBUFF,fp);
		fwrite(buffer,1,bcount,fpout);
	}
	fclose(fpout);
	fclose(fp);
        free(buffer);
        
}

#define PDBLOC "/data/pdb/pdb"
#define DOMAINLOC "/acrm/data/dompdb/"
#define PDBEXT ".ent"

char * giveMeStructureLocation(char * name) {
	char *filename=(char *)malloc((size_t)(sizeof(char)*(6+1)));
	char *result,*tbuffer,*tbuffer1,chainlabel;
	if(strlen(name) !=6) {
	fprintf(stderr,"the name does not conatin 6 char: %s\n",name);
		exit(1);
	} 
	filename=strcpy(filename,name);
	if (filename[5]!='0') {
		/*dompdb */
		result=multiappend("%s%s",DOMAINLOC,filename);		  
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
			return result;
		}
	}	  
		   
}
char * tmpdir (void) {
	return TMPDIR; 
}

#undef MAKEMESSAGE 
#undef TMPDIR
#undef MAXBUFF 

