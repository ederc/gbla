#ifndef __GB_ouvrir_H
#define __GB_ouvrir_H

FILE* ouvrir(char* a,char* b)
{
	if ( (strcmp(a,"-") == 0) && (strcmp(b,"r") == 0) )
		return stdin;
	else
	{
		FILE *f=fopen(a,b);
		if (!f)
		{
			fprintf(stderr,"Can't open %s\n",a);
			exit(1);
		}
		return f;
	}
}

#endif

/* vim: set ft=c: */
