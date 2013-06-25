#include <stdio.h>
#include <errno.h>
#include <sys/stat.h>
#include <fcntl.h>
//#include <sys/mman.h.>
#include <stdlib.h>
#include <string.h>


#define MAX_N 3
#define STRSIZE 62
#define N 'N'
#define SCAFFOLD_MARKER '>'


void flushbed (char* chr, int start, int end) {
	printf ("%s\t%d\t%d\n",chr,start,end);
}

int main (int argc, char** argv) {
	FILE* finput;
	FILE* fbed;
	char inbuf[STRSIZE];
	char chrname[STRSIZE];
	int regstart, regend, ncount,globalpos;
	char* carr;
	errno = 0;
//	int fd;
//	char* infile;
	char* outname;
	char* outfile;
	printf ("Opening in file %s\n", argv[1]);
	finput = fopen(argv[1],"r");
	if (finput == NULL)  {
		printf ("file open error #%d !\n", errno);
		return 1;
	}
	outname = malloc (strlen(argv[1])+5);
	memcpy (outname, argv[1], strlen(argv[1]));
	memcpy (&outname[strlen(argv[1])],".bed",5);
	printf ("Opening out file %s\n",outname);
	fbed = fopen(outname,"r");
	if (fbed !=NULL) {
		printf ("Warning! BED file exists!\n");
		fclose (fbed);
	}
	fbed = fopen(outname,"w");
	if (fbed == NULL)  {
		printf ("file open error #%d !\n", errno);
		return 1;
	}
	regstart=regend=globalpos=0;
	while (fgets(inbuf, STRSIZE, finput)) {
//		printf ("!");
		inbuf [strlen(inbuf)-1] = '\0';	//Carriage return
//		printf ("String: %s\n",inbuf);
		if (*inbuf == SCAFFOLD_MARKER) {
			if (globalpos) {
				flushbed (chrname,regstart,globalpos);
			}
			memcpy (chrname, inbuf+1, strlen(inbuf));
			printf ("Parsing scaffold '%s'\n", chrname);
			regstart=0;
			ncount=0;
			globalpos=0;
			continue;
		}
		if (*inbuf == '\0') continue; 	// Empty string
		carr = inbuf;
		while (*carr) {
			if (*carr == '\n') {
				printf ("Malformed data!");
			}
//			printf ("%c ",*carr);
			if (*carr== N ) {
				ncount++;
//				printf ("%d ",ncount);
			}
			else {
				if (ncount>=MAX_N) {
//					printf ("MAXN!");
					flushbed (chrname,regstart,globalpos-ncount);
					regstart=globalpos;
				}
				ncount=0;
//				printf ("! ");
			}
//		printf ("-> %s <-", inbuf);
			globalpos++;
			carr++;
		}
	}
	if (globalpos) {
		flushbed (chrname,regstart,globalpos);
        }
	return 0;
}
	
