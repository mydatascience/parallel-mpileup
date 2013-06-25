#include <stdio.h>
#include <errno.h>
#include <sys/stat.h>
#include <fcntl.h>
//#include <sys/mman.h.>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>

#define MAX_N 300
#define STRSIZE 100
#define N 'N'
#define SCAFFOLD_MARKER '>'
#define STD_ERROR 0.1

void flushbed (FILE* fout, char* chr, int start, int end) {
	if (fout == NULL) {
		printf("Bed file is closed\n");
	}
	fprintf (fout, "%s\t%d\t%d\n",chr,start,end);
}

size_t get_length (FILE* finput) {
	size_t length = 0;
	char inbuf[STRSIZE];
	while (fgets(inbuf, STRSIZE, finput)) {
		if (*inbuf == SCAFFOLD_MARKER) {
			continue;
		}
		length += (strlen(inbuf) - 1);
	}
	fseek(finput, 0, SEEK_SET);
	//printf("%u\n", length);
	return length;
}

FILE * open_bed(FILE* fbed, char * refname, int currfile, char * outname) {
	if (fbed != NULL) {
		fclose(fbed);
	}
	sprintf(outname, "%s.%i.bed", refname, currfile);
	return fopen(outname, "w");
}

void get_chrname(char * chrname, char inbuf[]) {
	char * delim = strchr(inbuf, ' ');
	int size = (delim == NULL) ? strlen(inbuf) - 1 : delim - inbuf - 1;
	strncpy (chrname, inbuf+1, size);
	chrname[size] = '\0';
}

int main (int argc, char** argv) {
	FILE* finput = NULL;
	FILE* fbed = NULL;
	char inbuf[STRSIZE];
	char chrname[STRSIZE];
	int ncount = 0;
	size_t regstart, currpos, bedstart, lastchr, last_n_end, globalpos;
	int currfile = 1;
	char* curr;
	errno = 0;
	char* outname;
	size_t ref_length = 0;
	size_t block_size = 0;
	int files_num = 4;
	char* refname;

	if (argc != 3) {
		printf("USAGE: ./group_divider reference.fasta files_number\n");
		return 0;
	}

	printf ("Opening in file %s\n", argv[1]);
	finput = fopen(argv[1],"r");
	if (finput == NULL)  {
		printf ("file open error #%d !\n", errno);
		return 1;
	}
	
	curr = basename(argv[1]);
	refname = malloc (strchr(curr, '.') - curr);
	strncpy(refname, curr, strchr(curr, '.') - curr);
	printf("%s\n", refname);

	outname = malloc (strlen(argv[1])+8);
	fbed = open_bed(fbed, refname, currfile, outname);

	if (argc == 3) {
		printf("Counting blocks size for %s threads\n", argv[2]);
		files_num = atoi(argv[2]);
	}
	ref_length = get_length(finput);
	block_size = ref_length / files_num;

	regstart = last_n_end = lastchr = currpos = globalpos = bedstart = 0;
	while (fgets(inbuf, STRSIZE, finput)) {
		inbuf [strlen(inbuf)-1] = '\0';	//Carriage return
		if ((curr = strchr(inbuf, '\r'))) {
			*curr = '\0';
		}

		if (*inbuf == SCAFFOLD_MARKER) {
			if (currpos) {
				flushbed (fbed, chrname,regstart,currpos);
				lastchr = globalpos;
			}
			get_chrname(chrname, inbuf);
			printf ("Parsing scaffold '%s'\n", chrname);
			regstart = ncount = currpos = 0;
			continue;
		}

		if (*inbuf == '\0') continue; 	// Empty string

		if (globalpos >= block_size * (currfile + STD_ERROR)) {
			currfile += 1;
			printf("Bed file #%i is started\n", currfile);
			if (lastchr > block_size * (currfile - 1 - STD_ERROR)) {
				fbed = open_bed(fbed, refname, currfile, outname);
			} else if (last_n_end > block_size * (currfile - 1 - STD_ERROR)) {
				flushbed(fbed, chrname, regstart, last_n_end);
				regstart = last_n_end;
			} else {
				flushbed(fbed, chrname, regstart, currpos - block_size * STD_ERROR);
				regstart = currpos - block_size * STD_ERROR;
				fbed = open_bed(fbed, refname, currfile, outname);
			}
		}

		curr = inbuf;
		while (*curr) {
			if (*curr == '\n') {
				printf ("Malformed data!\n");
			}
			if (*curr == N ) {
				ncount++;
			}
			else {
				if (ncount >= MAX_N) {
					if (regstart != currpos - ncount) {
						last_n_end = currpos - ncount / 2;
					}
				}
				ncount = 0;
			}
			currpos++;
			globalpos++;
			curr++;
		}
	}
	if (currpos) {
		flushbed (fbed, chrname,regstart,currpos);
    }
	return 0;
}
	
