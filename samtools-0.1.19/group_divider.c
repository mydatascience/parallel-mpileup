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
		fprintf (stderr,"Bed file is closed\n");
	}
	fprintf (fout, "%s\t%d\t%d\n",chr,start,end);
}


void echobed (char* chr, int start, int end) {
	fprintf (stderr,"[BED] %s\t%d\t%d\n",chr,start,end);
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
	//fprintf (stderr,"%u\n", length);
	return length;
}

//Generating array of BEDS to return;
typedef struct bed_list bed_list;
struct bed_list {
	char* bed;
	bed_list* next;
};
bed_list* beds_array = NULL;
bed_list* bed_tail;
void beds_array_add (char* newbed) {
	bed_list* bed = malloc (sizeof(beds_array));
	if (beds_array==NULL) {
		beds_array = bed;
		bed_tail = beds_array;
	}
	else {
		bed_tail->next = bed;
		bed_tail = bed_tail->next;
	}
	fprintf (stderr,"Adding BED %s...\n", newbed);
	bed_tail->bed = newbed;
}

int bed_list_count () {
	bed_list* bed = beds_array;
	int i=0;
	if (bed==NULL) {fprintf (stderr,"BUG!\n"); return 0;}
	while (bed) {
		i++;
		bed = bed->next;
	}
	return i;
}

char** bed_array_return () {
	char** result = malloc (sizeof(char*) * 10/*bed_list_count()*/);
	char** cur = result;
	bed_list* cbed = beds_array;
	if (cbed==NULL) fprintf (stderr,"BUG!\n");
	while (cbed) {
		*cur++ = cbed->bed;
		fprintf (stderr,"BED is %s!\n", cbed->bed);
		cbed=cbed->next;
	}
	return result;
}

FILE * open_bed(FILE* fbed, char * refname, int currfile, char * outname) {
	if (fbed != NULL) {
		fclose(fbed);
	}
	sprintf (outname, "%s.%i.bed", refname, currfile);
	beds_array_add (outname);
	return fopen(outname, "w");
}

void start_bed (int threadnum) {
	fprintf (stderr,"Opening thread #%d\n",threadnum);
}

void get_chrname(char * chrname, char inbuf[]) {
	char * delim = strchr(inbuf, ' ');
	int size = (delim == NULL) ? strlen(inbuf) - 1 : delim - inbuf - 1;
	strncpy (chrname, inbuf+1, size);
	chrname[size] = '\0';
}

int group_divider (char* filename, int threads, char** beds) {
	fprintf (stderr,"Group Divider!\n");
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

//	if (argc != 3) {
//		fprintf (stderr,"USAGE: ./group_divider reference.fasta files_number\n");
//		return 0;
//	}

	fprintf (stderr,"Opening in file %s\n", filename);
	finput = fopen(filename,"r");
	if (finput == NULL)  {
		fprintf (stderr,"file open error #%d !\n", errno);
		return 1;
	}


	curr = basename(filename);
	refname = malloc (strchr(curr, '.') - curr);
	strncpy(refname, curr, strchr(curr, '.') - curr);
	fprintf (stderr,"%s\n", refname);

	outname = malloc (strlen(filename)+8);
//	fprintf (stderr,"Outname = %s\n",outname);
	fbed = open_bed(fbed, refname, currfile, outname);

//	start_bed (currfile);


//	if (argc == 3) {
		fprintf (stderr,"Counting blocks size for %d threads\n", threads);
		files_num = threads;
//	}
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
//				echobed (chrname,regstart,currpos);
				lastchr = globalpos;
			}
			get_chrname(chrname, inbuf);
			fprintf (stderr,"Parsing scaffold '%s'\n", chrname);
			regstart = ncount = currpos = 0;
			continue;
		}

		if (*inbuf == '\0') continue; 	// Empty string

		if (globalpos >= block_size * (currfile + STD_ERROR)) {
			currfile += 1;
			fprintf (stderr,"Bed file #%i is started\n", currfile);
			if (lastchr > block_size * (currfile - 1 - STD_ERROR)) {
				fbed = open_bed(fbed, refname, currfile, outname);
//				start_bed (currfile);
			} else if (last_n_end > block_size * (currfile - 1 - STD_ERROR)) {
				flushbed(fbed, chrname, regstart, last_n_end);
//				echobed (chrname,regstart,last_n_end);
				regstart = last_n_end;
			} else {
				flushbed(fbed, chrname, regstart, currpos - block_size * STD_ERROR);
//				echobed (chrname,regstart,currpos - block_size * STD_ERROR);
				regstart = currpos - block_size * STD_ERROR;
				fbed = open_bed(fbed, refname, currfile, outname);
//				start_bed (currfile);
			}
		}

		curr = inbuf;
		while (*curr) {
			if (*curr == '\n') {
				fprintf (stderr,"Malformed data!\n");
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
//		echobed (chrname,regstart,currpos);
    }
	beds =  bed_array_return ();
	return 0;
}
	
