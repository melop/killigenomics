// Updated on 05/22/15

/*
 
   Program Make_InFile.cpp to make
   individual input files of nucleotide
   read quartets at every position of the
   reference genome

*/

#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
using namespace std;

int main(int argc, char *argv[])
{	
	char ref_file[100];
	int num_ind;
	int ig;
	char indiv[100];
	char in_file[100];
	char out_file[100];
	char ss1_out[10] = {"In_GFE_"};
	char ss2_out[10] = {".txt"};
	char ref_line[1000];
	int ref_num_lines;
	char ind_line[1000];
	int ind_num_read;
	char ind_scaf[100];
	string s_ind_scaf;
	int ind_site;
	int nA, nC, nG, nT;
	char header_scaf[10], header_site[10], header_ref[10];
	int ref_num_read;
	char ref_scaf[100];
	string s_ref_scaf;
	int ref_site;
	char ref_nuc[10];
	int match;
	
	FILE *instream_1;
	FILE *instream_2;
        FILE *outstream;
	
	// Read the name of the reference-genome file
	sscanf(argv[1], "%s", ref_file);

	// Read the number of individuals in the analysis
	sscanf(argv[2], "%d", &num_ind);
	printf("%d individuals\n", num_ind);

	if (argc != 3+2*num_ind) {	// A check to make sure the user has passed the right number of arguments.
		printf("Please check the input file.\n"); 
		exit(1);
	}	
	
	// Make the necessary file for each individual
	for (ig = 1; ig<=num_ind; ++ig) {
		// Read the input file name
		sscanf(argv[ig+2], "%s", in_file);

		// Read the individual ID
		sscanf(argv[ig+2+num_ind], "%s", indiv);

		// Output file name
		sprintf(out_file, "%s%s%s", ss1_out, indiv, ss2_out);

		// Open the output file
		outstream = fopen(out_file, "w");
		if (outstream == NULL ) {	// Exit on failure
			fprintf(stderr, "Cannot open %s for writing.\n", out_file); 
			exit(1);
		}

		// Open the reference input file for reading.
		instream_1 = fopen(ref_file, "r");
		if (instream_1 == NULL ) {       // Exit on failure
                        fprintf(stderr, "Cannot open %s for reading.\n", ref_file);
                        exit(1);
                }
		// Open the pro file of an individual
		instream_2 = fopen(in_file, "r");
		if (instream_2 == NULL ) {       // Exit on failure
                        fprintf(stderr, "Cannot open %s for reading.\n", in_file);
                        exit(1);
                }

		ref_num_lines = 0;
		while ( !feof(instream_2) ) {
			if(fgets(ind_line, 999, instream_2) == NULL)
                        	continue;
			ind_num_read = sscanf(ind_line, "%s %d %d %d %d %d", ind_scaf, &ind_site, &nA, &nC, &nG, &nT);
			s_ind_scaf = ind_scaf;
			match = 0;
			if (ind_num_read != 6) {
                        	fprintf(stderr, "wrong line format in line: %s ", ind_line);
                        	continue;
                	}
			while (match == 0) {
				if ( feof(instream_1) && match == 0) {
					fprintf(stderr, "Could not find site %d on %s in the reference file for %s.\n", ind_site, ind_scaf, in_file);
					fprintf(stderr, "Do not use %s.\n", out_file);
					exit(1);
				}
				if(fgets(ref_line, 999, instream_1) == NULL)
					continue;	
				if (ref_num_lines == 0){	// header line?
					ref_num_read = sscanf(ref_line, "%s %s %s", header_scaf, header_site, header_ref);
					fprintf(outstream, "%s\n", indiv);
					printf("%s\n", indiv);
					if (ref_num_read != 3) {
                                		fprintf(stderr, "wrong line format in line: %s ", ref_line);
                                		continue;
                        		}
					ref_num_lines = ref_num_lines + 1;
				} else {
					ref_num_read = sscanf(ref_line, "%s %d %s", ref_scaf, &ref_site, ref_nuc);
					s_ref_scaf = ref_scaf;
					if (ref_num_read != 3) {
						fprintf(stderr, "wrong line format in line: %s ", ref_line);
						continue;
					}
					ref_num_lines = ref_num_lines + 1;
					if (s_ref_scaf == s_ind_scaf && ref_site == ind_site) {
						fprintf(outstream, "%d/%d/%d/%d\n", nA, nC, nG, nT);
						// printf("%d/%d/%d/%d\n", nA, nC, nG, nT);
						match = 1;
					} else {
						fprintf(outstream, "0/0/0/0\n");
						// printf("0/0/0/0\n");
					}
				}
			}
		}
		while ( !feof(instream_1) ) {
			if(fgets(ref_line, 999, instream_1) == NULL)
        			continue;
			ref_num_read = sscanf(ref_line, "%s %d %s", ref_scaf, &ref_site, ref_nuc);
			ref_num_lines = ref_num_lines + 1;
        		if (ref_num_read != 3) {
        			fprintf(stderr, "wrong line format in line: %s ", ref_line);
                		continue;
        		}
			fprintf(outstream, "0/0/0/0\n");
        		// printf("0/0/0/0\n");
		}

		// Close the reference input file and individual file
		fclose(instream_1);
		fclose(instream_2);
	}  
	
	// Close the output file
	fclose(outstream);	
	return 0;
}
