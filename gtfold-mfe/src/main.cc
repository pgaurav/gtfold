/*
GTfold: compute minimum free energy of RNA secondary structure
Copyright (C) 2008  David A. Bader
http://www.cc.gatech.edu/~bader

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Authored by Amrita Mathuriya August 2007 - January 2009.*/
/* Modified by Sonny Hernandez May 2007 - Aug 2007. All comments added marked by "SH: "*/
/* Modified by Sainath Mallidi August 2009 - "*/

#include <errno.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "loader.h"
#include "algorithms.h"
#include "algorithms-partition.h"
#include "traceback.h"
#include "subopt_traceback.h"
#include "main.h"
#include "main-c.h"

using namespace std;

/* GLOBAL VARIABLES */
enum BOOL NOISOLATE;
enum BOOL BPP; // calculating base pair probabilities
enum BOOL USERDATA;
enum BOOL PARAMS;
enum BOOL LIMIT_DISTANCE;
enum BOOL VERBOSE;

int num_threads;
int LENGTH;
int contact_dist;

unsigned char *RNA1;
unsigned char *RNA; /* Contains RNA string in terms of 0, 1, 2, 3 for A, C, G and U respectively*/
int *structure; /* An array to contain the optimal structure */
int *V; /* int V[LENGTH][LENGTH]; */
int *W;
int **VBI; /* VBI(i,j) will contain the energy of optimal internal loop closed with (i,j) base pair */
int **VM; /* VM(i, j) will contain the energy of optimla multiloop closed with (i,j) base pair */
int **WM; /* This array is introduced to help multiloop calculations. WM(i,j) contains the optimal energy of string segment from si to sj if this forms part of a multiloop */
int *indx; /* This array is used to index V array. Here V array is mapped from 2D to 1D and indx array is used to get the mapping back.*/
int *constraints;

double **QB;  // QB[i][j] is the sum over all possible loops closed by (i,j),
// including the summed contributions of their subloops
double **Q;   // Q[i][j] in addition to the above quantity QB[i][j], Q[i][j]
// also includes all configurations with (i,j) not paired
double **QM;  // QM[i][j] is the sum of configuration energies from i to j,
// assuming that i,j are contained in a multiloop
double **P;   // P[i][j] The probability that nucleotides i and j form a basepair

/**
 * Print the hep message and quit.
 */
void help() {
    fprintf(stderr,
            "Usage: gtfold [OPTION]... FILE\n\n");

    fprintf(stderr,
            "  FILE is an RNA sequence file.  Single line or FASTA formats are accepted.\n\n");

    fprintf(stderr, "OPTIONS\n");
    fprintf(stderr,
            "   -c, --constraints FILE\n                        Load constraints from FILE.  See Constraint syntax below\n");
    fprintf(stderr,
            "   -d, --limitCD num    Set a maximum base pair contact distance to num. If no\n                        limit is given, base pairs can be over any distance\n");
    fprintf(stderr,
            "   -n, --noisolate      Prevent isolated base pairs from forming\n");
    fprintf(stderr,
            "   -o, --output FILE    Output to FILE (default output is to a .ct extension)\n");
    fprintf(stderr,
            "   -v, --verbose        Run in versbose mode\n");
    
	fprintf(stderr,
            "   -t, --threads num    Limit number of threads used\n");
    fprintf(stderr,
            "   -h, --help           Output help (this message) and exit\n");

    fprintf(stderr, "\nBETA OPTIONS\n");
    fprintf(stderr,
            "   --bpp                Calculate base pair probabilities\n");
    fprintf(stderr,
            "   --subopt range       Calculate suboptimal structures within 'range' kcal/mol\n");
    fprintf(stderr, "                        of the mfe\n");


    fprintf(stderr,
            "\nConstraint syntax:\n\tF i j k  # force (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) pairs\n\tP i j k  # prohibit (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) pairs\n\tP i 0 k  # make bases from i to i+k-1 single stranded bases.\n");
    exit(-1);
}

/* Function for calculating time */
double get_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;
}

/*Function for printing the sequence*/
void printSequence(int len) {
    int i = 1;
    for (i = 1; i <= len; i++) {
        if (RNA1[i] == 0)
            printf("A");
        else if (RNA1[i] == 1)
            printf("C");
        else if (RNA1[i] == 2)
            printf("G");
        else if (RNA1 [i] == 3)
            printf("T");
        else
            printf("N");
    }
    printf("\n");
}

/*Function for printing the input constraints*/
void printConstraints(int len) {
    int i = 1;
    for (i = 1; i <= len; i++) {
        if (constraints[i] > 0 && constraints[i] > i)
            printf("(");
        else if (constraints[i] > 0 && constraints[i] < i)
            printf(")");
        else if (constraints[i] < 0)
            printf("x");
        else
            printf(".");
    }
    printf("\n");
}

/*Function for printing the predicted structure*/
void printStructure(int len) {
    int i = 1;
    for (i = 1; i <= len; i++) {
        if (structure[i] > 0 && structure[i] > i)
            printf("(");
        else if (structure[i] > 0 && structure[i] < i)
            printf(")");
        else
            printf(".");
    }
    printf("\n");
}

/* Initialize global variables. */
void init_variables(int len) {

    int i;

    LENGTH = len + 1;

    RNA = (unsigned char *) malloc(LENGTH * sizeof(unsigned char));
    if (RNA == NULL) {
        perror("Cannot allocate variable 'RNA'");
        exit(-1);
    }
    RNA1 = (unsigned char *) malloc(LENGTH * sizeof(unsigned char));
    if (RNA1 == NULL) {
        perror("Cannot allocate variable 'RNA'");
        exit(-1);
    }
    structure = (int *) malloc(LENGTH * sizeof(int));
    if (structure == NULL) {
        perror("Cannot allocate variable 'structure'");
        exit(-1);
    }

    V = (int *) malloc(((LENGTH - 1) * (LENGTH) / 2 + 1) * sizeof(int));
    if (V == NULL) {
        perror("Cannot allocate variable 'V'");
        exit(-1);
    }

    W = (int *) malloc(LENGTH * sizeof(int));
    if (W == NULL) {
        perror("Cannot allocate variable 'W'");
        exit(-1);
    }

    VBI = (int **) malloc(LENGTH * sizeof(int *));
    if (VBI == NULL) {
        perror("Cannot allocate variable 'VBI'");
        exit(-1);
    }
    for (i = 0; i < LENGTH; i++) {
        VBI[i] = (int *) malloc(LENGTH * sizeof(int));
        if (VBI[i] == NULL) {
            perror("Cannot allocate variable 'VBI[i]'");
            exit(-1);
        }
    }

    VM = (int **) malloc(LENGTH * sizeof(int *));
    if (VM == NULL) {
        perror("Cannot allocate variable 'VM'");
        exit(-1);
    }
    for (i = 0; i < LENGTH; i++) {
        VM[i] = (int *) malloc(LENGTH * sizeof(int));
        if (VM[i] == NULL) {
            perror("Cannot allocate variable 'VM[i]'");
            exit(-1);
        }
    }

    WM = (int **) malloc(LENGTH * sizeof(int *));
    if (WM == NULL) {
        perror("Cannot allocate variable 'WM'");
        exit(-1);
    }
    for (i = 0; i < LENGTH; i++) {
        WM[i] = (int *) malloc(LENGTH * sizeof(int));
        if (WM[i] == NULL) {
            perror("Cannot allocate variable 'WM[i]'");
            exit(-1);
        }
    }

    indx = (int *) malloc(LENGTH * sizeof(int));
    if (indx == NULL) {
        perror("Cannot allocate variable 'indx'");
        exit(-1);
    }

    constraints = (int*) malloc((len + 1) * sizeof(int));
    if (constraints == NULL) {
        perror("Cannot allocate variable 'constraints'");
        exit(-1);
    }
}

/* deallocate global variables */
void free_variables() {
    int i;

    free(indx);
    for (i = 0; i < LENGTH; i++)
        free(WM[i]);
    free(WM);
    for (i = 0; i < LENGTH; i++)
        free(VM[i]);
    free(VM);
    for (i = 0; i < LENGTH; i++)
        free(VBI[i]);
    free(VBI);
    free(W);
    free(V);
    free(constraints);
    free(structure);
    free(RNA);
    free(RNA1);
}

void init_partition_function_variables(int length) {
    QB = mallocTwoD(length+1, length+1);
    if(QB == NULL) {
        fprintf(stderr,"Failed to allocate QB\n");
        exit(-1);
    }

    Q = mallocTwoD(length+1, length+1);
    if(Q == NULL) {
        fprintf(stderr,"Failed to allocate Q\n");
        exit(-1);
    }

    QM = mallocTwoD(length+1, length+1);
    if(QM == NULL) {
        fprintf(stderr,"Failed to allocate QM\n");
        exit(-1);
    }

    P = mallocTwoD(length+1, length+1);
    if(P == NULL) {
        fprintf(stderr,"Failed to allocate P\n");
        exit(-1);
    }
}

void free_partition_function_variables(int length) {
    freeTwoD(QB, length+1, length+1);
    freeTwoD(Q, length+1, length+1);
    freeTwoD(QM, length+1, length+1);
    freeTwoD(P, length+1, length+1);
}


/* main function - This calls
 *  1) Read command line arguments.
 *  2) populate() from loader.cc to read the thermodynamic parameters defined in the files given in data directory.
 *  3) Initialize variables
 *  4) Calls calculate function defined in algorithms.c for filling up energy tables.
 *  5) Calls trace function defined in trace.c file for tracing the optimal secondary structure
 *  6) Then it generates .ct file from the 1D array structure.
 *  */
int main(int argc, char** argv) {
    int i;
    ifstream cf;
    int length;
    string s, seq;
    int energy;
    double t1;
    NOISOLATE = FALSE;
    BPP = FALSE;
    VERBOSE = FALSE;
	int delta = -1;

    fprintf(stdout,
            "GTfold: A Scalable Multicore Code for RNA Secondary Structure Prediction\n");
    fprintf(
            stdout,
            "(c) 2007-2010  D.A. Bader, S. Mallidi, A. Mathuriya, C.E. Heitsch, S.C. Harvey\n");
    fprintf(stdout, "Georgia Institute of Technology\n\n");

    /* Reading command line arguments */
    if (argc < 2)
        help();

    int fileIndex = 0, consIndex = 0, dataIndex = 0, lcdIndex = 0;
    int rangeIndex = 0, thIndex = 0, outputFileIndex = 0;

    i = 1;
    while (i < argc) {
        if (argv[i][0] == '-') {
            if        (strcmp(argv[i], "--noisolate") == 0 ||
                    strcmp(argv[i], "-n") == 0) {
                NOISOLATE = TRUE;
            } else if (strcmp(argv[i], "--output") == 0 ||
                    strcmp(argv[i], "-o") == 0) {
                if (i < argc)
                    outputFileIndex = ++i;
                else
                    help();
            } else if (strcmp(argv[i], "--verbose") == 0 ||
                    strcmp(argv[i], "-v") == 0) {
                VERBOSE = TRUE;
            } else if (strcmp(argv[i], "--help") == 0 ||
                    strcmp(argv[i], "-h") == 0) {
                help();
            } else if (strcmp(argv[i], "--constraints") == 0 ||
                    strcmp(argv[i], "-c") == 0) {
                if (i < argc)
                    consIndex = ++i;
                else
                    help();
            } else if (strcmp(argv[i], "--threads") == 0 ||
                    strcmp(argv[i], "-t") == 0) {
                if (i < argc)
                    thIndex = ++i;
                else
                    help();
            } else if (strcmp(argv[i], "--datadir") == 0) {
                USERDATA = TRUE;
                if (i < argc)
                    dataIndex = ++i;
                else
                    help();
            } else if (strcmp(argv[i], "--limitCD") == 0 ||
                    strcmp(argv[i], "-d") == 0) {
                if (i < argc)
                    lcdIndex = ++i;
                else
                    help();
            } else if (strcmp(argv[i], "--bpp") == 0) {
                BPP = TRUE;
            } else if (strcmp(argv[i], "--subopt") == 0) {
				if (i<argc)
                    rangeIndex = ++i;
				else
					help();
            } else {
                fprintf(stderr, "Unrecognized option: `%s'\n", argv[i]);
                fprintf(stderr, "Try: `%s --help' for more information.\n", argv[0]);
                exit(1);
            }

        } else {
            fileIndex = i;
        }
        i++;
    }

    if (fileIndex == 0)
        help();

    fprintf(stdout, "File: %s\n", argv[fileIndex]);
    cf.open(argv[fileIndex], ios::in);
    if (cf == NULL) {
        fprintf(stdout, "File open failed.\n\n");
        exit(-1);
    }

    seq = "";
    s = "";

    //Handle FASTA input
    char ss[10000];
    cf.getline(ss, 10000);
    if (ss[0] != '>') {
        char *fline;
        fline = strtok(ss, " ");
        while (fline != NULL) {
            seq.append(fline);
            fline = strtok(NULL, " ");
        }
    }

    while (!cf.eof()) {
        cf >> s;
        seq.append(s);
        s = "";
    }
    s = seq;

    length = s.length();

    init_variables(length);

    fprintf(stdout, "Sequence length: %5d\n\n", length);

    cf.close();

    int **fbp = NULL, **pbp = NULL;
    int numfConstraints = 0, numpConstraints = 0;

    if (consIndex != 0) {
        int r = initialize_constraints(&fbp, &pbp, numpConstraints, numfConstraints, argv[consIndex]);
        if (r == ERR_OPEN_FILE) {
            free_variables();
            exit(-1);
        }
        else if (r == NO_CONS_FOUND) {
            fprintf(stdout, "constraints file empty! Quitting...\n");
            exit(1);
        }
    }

    if (thIndex > 0) {
        num_threads = atoi(argv[thIndex]);
    }

    if (handle_IUPAC_code(s, length)  == FAILURE) {
        free_variables();
        exit(0);
    }

    if(USERDATA==TRUE)
        populate(argv[dataIndex],true);
    else
        populate("Turner99",false); /* Defined in loader.cc file to read in the thermodynamic parameter values from the tables in the ../data directory. */

    initTables(length); /* Initialize global variables */

    ///// Run Configuration Output /////

    fprintf(stdout, "Run Configuration:\n");
    bool standardRun = TRUE;

    if (NOISOLATE == TRUE) {
        fprintf(stdout, "- preventing isolated base pairs\n");
        standardRun = FALSE;
    }

    if (consIndex != 0) {
        fprintf(stdout, "- using constraint file: %s\n", argv[consIndex]);
        standardRun = FALSE;
    }
	
	contact_dist = length+1;
    if (lcdIndex != 0) {
        contact_dist = atoi(argv[lcdIndex]);
        fprintf(stdout, "- maximum contact distance: %d\n", contact_dist);
        standardRun = FALSE;
    }

    if (BPP == TRUE) {
        fprintf(stdout, "+ calculating base pair probabilities\n");
        standardRun = FALSE;
    }

    if (rangeIndex != 0) {	
        delta = atoi(argv[rangeIndex]);
        fprintf(stdout, "+ calculating suboptimal structures within %d kcal/mol of MFE\n", delta);
        standardRun = FALSE;
	}

    if(standardRun)
        fprintf(stdout, "- standard\n");

    ///// end run config output /////


    fprintf(stdout,"\nComputing minimum free energy structure ... \n\n");
    fflush(stdout);

    t1 = get_seconds();
    energy = calculate(length, fbp, pbp, numfConstraints, numpConstraints);
    t1 = get_seconds() - t1;

    fprintf(stdout,"Minimum Free Energy: %9.2f\n", energy/100.00);
    fprintf(stdout,"MFE runtime (seconds): %9.6f\n\n", t1);

	if (rangeIndex != 0)
	{
		fprintf(stdout,"Computing suboptimal structures in range %d ...\n", delta);
		t1 = get_seconds();
		// Traces the optimal structure
		ss_map_t suboptdata =  subopt_traceback(length, delta); 
		t1 = get_seconds() - t1;
		fprintf(stdout,"Traceback runtime (seconds): %9.6f\n\n", t1);

		size_t pos = 0;
		std::string suboptfile;
		suboptfile += argv[fileIndex];
		
		if ((pos=suboptfile.find_last_of('/')) > 0) {
			suboptfile = suboptfile.substr(pos+1);
		}

		// if an extension exists, replace it with struct
		if(suboptfile.find(".") != string::npos)
			suboptfile.erase(suboptfile.rfind("."));
		suboptfile += ".struct";
		ofstream subopt_fp;
		
		fprintf(stdout, "Writing suboptimal structures: %s\n", suboptfile.c_str());
		subopt_fp.open (suboptfile.c_str());

		ss_map_t::iterator it;
		for (it = suboptdata.begin(); it != suboptdata.end(); ++it)
			subopt_fp << it->first << ' ' << it->second << std::endl;

		subopt_fp.close();

		exit(0);	
	}

    // only fill the partition function structures if they are needed for BPP
    if(BPP) {
        fprintf(stdout,"Filling Partition Function structure... \n\n");
        fflush(stdout);

        // malloc the arrays
        init_partition_function_variables(length);

        // fill the arrays
        fill_partition_fn_arrays(length, QB, Q, QM);

        fprintf(stdout,"Q[1][n]: %f\n\n", Q[1][length]);
    }


    t1 = get_seconds();
    trace(length); /* Traces the optimal structure*/
    t1 = get_seconds() - t1;

    std::stringstream ss1, ss2;

    ss1 << length;
    ss2 << energy/100.0;

    i = 0;

    // the output filename
    string outputfile;

    // use given output file
    if( outputFileIndex != 0)
        outputfile += argv[outputFileIndex];

    // or build off the input file
    else {
		size_t pos;
        outputfile += argv[fileIndex];
		
		if ((pos=outputfile.find_last_of('/')) > 0) {
			outputfile = outputfile.substr(pos+1);
		}

        // if an extension exists, replace it with ct
        if(outputfile.find(".") != string::npos)
            outputfile.erase(outputfile.rfind("."));

        outputfile += ".ct";
    }

    ofstream outfile;
    outfile.open ( outputfile.c_str() );

    fprintf(stdout, "Writing secondary structure: %s\n", outputfile.c_str());

    outfile << length << "\t  dG = " << energy/100.0;
    i = 1;
    while ( i <= length ) {
        outfile << endl << i << "\t" << s[i-1] << "\t" << i-1 << "\t" << (i+1)%(length+1) << "\t" << structure[i] << "\t" << i;
        i++;
    }
    outfile << endl;
    outfile.close();

    fprintf(stdout,"Traceback runtime (seconds): %9.6f\n\n", t1);

	printSequence(length);
	if (consIndex != 0) { 
		printConstraints(length);
	}
	printStructure(length);

    if(BPP) {
        fillBasePairProbabilities(length, structure, Q, QB, QM, P);

        printBasePairProbabilities(length, structure, P);

        free_partition_function_variables(length);
    }

    free_variables();

    return 0;
}


int initialize_constraints(int*** fbp, int ***pbp, 
						   int& numpConstraints, int& numfConstraints, 
						   const char* constr_file, enum BOOL verbose) {
    ifstream cfcons;

    fprintf(stdout, "Running with constraints\n");
    fprintf(stdout, "Opening constraint file: %s\n", constr_file);

    cfcons.open(constr_file, ios::in);
    if (cfcons != NULL)
        fprintf(stdout, "Constraint file opened.\n");
    else {
        fprintf(stderr, "Error opening constraint file\n\n");
        cfcons.close();
        return ERR_OPEN_FILE; //exit(-1);
    }

    char cons[100];

    while (!cfcons.eof()) {
        cfcons.getline(cons, 100);
        if (cons[0] == 'F' || cons[0] == 'f')
            numfConstraints++;
        if (cons[0] == 'P' || cons[0] == 'p')
            numpConstraints++;
    }
    cfcons.close();

    fprintf(stdout, "Number of Constraints given: %d\n\n", numfConstraints
            + numpConstraints);
    if (numfConstraints + numpConstraints != 0)
        fprintf(stdout, "Reading Constraints.\n");
    else {
        fprintf(stderr, "No Constraints found.\n\n");
        return NO_CONS_FOUND;
    }

    *fbp = (int**) malloc(numfConstraints * sizeof(int*));
    *pbp = (int**) malloc(numpConstraints * sizeof(int*));

    int fit = 0, pit = 0, it = 0;

    for (it = 0; it < numfConstraints; it++) {
        (*fbp)[it] = (int*) malloc(2* sizeof (int));
    }
    for(it=0; it<numpConstraints; it++) {
        (*pbp)[it] = (int*)malloc(2*sizeof(int));
    }
    cfcons.open(constr_file, ios::in);

    while(!cfcons.eof()) {
        cfcons.getline(cons,100);
        char *p=strtok(cons, " ");
        p = strtok(NULL, " ");
        if(cons[0]=='F' || cons[0]=='f') {
            int fit1=0;
            while(p!=NULL) {
                (*fbp)[fit][fit1++] = atoi(p);
                p = strtok(NULL, " ");
            }
            fit++;
        }
        if( cons[0]=='P' || cons[0]=='p') {
            int pit1=0;
            while(p!=NULL) {
                (*pbp)[pit][pit1++] = atoi(p);
                p = strtok(NULL, " ");
            }
            pit++;
        }
    }

	if (VERBOSE == TRUE)
	{
		fprintf(stdout, "Forced base pairs: ");
		for(it=0; it<numfConstraints; it++) {
			for(int k=1;k<= (*fbp)[it][2];k++)
				fprintf(stdout, "(%d,%d) ", (*fbp)[it][0]+k-1, (*fbp)[it][1]-k+1);
		}
		fprintf(stdout, "\nProhibited base pairs: ");
		for(it=0; it<numpConstraints; it++) {
			for(int k=1;k<= (*pbp)[it][2];k++)
				fprintf(stdout, "(%d,%d) ", (*pbp)[it][0]+k-1, (*pbp)[it][1]-k+1);
		}
		fprintf(stdout, "\n\n");
	}

    return SUCCESS;
}

int handle_IUPAC_code(const std::string& s, const int length) {
	int* stack_unidentified_base;
	int stack_count=0;
	bool unspecd=0;
	stack_unidentified_base=new int[length];

	// Conversion of the sequence to numerical values. 
	for(int i = 1; i <= length; i++) {
		RNA[i] = getBase(s.substr(i-1,1));
		RNA1[i] = getBase1(s.substr(i-1,1));
		if (RNA[i]=='X') {
			fprintf(stderr,"ERROR: Base unrecognized\n");
			return FAILURE; 
		}
		else if(RNA[i]!='X' && RNA1[i]=='N'){
			unspecd=1;
			stack_unidentified_base[stack_count]=i;
			stack_count++;
		}

	}
	if (unspecd) {
		printf("IUPAC codes have been detected at positions: ");

		for(int i = 0; i < stack_count; i++) {
			printf("%d",stack_unidentified_base[i]);
			if (i < stack_count -1) printf(", ");
		}
		printf("\n");

		printf("You may wish to resubmit the sequence with fully specified positions.\n");
		printf("Aborting the current computation.\n\n");

		return FAILURE;
	}
	return SUCCESS;
}


//void force_noncanonical_basepair(const char* ncb, int len) {
//    if (ncb == 0 || ncb[0] == '\0') return;
//
//    printf("Permitted noncanonical base pairs : \n");	
//
//   std::string ncb1(ncb);
//
//   for (unsigned int i =0 ; i < ncb1.size(); ++i) {
//       ncb1[i] = toupper(ncb1[i]);
//   }
//
//    std::vector<std::string> tokens;
//    tokenize(ncb1, tokens, ",");	
//
//    for (unsigned int i = 0; i < tokens.size(); ++i) {
//       trim_spaces(tokens[i]);
//       if (tokens.size() != 3 && tokens[i][1] != '-')
//           // ignore
//           continue;
//
//      char b1 = getBase(tokens[i].substr(0,1));
//      char b2 = getBase(tokens[i].substr(2,1));
//
//      int r1=0;
//      r1 = update_chPair(b1, b2);			
//     if (r1 == 1)
//         printf("(%c,%c) ",  tokens[i][0], tokens[i][2]) ;
// }
// printf("\n\n");
//}  
