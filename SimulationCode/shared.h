
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>


using namespace std;

struct run_params {
	int bottle;
	int rstart;  //Flag to start simulation using the diversity measurements from a random point in the data.  If set to zero data from the first time point is used.
    int multi;  //Number of repeat simulations to run from the same start point
	double min_q;  //Imposed minimum allele frequency
    int reps; //Number of simulations to run from the same start point
    int process;  //Flag to process previous data: Requires 200 outputs from --multi 1 flagged code
    int cut_variants; //Flag to remove variation from the initial population using statistics from the analysis of SureSelect replicate samples
    int add_ld; //Flag to impose linkage disequilibrium within the initial population
};

struct var {  //Unit for sparse storage of mutations
	int loc;
	int nuc;
	int ns;
	double fit;
};

struct haplo {
	vector<var> seq;
	double fitness;
};

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>

//io_new
void GetVariantData (vector< vector < vector<double> > >& init_var);
void GetRandomVariantData (vector< vector < vector<double> > >& init_var, gsl_rng *rgen);
void PrintInitialFreqs (vector< vector< vector<double> > >& afs_init);
void PrintRepFreqs (int rep, vector< vector< vector<double> > >& afs);
void ReadInitFreqs (vector< vector< vector<double> > >& afs);
void ReadRepFreqs (int rep, vector< vector< vector<double> > >& afs);


//utilities
void GetOptions (run_params& p, int argc, const char **argv);
void ProcessData (run_params p, int& N_b, gsl_rng *rgen);
void MakePermutation (int n, int k, vector<int> aa, vector<int>& p, gsl_rng *rgen);
void GetPrimesList (run_params p, int N_b, vector<int>& prs);
void GetConsensus (vector< vector< vector< double> > >& init_var, vector< vector<int> >& consensus);
void SetupPopulation(int N_max, vector< vector<haplo> >& population);
void AddInitialVariants (run_params p, int N_b, vector<int>& prs, vector< vector<int> >& consensus, vector< vector< vector< double> > >& init_var, vector< vector<haplo> >& population, gsl_rng *rgen);
void CutInitialFrequencies(double& q, gsl_rng *rgen);
void GetPrimePermutation (int& chk, int N_b, int k, vector<int>& perm, vector<int>& prs, gsl_rng *rgen);
void CalculateAlleleFrequencies (double size, vector< vector<int> > consensus, vector< vector<haplo> >& population, vector< vector< vector<double> > >& afs);
void ProcFreqs (run_params p, vector< vector< vector<double> > >& afs);
void CalcDist (vector< vector< vector<double> > >& afs, vector< vector< vector<double> > >& afs_init);




