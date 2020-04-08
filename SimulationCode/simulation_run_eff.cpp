//Program to calculate Wright-Fisher simulations of a single generation

//This code uses a pseudorandom method to generate the initial population in an efficient way

#include <vector>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;

#include "shared.h"

vector<double> p_select;

int main(int argc, const char **argv){
    
    run_params p;
    GetOptions(p,argc,argv);
    
    //Set up random number generator
    int seed=atoi(argv[1]);
    gsl_rng_env_setup();
    gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rgen, seed);
	
	//Use a constant population size
    int N_b=p.bottle; /* bottleneck size theg*/

    if (p.multi==1&&p.process==1) {
        ProcessData(p,N_b,rgen);
        return 0;
    }
    
	//Setup : Get a list of prime numbers bigger than N_b
	vector<int> prs;
	GetPrimesList(p,N_b,prs);

    //Get consensus from variant data
	cout << "Get consensus\n";
    vector< vector< vector< double> > > init_var;
	if (p.rstart==1) {
		GetRandomVariantData(init_var,rgen);
	} else {
		GetVariantData(init_var);
	}
	vector< vector<int> > consensus;
	GetConsensus(init_var,consensus);
	
    //Generate population - eight segments
    vector< vector<haplo> > population;
    cout << "Initialising population\n";
	SetupPopulation(N_b,population);

	//Add in initial variants
    //Note : This step incorporates the reduction of variants based upon the SureSelect method statistics.  This reduction can be removed using the --cut_variants 0 flag.
    AddInitialVariants (p,N_b,prs,consensus,init_var,population,rgen);
    
	cout << "Initial population size is " << population[0].size() << "\n";

	vector< vector< vector<double> > > afs_init;
	CalculateAlleleFrequencies (N_b,consensus,population,afs_init);
    if (p.multi==1) {
        PrintInitialFreqs (afs_init);
    }
	vector< vector<haplo> > population_init=population;
    
	//Multiple replicates from same start point
    if (p.multi==0) {
        p.reps=1;
    }
	for (int rep=0;rep<p.reps;rep++) {
		//Single neutral generation - Bernoulli sampling
		for (int i=0;i<N_b;i++) {
			int s=floor(gsl_rng_uniform(rgen)*N_b);
			for (int seg=0;seg<8;seg++) {
				population[seg][i]=population_init[seg][s];
			}
		}
		
		vector< vector< vector<double> > > afs;
		CalculateAlleleFrequencies (N_b,consensus,population,afs);
        if (p.multi==1) {
            PrintRepFreqs (rep,afs);
        } else {
            //Apply a minimum frequency to the allele frequencies.  Values below the minimum set to the minimum freuqency
            //This can be controlled with the flag --min_q q in the command line
            ProcFreqs(p,afs);
                if (rep==0) {
                    ProcFreqs(p,afs_init);
                }
            //Calculate allele frequency differences
            CalcDist(afs,afs_init);
        }
	}
	
    return 0;
}


