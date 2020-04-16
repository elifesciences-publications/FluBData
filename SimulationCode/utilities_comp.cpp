#include "shared.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.bottle=100;
	p.rstart=0;
    	p.multi=0;
	p.min_q=0.001;
	p.reps=10;
	p.process=0;
	p.cut_variants=1;
	p.nperms=50000;
	int x=2;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--bottle")==0) {
			x++;
			p.bottle=atoi(argv[x]);
		} else if (p_switch.compare("--rstart")==0) {
			x++;
			p.rstart=atoi(argv[x]);
        } else if (p_switch.compare("--nperms")==0) {
            x++;
            p.nperms=atoi(argv[x]);
        } else if (p_switch.compare("--multi")==0) {
            x++;
            p.multi=atoi(argv[x]);
        } else if (p_switch.compare("--process")==0) {
            x++;
            p.process=atoi(argv[x]);
        } else if (p_switch.compare("--cut_variants")==0) {
            x++;
            p.cut_variants=atoi(argv[x]);
		} else if (p_switch.compare("--minq")==0) {
			x++;
			p.min_q=atof(argv[x]);
        } else if (p_switch.compare("--reps")==0) {
            x++;
            p.reps=atoi(argv[x]);
		} else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void ProcessData (run_params p, int& N_b, gsl_rng *rgen) {
    N_b=N_b/1000000;
    //Read in initial frequencies
    vector< vector< vector<double> > > afs_init;
    ReadInitFreqs(afs_init);
    for (int r=0;r<p.reps;r++) {
        //Read in final frequencies
        vector< vector< vector< vector< double> > > > afs_sets;
        vector< vector< vector< double> > > afs_mean;
        //Use permutation to get files
        vector<int> aa;
        for (int i=0;i<200;i++) {
            aa.push_back(i);
        }
        vector<int> pp;
        MakePermutation(200,N_b,aa,pp,rgen);
        for (int i=0;i<N_b;i++) {
            vector< vector< vector< double> > > afs;
            int rep=pp[i];
            ReadRepFreqs(rep,afs);
            afs_sets.push_back(afs);
        }
        for (int i=0;i<N_b;i++) {
            vector< vector< vector< double> > > afs;
            int rep=pp[i];
            ReadRepFreqs(rep,afs);
            afs_sets.push_back(afs);
        }

        //Find mean allele frequency spectrum
        afs_mean=afs_sets[0];
        for (int i=1;i<N_b;i++) {
            for (int j=0;j<afs_sets[i].size();j++) {
                for (int k=0;k<afs_sets[i][j].size();k++) {
                    for (int l=0;l<afs_sets[i][j][k].size();l++) {
                        afs_mean[j][k][l]=afs_mean[j][k][l]+afs_sets[i][j][k][l];
                    }
                }
            }
        }
        for (int i=0;i<afs_mean.size();i++) {
            for (int j=0;j<afs_mean[i].size();j++) {
                for (int k=0;k<afs_mean[i][j].size();k++) {
                    afs_mean[i][j][k]=afs_mean[i][j][k]/N_b;
                }
            }
        }
        ProcFreqs(p,afs_init);
        ProcFreqs(p,afs_mean);
        CalcDist(afs_mean,afs_init);
    }
}

void MakePermutation (int n, int k, vector<int> aa, vector<int>& p, gsl_rng *rgen) {
    p.clear();
    for (int i=0;i<k;i++) {
        p.push_back(-1);
        int r=floor(gsl_rng_uniform(rgen)*aa.size());
        p[i]=aa[r];
        aa.erase(aa.begin()+r);
    }
}

void GetPrimesList (run_params p, int N_b, vector<int>& prs) {
    int lp=0;
    int val=N_b+1;
    int sqv=floor(sqrt(val+0.1))+1;
    while (lp<5000) {
        int isp=1;
        for (int i=2;i<=sqv;i++) {
            if (val%i==0) {
                isp=0;
                break;
            }
        }
        if (isp==1) {
            lp++;
            prs.push_back(val);
        }
        val++;
    }
}

void GetPrimeFactors (int phi, vector<int>& factors) {
    int found=0;
    while (phi%2==0) {
        phi=phi/2;
        if (found==0) {
            factors.push_back(2);
            found=1;
        }
    }
    int sqv=floor(sqrt(phi+0.1))+1;
    for (int i=3;i<=sqv;i=i+2) {
        found=0;
        while (phi%i==0) {
            phi=phi/i;
            if (found==0) {
                factors.push_back(i);
                found=1;
            }
        }
    }
    if (phi>1) {
        factors.push_back(phi);
    }
}

int GetPrimitiveRoot (int pr, int phi, vector<int> factors) {
    for (int i=2;i<=phi;i++) {
        int found=1;
        for (int j=0;j<factors.size();j++) {
            int r=FindPower(i,phi/factors[j],pr);
            if (r==1) {
                found=0;
                break;
            }
        }
        if (found==1) {
            return i;
        }
    }
    return -1;
}

void GetPrimitiveRoots (int pr, int phi, vector<int> factors, vector<int>& roots) {
    for (int i=2;i<=phi;i++) {
        int found=1;
        for (int j=0;j<factors.size();j++) {
            int r=FindPower(i,phi/factors[j],pr);
            if (r==1) {
                found=0;
                break;
            }
        }
        if (found==1) {
            roots.push_back(i);
        }
    }
}

int FindPower (int x, int y, int p) {
    int res=1;
    x=x%p;
    while (y>0) {
        if (y & 1)
            res=(res*x)%p;
        
        y = y >> 1;
        x=(x*x)%p;
    }
    return res;
}

void GetConsensus (vector< vector< vector< double> > >& init_var, vector< vector<int> >& consensus) {
    for (int i=0;i<8;i++) {
        vector<int> con;
        for (int j=0;j<init_var[i].size();j++) {
            double max=init_var[i][j][0];
            int c=0;
            if (init_var[i][j][1]>max) {
                max=init_var[i][j][1];
                c=1;
            }
            if (init_var[i][j][2]>max) {
                max=init_var[i][j][2];
                c=2;
            }
            if (init_var[i][j][3]>max) {
                max=init_var[i][j][3];
                c=3;
            }
            con.push_back(c);
        }
        consensus.push_back(con);
    }
}

void SetupPopulation(int N_max, vector< vector<haplo> >& population) {
    vector<haplo> pop;
    for (int i=0;i<N_max;i++) {
        haplo h;
        h.fitness=1.0;
        pop.push_back(h);
    }
    for (int i=0;i<8;i++) {
        population.push_back(pop);
    }
}

void AddInitialVariants (run_params p, int N_b, vector<int>& prs, vector< vector<int> >& all_roots, vector< vector<int> >& consensus, vector< vector< vector< double> > >& init_var, vector< vector<haplo> >& population, gsl_rng *rgen) {
    cout << "Add initial variants\n";
    cout << init_var.size() << "\n";
    cout << population.size() << "\n";
    for (int i=0;i<8;i++) {
        cout << "Segment " << i+1 << " " << init_var[i].size() << " "  << consensus[i].size() << "\n";
        for (int j=0;j<init_var[i].size();j++) {
            vector<int> Nm_vals;
            vector<int> k_vals;
            double tot_Nm=0;
            for (int k=0;k<4;k++) {
                if (k!=consensus[i][j]) {
                    double q=init_var[i][j][k];
                    //Filtering - account for sequencing noise
                    if (p.cut_variants==1) {
                        CutInitialFrequencies(q,rgen);
                    }
                    int N_m=gsl_ran_binomial(rgen,q,N_b);
                    if (N_m>0) {
                        Nm_vals.push_back(N_m);
                        k_vals.push_back(k);
                        tot_Nm=tot_Nm+N_m;
                    }
                }
            }
            if (tot_Nm>0) {
                vector<int> p;
                GetPrimePermutation (N_b,tot_Nm,p,prs,all_roots,rgen);
                int st=0;
                for (int k=0;k<Nm_vals.size();k++) {
                    var v;
                    v.loc=j;
                    v.nuc=k_vals[k];
                    for (int l=st;l<st+Nm_vals[k];l++) {
                        population[i][p[l]].seq.push_back(v);
                    }
                    st=st+Nm_vals[k];
                }
            }
        }
    }

}

void CutInitialFrequencies(double& q, gsl_rng *rgen) {
	if (q>=0.02) {
		double rn=gsl_rng_uniform(rgen);
		if (rn>0.999175) {
			q=0;
		}
	}
	if (q<0.02&&q>=0.015) {
		double rn=gsl_rng_uniform(rgen);
		if (rn>0.990431) {
			q=0;
		}
	}
	if (q<0.015&&q>=0.01) {
		double rn=gsl_rng_uniform(rgen);
		if (rn>0.98364) {
			q=0;
		}
	}
	if (q<0.01&&q>=0.008) {
		double rn=gsl_rng_uniform(rgen);
		if (rn>0.977723) {
			q=0;
		}
	}
	if (q<0.008&&q>=0.006) {
		double rn=gsl_rng_uniform(rgen);
		if (rn>0.931298) {
			q=0;
		}
	}
	if (q<0.006&&q>=0.004) {
		double rn=gsl_rng_uniform(rgen);
		if (rn>0.891073) {
			q=0;
		}
	}
	if (q<0.004&&q>=0.003) {
		double rn=gsl_rng_uniform(rgen);
		if (rn>0.837442) {
			q=0;
		}
	}
	if (q<0.003&&q>=0.002) {
		double rn=gsl_rng_uniform(rgen);
		if (rn>0.804382) {
			q=0;
		}
	}
	if (q<0.002&&q>=0.001) {
		double rn=gsl_rng_uniform(rgen);
		if (rn>0.302934) {
			q=0;
		}
	}
	if (q<0.001&&q>=0.0005) {
		double rn=gsl_rng_uniform(rgen);
		if (rn>0.200889) {
			q=0;
		}
	}
	if (q<0.0005) {
		q=0;
	}
}

void GetPrimePermutation (int N_b, int k, vector<int>& perm, vector<int>& prs, vector< vector<int> >& all_roots, gsl_rng *rgen) {
    //Random prime number
    int rp=floor(gsl_rng_uniform(rgen)*all_roots.size());
    int pp=prs[rp]; //Prime number
    //Find random primitive root
    int rs=floor(gsl_rng_uniform(rgen)*all_roots[rp].size());
    long long r=all_roots[rp][rs]; //Primitive root
    long long r_orig=r;
   //cout << "Orig " << r_orig << "\n";
    vector<int> ps;
    int index=0;
    if (r<=N_b) {
        perm.push_back(r-1);
        index++;
    }
    int rstore=-1;
    while (index<k) {
        r=(r*r_orig)%pp;
        if (r<0) {
            r=r+pp;
        }
        if (r<=N_b) {  //Some numbers generated in prime permutation will be too large.  Remove these
            perm.push_back(r-1);
            index++;
        }
        rstore=r;
    }
}

void CalculateAlleleFrequencies (double size, vector< vector<int> > consensus, vector< vector<haplo> >& population, vector< vector< vector<double> > >& afs) {
    afs.clear();
    vector<double> a (4,0);
    for (int g=0;g<population.size();g++) {
        vector< vector<double> > b;
        for (int j=0;j<consensus[g].size();j++) {
            b.push_back(a);
        }
        afs.push_back(b);
    }
    // cout << afs.size();
    //Generate provisional AFS from the population
    double ind=1/(size+0.);
    for (int j=0;j<population.size();j++) {
        for (int k=0;k<population[j].size();k++) {
            for (int l=0;l<population[j][k].seq.size();l++) {
                if (population[j][k].seq[l].loc<afs[j].size()) {
                    afs[j][population[j][k].seq[l].loc][population[j][k].seq[l].nuc]=afs[j][population[j][k].seq[l].loc][population[j][k].seq[l].nuc]+ind;
                }
            }
        }
    }
    
    for (int i=0;i<afs.size();i++) {
        for (int j=0;j<afs[i].size();j++) {
            double tot=0;
            for (int k=0;k<4;k++) {
                tot=tot+afs[i][j][k];
            }
            tot=1-tot;
            for (int k=0;k<4;k++) {
                if (consensus[i][j]==k) {
                    afs[i][j][k]=tot;
                }
            }
        }
    }
}

void ProcFreqs (run_params p, vector< vector< vector<double> > >& afs) {
    for (int i=0;i<afs.size();i++) {
        for (int j=0;j<afs[i].size();j++) {
            //Find max and impose min
            double tot=0;
            int index=-1;
            double max=0;
            for (int k=0;k<4;k++) {
                if (afs[i][j][k]<p.min_q) {
                    afs[i][j][k]=p.min_q;
                }
                if (afs[i][j][k]>max) {
                    max=afs[i][j][k];
                    index=k;
                }
                tot=tot+afs[i][j][k];
            }
            afs[i][j][index]=afs[i][j][index]-(tot-1);
        }
    }
}

void CalcDist (vector< vector< vector<double> > >& afs, vector< vector< vector<double> > >& afs_init) {
    double diff=0;
    for (int i=0;i<afs.size();i++) {
        for (int j=0;j<afs[i].size();j++) {
            double d=0;
            for (int k=0;k<4;k++) {
                d=d+abs(afs[i][j][k]-afs_init[i][j][k]);
            }
            diff=diff+(d/2);
        }
    }
    cout << "Distance " << diff << "\n";
}
