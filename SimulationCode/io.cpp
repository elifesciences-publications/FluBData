#include "shared.h"
#include <iostream>
#include <string>
#include <sstream>

void GetVariantData (vector< vector < vector<double> > >& init_var) {
	for (int j=1;j<=8;j++) {
		vector< vector<double> > vars;
		stringstream ss_j;
		ss_j << j;
		string s_j = ss_j.str();
		string s="../Data/Seg"+s_j+"/Variants0.out";
		cout << s << "\n";
		ifstream in_file;
		in_file.open(s.c_str());
		int x;
		do {
			vector<double> var;
			if (!(in_file >> x)) break;
			var.push_back(x);
			if (!(in_file >> x)) break;
			var.push_back(x);
			if (!(in_file >> x)) break;
			var.push_back(x);
			if (!(in_file >> x)) break;
			var.push_back(x);
			if (!(in_file >> x)) break;
			var.push_back(x);
			if (!(in_file >> x)) break;
			var.push_back(x);
			vector<double> v;
			double q=(var[1]+0.)/(var[5]+0.);
			v.push_back(q);
			q=(var[2]+0.)/(var[5]+0.);
			v.push_back(q);
			q=(var[3]+0.)/(var[5]+0.);
			v.push_back(q);
			q=(var[4]+0.)/(var[5]+0.);
			v.push_back(q);
			vars.push_back(v);
		} while (1==1);
		init_var.push_back(vars);
	}
}

void GetRandomVariantData (vector< vector < vector<double> > >& init_var, gsl_rng *rgen) {
	for (int j=1;j<=8;j++) {
		vector< vector<double> > vars;
		stringstream ss_j;
		ss_j << j;
		string s_j = ss_j.str();
        
        //Note : These lines are dependent upon the input data: How many variants files are available to the code will depend upon the dataset
		int v=7;
		while (v==7||v==6) {
			 v=floor(gsl_rng_uniform(rgen)*9); 
		}
        //End of note
        
		stringstream ss_k;
		ss_k << v;
		string v_s = ss_k.str();
		string s="../Data/Seg"+s_j+"/Variants"+v_s+".out";
		cout << s << "\n";
		ifstream in_file;
		in_file.open(s.c_str());
		int x;
		do {
			vector<double> var;
			if (!(in_file >> x)) break;
			var.push_back(x);
			if (!(in_file >> x)) break;
			var.push_back(x);
			if (!(in_file >> x)) break;
			var.push_back(x);
			if (!(in_file >> x)) break;
			var.push_back(x);
			if (!(in_file >> x)) break;
			var.push_back(x);
			if (!(in_file >> x)) break;
			var.push_back(x);
			vector<double> v;
            //Note : Variants are only read in where the read depth is greater than zero
			if (var[5]>0) {
				double q=(var[1]+0.)/(var[5]+0.);
				v.push_back(q);
				q=(var[2]+0.)/(var[5]+0.);
				v.push_back(q);
				q=(var[3]+0.)/(var[5]+0.);
				v.push_back(q);
				q=(var[4]+0.)/(var[5]+0.);
				v.push_back(q);
			} else {
				v.push_back(1);
				v.push_back(0);
				v.push_back(0);
				v.push_back(0);
			}
			vars.push_back(v);
		} while (1==1);
		init_var.push_back(vars);
	}
}


void PrintInitialFreqs (vector< vector< vector<double> > >& afs_init) {
	for (int i=0;i<afs_init.size();i++) {
		stringstream ss_j;
		ss_j << i;
		string s_i = ss_j.str();
		string s="Seg_"+s_i+"_Initial.out";
		
		ofstream afs_file;
		afs_file.open(s.c_str());
		for (int j=0;j<afs_init[i].size();j++) {
			afs_file << j << " ";
			for (int k=0;k<4;k++) {
				afs_file << afs_init[i][j][k] << " ";
			}
			afs_file << "\n";
		}
	}
}

void PrintRepFreqs (int rep, vector< vector< vector<double> > >& afs) {
    for (int i=0;i<afs.size();i++) {
        stringstream ss_i;
        stringstream ss_j;
        ss_i << i;
        string s_i = ss_i.str();
        ss_j << rep;
        string s_j = ss_j.str();
        string s="Seg_"+s_i+"_"+s_j+"_Final.out";
        cout << s << "\n";
        ofstream afs_file;
        afs_file.open(s.c_str());
        for (int j=0;j<afs[i].size();j++) {
            afs_file << j << " ";
            for (int k=0;k<4;k++) {
                afs_file << afs[i][j][k] << " ";
            }
            afs_file << "\n";
        }
    }
}

//Routine used in data processing
void ReadInitFreqs (vector< vector< vector<double> > >& afs){
    for (int i=0;i<8;i++) {
        vector< vector<double> > af;
        stringstream ss_j;
        ss_j << i;
        string s_i = ss_j.str();
        string s="Seg_"+s_i+"_Initial.out";
        ifstream afs_file;
        afs_file.open(s.c_str());
        int n;
        double x;
        for (int j=0;j<1000000;j++) {
            vector<double> a;
            if(!(afs_file >> n)) break;
            if(!(afs_file >> x)) break;
            a.push_back(x);
            if(!(afs_file >> x)) break;
            a.push_back(x);
            if(!(afs_file >> x)) break;
            a.push_back(x);
            if(!(afs_file >> x)) break;
            a.push_back(x);
            af.push_back(a);
        }
        afs.push_back(af);
        afs_file.close();
    }
}

void ReadRepFreqs (int rep, vector< vector< vector<double> > >& afs){
    for (int i=0;i<8;i++) {
        vector< vector<double> > af;
        stringstream ss_i;
        stringstream ss_j;
        ss_i << i;
        string s_i = ss_i.str();
        ss_j << rep;
        string s_j = ss_j.str();
        string s="Seg_"+s_i+"_"+s_j+"_Final.out";
        ifstream afs_file;
        afs_file.open(s.c_str());
        int n;
        double x;
        for (int j=0;j<1000000;j++) {
            vector<double> a;
            if(!(afs_file >> n)) break;
            if(!(afs_file >> x)) break;
            a.push_back(x);
            if(!(afs_file >> x)) break;
            a.push_back(x);
            if(!(afs_file >> x)) break;
            a.push_back(x);
            if(!(afs_file >> x)) break;
            a.push_back(x);
            af.push_back(a);
        }
        afs.push_back(af);
        afs_file.close();
    }
}


