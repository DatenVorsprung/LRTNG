/*
 * processConfigModel.h
 *
 ** Lagrangian Reachtubes: The next Generation
 *
 *      Authors: Sophie Gruenbacher
 *      Contact: sophie.gruenbacher@tuwien.ac.at
 */

#ifndef PROCESS_INPUTFILES_H_
#define PROCESS_INPUTFILES_H_

#include <iostream>
#include <fstream>
#include <string>
#include "stdlib.h"
#include <math.h>

/* IBEX */
#include "ibex.h"

using namespace std;
using namespace ibex;

//process model and config file
IntervalVector processInput(
    double &totalTime, double &h, int &order, int &dim, double &rad0, int &exactVars,
    bool &computeVolume, bool &timeModel, bool &output, string &fdyn_file,
    string model){
	string ignoreLine;
	string parameters;
	ifstream fin;

/******** Process config file ********/
	string directory = "./Benchmarks/";

	string prefix_fdyn = "fdyn_";
    string prefix_init = "init_";
	string suffix_txt = ".txt";
    string suffix_fdyn = "_fdyn.txt";
    string suffix_init = "_init.txt";

    //TODO check if file exists!

	string model_file = directory + model + suffix_init;
	fdyn_file = directory + model + suffix_fdyn;

	cout << "Used file with initial values: " << model_file << endl;
	cout << "Used file with dynamics: " << fdyn_file << endl;


/******** Process model file ********/
	fin.open(model_file);

	fin >> model >> ignoreLine >> ignoreLine;//name of model
	cout << "Model: " << model << endl;

	/* get dimension of model and initial radius of box */
	fin >> dim >> ignoreLine >> ignoreLine;
	cout << "Dimension: " << dim << endl;

	fin >> rad0 >> ignoreLine >> ignoreLine;

	fin >> timeModel >> ignoreLine >> ignoreLine;
	cout << "time-model (boolean): " << timeModel << endl;

	fin >> exactVars >> ignoreLine >> ignoreLine;
	cout << "Control Variables with rad=0 (excluding time): " << exactVars << endl;

	if(timeModel)
		exactVars += 1;

	fin >> parameters;

	/* initial center */
	double temp_cx_i; //temp value for center of a variable‚
	IntervalVector cx(dim);

	for(int i=0; i<dim; i++){
	    	fin >> temp_cx_i;
	    	cx[i]=Interval(temp_cx_i);
	}

	cout << "Initial center: " << endl << endl;
    cout << parameters << endl << endl;
	cout << cx << endl << endl;

	fin.close();

	return cx;
}

//get symbolic differentiation from fdyn acording to integration order
void symbolicDifferentiation(Function &fdyn, std::vector<ibex::Function> &d2fdyn,
							 std::vector<std::vector<ibex::Function>> &d3fdyn,
							 std::vector<std::vector<std::vector<ibex::Function>>> &d4fdyn,
							 std::vector<std::vector<std::vector<std::vector<ibex::Function>>>> &d5fdyn, int &order, int &dim)
{
	 for(int i=0; i<dim; i++){
	    	ibex::Function df_i(fdyn[i], Function::DIFF);
	    	ibex::Function d2f_i(df_i, Function::DIFF);
	    	d2fdyn.push_back(d2f_i);

	        if(order > 1){
	        	std::vector<ibex::Function> d3f_i;
	        	std::vector<std::vector<ibex::Function>> d4f_i;
	        	std::vector<std::vector<std::vector<ibex::Function>>> d5f_i;
	        	for(int j=0; j<dim; j++){
	        		ibex::Function d2f_i_j(df_i[j], Function::DIFF);
	        		ibex::Function d3f_i_j(d2f_i_j, Function::DIFF);
	        		d3f_i.push_back(d3f_i_j);

	        		if(order > 2){
	        			std::vector<ibex::Function> d4f_i_j;
	        			std::vector<std::vector<ibex::Function>> d5f_i_j;
	        			for(int k=0; k<dim; k++){
	        				ibex::Function d3f_i_j_k(d2f_i_j[k], Function::DIFF);
	        				ibex::Function d4f_i_j_k(d3f_i_j_k, Function::DIFF);
	        				d4f_i_j.push_back(d4f_i_j_k);

	        				std::vector<ibex::Function> d5f_i_j_k;
	        				for(int l=0; l<dim; l++){
	        					ibex::Function d4f_i_j_k_l(d3f_i_j_k[l], Function::DIFF);
	        					ibex::Function d5f_i_j_k_l(d4f_i_j_k_l, Function::DIFF);
	        					d5f_i_j_k.push_back(d5f_i_j_k_l);
	        				}
	        				d5f_i_j.push_back(d5f_i_j_k);
	        			}
	        			d4f_i.push_back(d4f_i_j);
	        			d5f_i.push_back(d5f_i_j);
	        		}
	        	}
	        	d3fdyn.push_back(d3f_i);
	        	d4fdyn.push_back(d4f_i);
	        	d5fdyn.push_back(d5f_i);
	        }
	    }
}

#endif /* PROCESS_INPUTFILES_H_ */
