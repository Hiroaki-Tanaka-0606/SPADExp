// calcPSF main program

// include libraries from include path
#include <iostream>
#include <unistd.h>
#include <cstring>

// include libraries from current directory
#include "variables.hpp"
#include "setup.hpp"
#include "load_input.hpp"
#include "validation_Thomas_Fermi.hpp"
#include "calculation_Thomas_Fermi.hpp"
#include "validation_atomic_wfn.hpp"
#include "calculation_atomic_wfn.hpp"
#include "validation_psf.hpp"
#include "calculation_psf.hpp"
#include "log.hpp"

//namespace
using namespace std;

// main function
int main(int argc, const char** argv){
	cout << "Starting calcPSF program..." << endl;
	
	// stdin, stdout, stderr are terminal or not (file)
	int stdin_is_terminal=isatty(fileno(stdin));
	int stdout_is_terminal=isatty(fileno(stdout));
	int stderr_is_terminal=isatty(fileno(stderr));

	cout << "----IO information----" << endl;
	cout << "stdin is " << (stdin_is_terminal ? "" : "NOT ") << "a terminal" << endl;
	cout << "stdout is " << (stdout_is_terminal ? "" : "NOT ") << "a terminal" << endl;
	cout << "stderr is " << (stderr_is_terminal ? "" : "NOT ") << "a terminal" << endl;
	cout << endl;

	// stdin should not be a terminal
	if(stdin_is_terminal==1){
		cout << "Error: stdin should be a file" << endl;
		printf("Usage: %s < input_file\n", argv[0]);
		return 0;
	}

	// initialize: load default values and constants
	initialize();
	
	// load input
	cout << "----Load input----" << endl;
	if(load_input()!=1){
		goto FINALIZATION;
	}
	cout << endl;

	// open output file
	if(Output_file_set==false){
		cout << "Error: invalid Output_file" << endl;
		goto FINALIZATION;
	}
	Output_file_obj=fopen(Output_file, "w");
	if(Output_file_obj==NULL){
		cout << "Error: could not open the output file" << endl;
		goto FINALIZATION;
	}
	
	// open log file if necessary
	if(Log_file_set){
		Log_file_obj=fopen(Log_file, "w");
		if(Log_file_obj==NULL){
			cout << "Error: could not open the log file" << endl;
			goto FINALIZATION;
		}
	}

	// perform calculation
	cout << "----Perform calculation----" << endl;
	int validation_result;
	if(strcmp(Calculation, "Thomas-Fermi")==0){
		// Thomas-Fermi potential calculation
		write_log((char*)"----Thomas-Fermi potential calculation----");
		validation_result=validation_Thomas_Fermi();
		if(validation_result==1){
			if(TF_test){
				test_Thomas_Fermi();
			}else{
				calc_Thomas_Fermi();
			}
		}
	}else if(strcmp(Calculation, "Atomic-wfn")==0){
		// Atomic wavefunction calculation
		write_log((char*)"----Atomic wavefunction calculation----");
		validation_result=validation_atomic_wfn();
		if(validation_result==1){
			sequence_atomic_wfn();
		}
	}else if(strcmp(Calculation, "PSF")==0){
		// Photoemission structure factor calculation
		write_log((char*)"----Photoemission structure factor calculation----");
		validation_result=validation_PSF();
		if(validation_result==1){
			PSF_test();
		}
	}else{
		cout << "Error: invalid calculation" << endl;
		goto FINALIZATION;
	}

	goto FINALIZATION;
	
 FINALIZATION:
	// finalize: free memory allocated by "new" in initialize()
	finalize();
}
