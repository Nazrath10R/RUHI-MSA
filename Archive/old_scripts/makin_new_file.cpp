
include <iostream>
#include <fstream>
#include <vector>
#include "string"
#include <limits>
#include <algorithm>
#include <math.h>
#include <chrono>
#include <tuple>
#include <sstream>
#include <cstdlib>

using namespace std;

float* openFile(string path, int length) {
	ifstream inputFile;
	inputFile.open(path);
	float *array = new float[length];


	// put into function
	// if file is not found
	if (!inputFile) {
		cout << "no file" << "\n" << "double check: " << path << "\n";
		exit(-1);
	}

	// masses array
	int count = 0;
	while (!inputFile.eof()) {
		inputFile >> array[count];
		count++;
//        cout << masses[count-1]<<endl;
	}

int main(int argc, char** argv) {
	
	string input_path {"/Users/pedrocardoso/Documents/PTMs_project/RUHI-MSA/input/psm_Crystall.txt"};
	
	string output_file_path = "/Users/pedrocardoso/Documents/PTMs_project/RUHI-MSA/input/psm_Crystall_changed.txt"
	
	
	
	
	ifstream inputPeptideData;
	inputPeptideData.open(input_path);
	

	string peptide;
	float mass_shift;
	float peptide_mass;

	if (!inputPeptideData) { //testing if the file is being open correctly
		throw std::runtime_error("Error opening file");
	}
	else
	while (!inputPeptideData.eof()) {
		//================================//
		
		inputPeptideData >> peptide >> mass_shift >> peptide_mass;
		
		
		
		
		
		ofstream output_file {};
		output_file.open(output_file_path,  ios::out | ios::app);

		output_file << "\n" << peptide << "\n" << mass_shift << "\n" << peptide_mass;
		output_file.close();
		
	}