
/*  main.cpp
    RUHI-MSA

    Created by Nazrath Nawaz on 28/08/2019.
    Developed by Pedro Cardoso on 09/06/2020.
    Copyright © 2019 Nazrath Nawaz. All rights reserved.
*/
#include "main.h"


//=======================================================================================================//
using namespace std;


int main(int argc, char** argv) {

	/* Main arguments:
	 * first = sample name
	 * second = tolerance
	 */ 
	auto start = chrono::steady_clock::now();
	float tolerance {};
	
	

	//setting user argument options (file to analyse and tolerance)
	string sample {};
	if (argc ==1) {
		cout << "at least one argument (sample name) is necessary" << endl;
		exit(1);
	}
	else if (argc < 3) {
		tolerance = 10.0; // tolerance default is 10
		sample = argv[1];
	} 
	else if (argc < 4) {
		stringstream convert {argv[2]};
		if (!(convert >> tolerance)) { // converting agrv[2] to float
			cout << "Please insert a valid tolerance (must be an integer)" << endl;
		}
		sample = argv[1];
	}
	else {
		cout << "This program only takes 2 arguments" << endl;
		exit(1);
	}


	
//printing logo	
	ifstream f("../logo.txt");
    if (f.is_open())
        std::cout << f.rdbuf();

	cout << endl << "Tolerance: " << tolerance << " \nFilename: " << sample << endl;
	//--------------------------Setting paths------------------------------------------//
	
	//string sample = "JOHNPP_l1";
	string dir_path = "../";
	string input_path = dir_path + "input/" + sample + ".txt";
	string libs = "libs/"; 
	//float tol = tolerance;
	//string tole = to_string(tol); //just for testing different tolerances


	//================================| Input |================================//
	
                                                       
	
	//// PTM masses and probabilities

	// mass shift array
	const int PTM_list = 1498;  //1497 PTMs
	const int list_end = PTM_list; //last element to create null prob for no matches. Prob list will have one more element
	const int nat_list = 90; //89 natural PTMs
	const int nat_end = nat_list;

	//// UniMod masses and probabilities

	// mass shift array // PTMs array
	float* masses = openFile(dir_path + libs + "masses.txt", PTM_list);

	// Probabilities array
	float* prob = openFile(dir_path + libs + "prob.txt", (PTM_list+1));

	// null probability for no matches (last new vector position)
	prob[list_end] = 0;
	//cout << masses[0] << " and " << masses[PTM_list-1] << " and " << prob[PTM_list-1] << " " << prob[PTM_list] << endl;
	
	//// Natural frequency masses and probabilities

	// Nat probabilities array
	float* nat_prob = openFile(dir_path + libs + "nat_prob.txt", (nat_list+1));

	nat_prob[nat_end] = 0;

	// Nat masses array
	float* nat_masses = openFile(dir_path + libs + "nat_masses.txt", nat_list);


	// Nat names array
	string nat_names[nat_list];
	ifstream fakeNames;
	fakeNames.open(dir_path + libs + "nat_names.txt");

	int names_count = 0;
	while (!fakeNames.eof()) {
		fakeNames >> nat_names[names_count];
		names_count++;
	}

	// sizes
	
	int nat_prob_length = nat_list;

	//-------------------------------------------------------//

	//-------------------------OUTPUT file--------------------------------//


	string output_file_path = dir_path + "results/" + sample + "_output.txt";
	
	ofstream output_file;
	output_file.open(output_file_path, ios::out);

	output_file << "Peptide" << "\t" << "Mass_shift" << "\t" << "Peptide Mass" << "\t" << "Mass of PTMs" << "\t" << "Score" << "\t" <<  "PTM1" << "\t"<< "PTM2" << "\t"<< "PTM3" << "\t"<< "PTM4" << "\t" << "PTM5";
	output_file.close();
	
	
	//================================| Peptide - Mass Shift |================================//
	//reading Peptide data: sequence, Mass_shift and peptide Mass
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
		

        // Defining mass error

        //float tolerance = 10; // this has been defined on the  .main arguments
        float mass_error = peptide_mass * tolerance / 1000000;
        //cout << mass_error << endl;
        float mass_shift_upper = mass_shift + mass_error;
        float mass_shift_lower = mass_shift - mass_error;
        //cout << mass_shift_lower << "-" << mass_shift_upper << endl;
        mass_shift_upper = roundf(mass_shift_upper*10000)/10000; //check if this works without rounding. If we round we are losing precision
        mass_shift_lower = roundf(mass_shift_lower*10000)/10000;
		//cout << mass_shift_lower << "-" << mass_shift_upper << endl;
        

		//================================//

		// counters
		int combination = 0;
		int more_than_3 = 0;
		//=======================================================================================================//


		// null vector of probability positions to populate
		vector<int> highest_prob_pos (5,list_end);
		//	vector<int> *highest_prob_pos_1 = new vector<int>();

		// number of loops to attempt counter
		int number_of_possible_ptms = 1;
		//cout << "Searching combinations of:" << "\n";

		//--------------------------------------------------------------------------------------------------//

		//––––––––––––––––––––––––––//
		//––––––––| 1 PTM |–––––––––//
		//––––––––––––––––––––––––––//

		//Search PTMs matching mass shift range
		auto [highest_prob_pos_1_ptr, position_combination1_ptr, pos_length] = finding1PTM(masses, prob, PTM_list, mass_shift_lower, mass_shift_upper);

		int combination_1 = *pos_length;
		float prob_average_1 = 0.0;
		
		// storing highest probability and index of matched PTM 
		if (combination_1 > 0) {
			prob_average_1 = prob[(*highest_prob_pos_1_ptr).at(0)];
			//cout << "\n" <<  "The mass shift can be explained by 1 PTM  - ";
			//cout << "with score: " << prob_average_1 << "\n";
			highest_prob_pos[0] = (*highest_prob_pos_1_ptr).at(0);
			combination = combination_1;
			more_than_3 = 1;

		} else { number_of_possible_ptms++; } 


		//--------------------------------------------------------------------------------------------------//
		//auto end2 = chrono::steady_clock::now();
		//––––––––––––––––––––––––––//
		//––––––––| 2 PTMs |––––––––//
		//––––––––––––––––––––––––––//
		
		//    	Search combination of 2 and 3 PTMs that match mass shift
		
		auto [highest_prob_pos_2, position_combination2, combination_2, highest_prob_pos_3, combination_3] = findingCombinationOf2and3PTM(masses, prob, PTM_list, mass_shift_lower, mass_shift_upper);
		
		
		//auto [highest_prob_pos_2, combination_2, highest_prob_pos_3, combination_3, highest_prob_pos_4, combination_4, highest_prob_pos_5, combination_5] = findingCombinationOf2_3_4_5_PTM(masses, prob, PTM_list, mass_shift_lower, mass_shift_upper);
		
//    	Search combination of 2 PTMs that match mass shift
		
		//auto [highest_prob_pos_2, position_combination2, combination_2] = findingCombinationOf2PTM(masses, prob, PTM_list, mass_shift_lower, mass_shift_upper);

		float prob_average_2 = 0.0;

		// storing highest probability matches and its index
		if (combination_2 > 0) {

			prob_average_2 = (prob[highest_prob_pos_2[0]] + prob[highest_prob_pos_2[1]]) / 2;
			//cout << "\n" << "The mass shift can be explained by 2 PTMs - ";
			//cout << "with score: " << prob_average_2 << "\n";
			if (prob_average_2 > prob_average_1) {
				//cout << "2 PTMs have a higher weighting" << "\n";
				highest_prob_pos[0] = highest_prob_pos_2[0]; highest_prob_pos[1] = highest_prob_pos_2[1];
				combination = combination_2;
				more_than_3 = 2;
			}
		} else { number_of_possible_ptms++; }


		//--------------------------------------------------------------------------------------------------//

		//––––––––––––––––––––––––––//
		//––––––––| 3 PTMs |––––––––//
		//––––––––––––––––––––––––––//

//     Search combination of 3 PTMs that match mass shift

		//auto [highest_prob_pos_3_ptr, combination_3_ptr] =findingCombinationOf3PTM ( masses ,prob , PTM_list, mass_shift_lower, mass_shift_upper);

		float prob_average_3 = 0.0;
		// storing highest probability matches and its index
		if (combination_3 > 0) {

			prob_average_3 = (prob[highest_prob_pos_3[0]] + prob[highest_prob_pos_3[1]] + prob[highest_prob_pos_3[2]]) / 3;
			//cout << "\n" << "The mass shift can be explained by 3 PTMs - ";
			//cout << "with score: " << prob_average_3 << "\n";

			if ( prob_average_3 > prob_average_1 && prob_average_3 > prob_average_2 ) {
				//cout << "3 PTMs have a higher weighting" << "\n";
				highest_prob_pos[0] = highest_prob_pos_3[0]; highest_prob_pos[1] = highest_prob_pos_3[1];  highest_prob_pos[2] = highest_prob_pos_3[2];
				combination = combination_3;
				more_than_3 = 3;
			}
		} else { number_of_possible_ptms++; }


		//--------------------------------------------------------------------------------------------------//

		// highest values of all probabilities till 3
		vector<float> prob_vector_1 = {prob_average_1, prob_average_2, prob_average_3};
		float prob_vector_max_1 = *max_element(prob_vector_1.begin(), prob_vector_1.end());

//		auto end3 = chrono::steady_clock::now();
		//--------------------------------------------------------------------------------------------------//

		//––––––––––––––––––––––––––//
		//––––––––| 4 PTMs |––––––––//
		//––––––––––––––––––––––––––//

//    Search for combination of 4 natural PTMs that have an higher probability

		//float sum_of_probs_4 = prob_vector_max_1 * 4; // prob max is an average so we need to find the actual value of the probs before 
		//auto [highest_prob_pos_4, combination_4] = findingCombinationOf4NaturalPTM ( nat_prob, nat_masses, prob, nat_prob_length, mass_shift_lower, mass_shift_upper, sum_of_probs_4);
	//--------------------------------------------------------------------------------------------------//
//		float prob_average_4 = 0.0;
		// storing highest probability matches and its index
//		if (combination_4 > 0) {
//
//			prob_average_4 = ((nat_prob[highest_prob_pos_4[0]] + nat_prob[highest_prob_pos_4[1]] + nat_prob[highest_prob_pos_4[2]] + nat_prob[highest_prob_pos_4[3]]) / 4);
//			cout << "\n" << "The mass shift can be explained by 4 PTMs - ";
//			cout << "with score: " << prob_average_4 << "\n";
//
//			if ( prob_average_4 > prob_average_3 && prob_average_4 > prob_average_2 && prob_average_4 > prob_average_1) {
//				cout << "  >>> 4 PTMs have a higher weighting  <<< " << "\n" << "\n"; // this will always happen because we are only checking the combinations with higher prob than the rest 
//				highest_prob_pos[0] = highest_prob_pos_4[0]; highest_prob_pos[1] = highest_prob_pos_4[1];
//				highest_prob_pos[2] = highest_prob_pos_4[2]; highest_prob_pos[3] = highest_prob_pos_4[3];
//				combination = combination_4;
//				more_than_3 = 4;
//			}
//
//		} else { number_of_possible_ptms++; }
//
//
// highest values of all probabilities
//vector<float> prob_vector_2 = {prob_average_1, prob_average_2, prob_average_3, prob_average_4};
//float prob_vector_max_2 = *max_element(prob_vector_2.begin(), prob_vector_2.end());


		//--------------------------------------------------------------------------------------------------//

		//––––––––––––––––––––––––––//
		//––––––––| 4 and 5 Natural PTMs |––––––––//
		//––––––––––––––––––––––––––//
		float sum_of_probs_4 = prob_vector_max_1 * 4; // prob max is an average so we need to find "raw" probablility.
		float prob_average_4 = 0.0;
//    Search for combination of 4 and/or 5 natural PTMs that have higher probability
		auto [highest_prob_pos_4, combination_4, highest_prob_pos_5, combination_5, new_prob_average_4] = findingCombinationOf4and5NaturalPTM (nat_prob,nat_masses,prob, nat_prob_length, mass_shift_lower, mass_shift_upper, sum_of_probs_4, prob_average_1, prob_average_2, prob_average_3, prob_average_4, prob_vector_max_1);

		prob_average_4 = new_prob_average_4;
		//---------------------------------------------------------------------------------------------------------//
		
		// storing highest probability matches and its index for combination of 4
		if (combination_4 > 0) {
			//prob_average_4 has been calculated during findingCombinationOf4and5NaturalPTM function.
			//prob_average_4 = (nat_prob[highest_prob_pos_4[0]] + nat_prob[highest_prob_pos_4[1]] + nat_prob[highest_prob_pos_4[2]] + nat_prob[highest_prob_pos_4[3]]) / 4;
	
			//cout << "\n" << "The mass shift can be explained by 4 PTMs - ";
			//cout << "with score: " << prob_average_4 << "\n";

			if ( prob_average_4 > prob_average_3 && prob_average_4 > prob_average_2 && prob_average_4 > prob_average_1) {
				//cout << "  >>> 4 PTMs have a higher weighting  <<< " << "\n" << "\n"; // this will always happen because we are only checking the combinations with higher prob than the rest 
				highest_prob_pos[0] = highest_prob_pos_4[0]; highest_prob_pos[1] = highest_prob_pos_4[1];
				highest_prob_pos[2] = highest_prob_pos_4[2]; highest_prob_pos[3] = highest_prob_pos_4[3];
				combination = combination_4;
				more_than_3 = 4;
			}

		} else { number_of_possible_ptms++; }
		//--------------------------------------------------------------------------------------------------//

		float prob_average_5 = 0.0;
		// storing highest probability matches and its index for combination of 5
		if (combination_5 > 0) {

			prob_average_5 = (nat_prob[highest_prob_pos_5[0]] + nat_prob[highest_prob_pos_5[1]] + nat_prob[highest_prob_pos_5[2]] + nat_prob[highest_prob_pos_5[3]] + nat_prob[highest_prob_pos_5[4]]) / 5;

			//cout << "\n" << "The mass shift can be explained by 5 PTMs - ";
			//cout << "with score: " << prob_average_5 << "\n";


			if (prob_average_5 > prob_average_4 && prob_average_5 > prob_average_3 && prob_average_5 > prob_average_2 && prob_average_5 > prob_average_1) {
				//cout << "  >>> 5 PTMs have a higher weighting  <<< " << "\n" << "\n";
				highest_prob_pos[0] = highest_prob_pos_5[0]; highest_prob_pos[1] = highest_prob_pos_5[1];
				highest_prob_pos[2] = highest_prob_pos_5[2]; highest_prob_pos[3] = highest_prob_pos_5[3]; highest_prob_pos[4] = highest_prob_pos_5[4];
				combination = combination_5;
				more_than_3 = 5;
			}

		} else { number_of_possible_ptms++; }
		//auto end4 = chrono::steady_clock::now();
		//------------------------------------------------------------------------------//


		//––––––––––––––––––––––––––//
		//––––––| 4 PTMs alt |––––––//
		//––––––––––––––––––––––––––//

		if (combination == 0 && number_of_possible_ptms == 6) {

//        cout << "4 PTMs (Alternative)...";
			//
			//auto [highest_prob_pos_4, combination_4] = findingCombinationOf4PTM_inprogress(masses, prob, PTM_list, mass_shift_lower, mass_shift_upper);
			
			//auto [highest_prob_pos_4, combination_4, highest_prob_pos_5, combination_5] = findingCombinationOf4and5PTM(masses, prob, PTM_list, mass_shift_lower, mass_shift_upper);
			
			auto [highest_prob_pos_4, combination_4, highest_prob_pos_5, combination_5] = findingCombinationOf4and5PTM_Best(masses, prob, PTM_list, mass_shift_lower, mass_shift_upper);
			
			number_of_possible_ptms = 4;
			// cout << "4 PTMs...";



			float prob_average_4 = 0.0;

			if (combination_4 > 0) {

				prob_average_4 = (prob[highest_prob_pos_4[0]] + prob[highest_prob_pos_4[1]] + prob[highest_prob_pos_4[2]] + prob[highest_prob_pos_4[3]]) / 4;

				//cout << "\n" << "The mass shift can be explained by 4 PTMs - ";
				//cout << "with score: " << prob_average_4 << "\n";

				if ( prob_average_4 > prob_average_3 && prob_average_4 > prob_average_2 && prob_average_4 > prob_average_1) {
					//cout << "  >>> 4 PTMs have a higher weighting  <<< " << "\n" << "\n";
					highest_prob_pos[0] = highest_prob_pos_4[0]; highest_prob_pos[1] = highest_prob_pos_4[1];
					highest_prob_pos[2] = highest_prob_pos_4[2]; highest_prob_pos[3] = highest_prob_pos_4[3];
					combination = combination_4;
					more_than_3 = 6;
				}

				//		break;

			} else {
				number_of_possible_ptms = 6;
			}
			
			
			float prob_average_5 = 0.0;

			if (combination_5 > 0) {

				prob_average_5 = (prob[highest_prob_pos_5[0]] + prob[highest_prob_pos_5[1]] + prob[highest_prob_pos_5[2]] + prob[highest_prob_pos_5[3]] +prob[highest_prob_pos_5[4]]) / 5;

				//cout << "\n" << "The mass shift can be explained by 4 PTMs - ";
				//cout << "with score: " << prob_average_4 << "\n";

				if ( prob_average_5 > prob_average_4 && prob_average_5 > prob_average_3 && prob_average_5 > prob_average_2 && prob_average_5 > prob_average_1) {
					//cout << "  >>> 4 PTMs have a higher weighting  <<< " << "\n" << "\n";
					highest_prob_pos[0] = highest_prob_pos_5[0]; highest_prob_pos[1] = highest_prob_pos_5[1];
					highest_prob_pos[2] = highest_prob_pos_5[2]; highest_prob_pos[3] = highest_prob_pos_5[3]; highest_prob_pos[4] = highest_prob_pos_5[4];
					combination = combination_5;
					more_than_3 = 6;
				}

				//		break;

			} else {
				number_of_possible_ptms = 6;
			}
		}
//		auto end5 = chrono::steady_clock::now();

		//=======================================================================================================//

		// No results
		if (combination == 0 && number_of_possible_ptms == 6) {
			cout << "\n" << "The mass shift of " <<  peptide << " can NOT be explained till 5 PTMs" << "\n" << "\n";

			//ofstream output_file;
			output_file.open(output_file_path, ios::out | ios::app);

			output_file << peptide << "\t" << mass_shift << "\t" << peptide_mass << "\t" << "" << "\t" << ""<< "\t" << ""<< "\t" << ""<< "\t" << "" << "\n";
			output_file.close();


			//// skip to next peptide
			continue;
		}

		//=======================================================================================================//


		//================================| Output |================================//

		//// print output combination stats
//		cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << "\n\n";
//		cout << "peptide:  " << peptide << "\n\n";
//		cout << "PTM list:   " << length << ";  number of PTMs: " << number_of_possible_ptms << "\n";
//		cout << "mass shift: " << mass_shift << " ± " << mass_error << ";  combinations: " << combination << "\n" << "\n";
//		cout << "=======================================================================" << "\n" << "\n";


//      Getting names of founded PTMs

		vector<string> out_PTM_name;
		float sum_of_masses {}; // mass of all PTMS
		float prob_PTMs {}; // score of PTMs
		//// 1-3 PTMs
		if (more_than_3 <= 3) {
			//// print PTM names from file
			fstream inputPTMnames;
			inputPTMnames.open(dir_path + libs + "names.txt");
			//string sline;
			
			
			for (int i = 0; i < 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					GotoLine(inputPTMnames, 1 + highest_prob_pos[i]);
					string line;
					inputPTMnames >> line;
					//cout << line << "\t";
					out_PTM_name.push_back(line);
					sum_of_masses += masses[highest_prob_pos[i]];
					prob_PTMs += prob[highest_prob_pos[i]];
				}
			}


		//// 4-5 PTMs
		} else if (more_than_3 <= 5) {
			
			fstream Natural_PTM_Names;
			Natural_PTM_Names.open(dir_path + libs + "nat_names.txt");
			//cout << "\nPTM combination: ";
			
			
			for (int i = 0; i <= 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					GotoLine(Natural_PTM_Names, 1 + highest_prob_pos[i]);
					string line;
					Natural_PTM_Names >> line;
					//cout << line << "\t";
					out_PTM_name.push_back(line);
					sum_of_masses += nat_masses[highest_prob_pos[i]];
					prob_PTMs += nat_prob[highest_prob_pos[i]];
				}
			}


		} else if (more_than_3 == 6) {

			//// print PTM names from file
			fstream inputPTMnames;
			inputPTMnames.open(dir_path + libs+  "names.txt");
			//string sline;

			
			
			for (int i = 0; i <= 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					GotoLine(inputPTMnames, 1 + highest_prob_pos[i]);
					string line;
					inputPTMnames >> line;
					//cout << line << "\t";
					out_PTM_name.push_back(line);
					sum_of_masses += masses[highest_prob_pos[i]];
					prob_PTMs += prob[highest_prob_pos[i]];
				}
			}
		}


		// print run-time
		//cout << "\n\n" << "execution time = " << elapsed_time << "\n\n";
		//cout << "=======================================================================" << "\n";


		//--------------------------------------------------------------------------------------------------/


		//================================| Write Output |================================//

		//// write output into file
		//ofstream output_file;
		output_file.open(output_file_path, ios::out | ios::app); //ios::app to open file in append mode
		output_file << "\n" << peptide << "\t" << mass_shift << "\t" << peptide_mass << "\t" 
							 << sum_of_masses << "\t" << prob_PTMs << "\t";
		for (int i = 0; i < out_PTM_name.size(); i++) {
			output_file << out_PTM_name[i] << "\t";
		}
		output_file.close();


		//--------------------------------------------------------------------------------------------------/


		//================================| Alt Output |================================//
		
		//// write alternative outputs into files

//		string alt_output_file_path = dir_path + "results/" + sample + "/" + peptide + "_alternative.txt";
//
//		ofstream alt_output_file;
//		alt_output_file.open(alt_output_file_path, ios::out | ios::app);
//
//
//
//
//		if (combination_1 > 0 && prob_average_1 > 0) {
//
//			alt_output_file << "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t";
//
//			for (int i = 1; i < 1; i++) {
//				alt_output_file << i << "\t";
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
//			for (int i = 0; i < 1; i++) {
////                alt_output_file << names[highest_prob_pos_1[i]] << "\t";
//				alt_output_file << "\t";
//			}
//
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
//			for (int i = 0; i < 1; i++) {
//				if (highest_prob_pos_1[i] != list_end) {
//					alt_output_file << masses[highest_prob_pos_1[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
//			for (int i = 0; i < 1; i++) {
//				if (highest_prob_pos_1[i] != list_end) {
//					alt_output_file << prob[highest_prob_pos_1[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
//			for (int i = 0; i < 1; i++) {
//				if (highest_prob_pos_1[i] != list_end) {
//					alt_output_file << 1 + highest_prob_pos_1[i] << "\t";
//				}
//			}
//
//
//		} else if (combination_2 > 0 && prob_average_2 > 0) {
//
//			alt_output_file << "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t";
//
//			for (int i = 1; i < 2; i++) {
//				alt_output_file << i << "\t";
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
//			for (int i = 0; i < 2; i++) {
////                alt_output_file << names[highest_prob_pos_1[i]] << "\t";
//				alt_output_file << "\t";
//			}
//
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
//			for (int i = 0; i < 2; i++) {
//				if (highest_prob_pos_2[i] != list_end) {
//					alt_output_file << masses[highest_prob_pos_2[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
//			for (int i = 0; i < 2; i++) {
//				if (highest_prob_pos_2[i] != list_end) {
//					alt_output_file << prob[highest_prob_pos_2[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
//			for (int i = 0; i < 2; i++) {
//				if (highest_prob_pos_2[i] != list_end) {
//					alt_output_file << 1 + highest_prob_pos_2[i] << "\t";
//				}
//			}
//
//		} else if (combination_3 > 0 && prob_average_3 > 0) {
//
//			alt_output_file << "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t";
//
//			for (int i = 1; i < 3; i++) {
//				alt_output_file << i << "\t";
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
//			for (int i = 0; i < 3; i++) {
//				// alt_output_file << names[highest_prob_pos_1[i]] << "\t";
//				alt_output_file << "\t";
//			}
//
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
//			for (int i = 0; i < 3; i++) {
//				if (highest_prob_pos_3[i] != list_end) {
//					alt_output_file << masses[highest_prob_pos_3[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
//			for (int i = 0; i < 3; i++) {
//				if (highest_prob_pos_3[i] != list_end) {
//					alt_output_file << prob[highest_prob_pos_3[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
//			for (int i = 0; i < 3; i++) {
//				if (highest_prob_pos_3[i] != list_end) {
//					alt_output_file << 1 + highest_prob_pos_3[i] << "\t";
//				}
//			}
//		}
//
//		else if ( more_than_3 == 6 && combination_4 > 0 && prob_average_4 > 0) {
//
//			alt_output_file << "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t";
//
//			for (int i = 1; i < 4; i++) {
//				alt_output_file << i << "\t";
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
//			for (int i = 0; i < 4; i++) {
//				// alt_output_file << names[highest_prob_pos_1[i]] << "\t";
//				alt_output_file << "\t";
//			}
//
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
//			for (int i = 0; i < 4; i++) {
//				if (highest_prob_pos_4[i] != list_end) {
//					alt_output_file << masses[highest_prob_pos_4[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
//			for (int i = 0; i < 4; i++) {
//				if (highest_prob_pos_4[i] != list_end) {
//					alt_output_file << prob[highest_prob_pos_4[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
//			for (int i = 0; i < 4; i++) {
//				if (highest_prob_pos_4[i] != list_end) {
//					alt_output_file << 1 + highest_prob_pos_4[i] << "\t";
//				}
//			}
//		}
//
//
//
//
//
//		if (more_than_3 == 4 && combination_4 > 0) {
//			alt_output_file <<  "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t" << "1" << "\t" << "2" << "\t" << "3" << "\t" << "4" << "\t";
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
//			for (int i = 0; i < 4; i++) {
//				if (highest_prob_pos_4[i] != list_end) {
//					alt_output_file << nat_names[highest_prob_pos_4[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
//			for (int i = 0; i < 4; i++) {
//				if (highest_prob_pos_4[i] != list_end) {
//					alt_output_file << nat_masses[highest_prob_pos_4[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
//			for (int i = 0; i < 4; i++) {
//				if (highest_prob_pos_4[i] != list_end) {
//					alt_output_file << nat_prob[highest_prob_pos_4[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
//			for (int i = 0; i < 4; i++) {
//				if (highest_prob_pos_4[i] != list_end) {
//					alt_output_file << 1 + highest_prob_pos_4[i] << "\t";
//				}
//			}
//
//		} else if (more_than_3 == 5 && combination_5 > 0) {
//			alt_output_file << "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t" << "1" << "\t" << "2" << "\t" << "3" << "\t" << "4" << "\t" << "5" << "\t";
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
//
//			for (int i = 0; i < 5; i++) {
//				if (highest_prob_pos_5[i] != list_end) {
//					alt_output_file << nat_names[highest_prob_pos_5[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
//			for (int i = 0; i < 5; i++) {
//				if (highest_prob_pos_5[i] != list_end) {
//					alt_output_file << nat_masses[highest_prob_pos_5[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
//			for (int i = 0; i < 5; i++) {
//				if (highest_prob_pos_5[i] != list_end) {
//					alt_output_file << nat_prob[highest_prob_pos_5[i]] << "\t";
//				}
//			}
//
//			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
//			for (int i = 0; i < 5; i++) {
//				if (highest_prob_pos_5[i] != list_end) {
//					alt_output_file << 1 + highest_prob_pos_5[i] << "\t";
//				}
//			}
//		}
//
//
//		alt_output_file << "\n";
//		alt_output_file.close();



	}
	
	delete [] masses; delete[]  prob; delete [] nat_masses; delete [] nat_prob; //deleting array created using memory allocation
	
	auto end = chrono::steady_clock::now();
		auto diff = end - start;
		string elapsed_time;

		if (chrono::duration <double, milli> (diff).count() < 1) {
			elapsed_time = to_string(chrono::duration <double, milli> (diff).count()) + " ms";
		} else {
			elapsed_time = to_string(chrono::duration <double, milli> (diff).count() / 1000) + " s";
		}
		
	
	
	output_file.open(output_file_path, ios::out | ios::app); 
	
	output_file << "\n" << "Elapsed time is : " << elapsed_time << "\n";
	output_file.close();
		
	return 0;

}



//a->push_back(3);
//cout << "here " << a->at(0);


