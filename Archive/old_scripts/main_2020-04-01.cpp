
//  main.cpp
//  ptm
//
//  Created by Nazrath Nawaz on 28/08/2019.
//  Copyright © 2019 Nazrath Nawaz. All rights reserved.
//



#include "main.h"

//=======================================================================================================//

//template<typename T, size_t n>
//void print_array(T const(& arr)[n])
//{
//    for (size_t i = 0; i < n; i++)
//        std::cout << arr[i] << ' ';
//}


using namespace std;


// function to extract line of file based on position specified
void GotoLine(fstream& file, unsigned int num) {
	file.seekg(ios::beg);
	for (int i = 0; i < num - 1; ++i) {
		file.ignore(numeric_limits<streamsize>::max(), '\n');
	}
}


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

	return array;
}


//=======================================================================================================//



int main( ) {
	auto start = chrono::steady_clock::now();
	string sample = "JOHNPP_l1";

	string dir_path = "/Users/pedrocardoso/Documents/RUHI-MSA/";
	string input_path = dir_path + "input/" + sample + ".txt";



	//-------------------------------- Read Input files ------------------------------------------//

	// number of lines in input file
	//
	int lines;
	string s;

	ifstream inputFile;
	inputFile.open(input_path);

	if (!inputFile) {
		cout << "no file" << "\n" << "double check in: " << input_path << "\n";
		exit(-1);
	}

	while (!inputFile.eof()) {
		getline(inputFile, s);
		lines++;
	}

    cout <<  "No of Peptides: " << lines << "\n";

	//--------------------------------//



	//================================| Input |================================//


	//// PTM masses and probabilities

	// mass shift array
	const int PTM_list = 1498; 
	const int list_end = PTM_list -1;
	const int nat_list = 90;
	const int nat_end = nat_list -1;

	//// UniMod masses and probabilities

	// mass shift array
	float* masses = openFile(dir_path + "new_masses.txt", PTM_list);

	// Probabilities array
	float* prob = openFile(dir_path + "new_prob.txt", PTM_list);

	// null probability for no matches (last new vector position)
	prob[list_end] = 0;


	//// Natural frequency masses and probabilities

	// Nat probabilities array
	float* nat_prob = openFile(dir_path + "new_nat_prob.txt", nat_list);

	nat_prob[nat_end] = 0;

	// Nat masses array
	float* nat_masses = openFile(dir_path + "new_nat_masses.txt", nat_list);


	// Nat names array
	string nat_names[nat_list];
	ifstream fakeNames;
	fakeNames.open(dir_path + "new_nat_names.txt");

	int names_count = 0;
	while (!fakeNames.eof()) {
		fakeNames >> nat_names[names_count];
		names_count++;
	}

	// sizes
	int length = PTM_list;
	int nat_prob_length = nat_list;

	//-------------------------------------------------------//




	//================================| Peptide - Mass Shift |================================//


	for (int x = 1; x <= lines; x += 3) {

		//================================//

		fstream inputPeptides;
		inputPeptides.open(input_path);
		GotoLine(inputPeptides, x);
		string peptide;
		inputPeptides >> peptide;

		fstream inputMassShift;
		inputMassShift.open(input_path);
		GotoLine(inputMassShift, x + 1);
		float mass_shift;
		inputMassShift >> mass_shift;

        fstream inputPeptideMass;
        inputPeptideMass.open(input_path);
        GotoLine(inputPeptideMass, x + 2);
        float peptide_mass;
        inputPeptideMass >> peptide_mass;

        
        // mass error
//        float peptide_mass = 1000;
        float tolerance = 10;
        float mass_error = peptide_mass * tolerance / 1000000;
        
        float mass_shift_upper = mass_shift + mass_error;
        float mass_shift_lower = mass_shift - mass_error;
        
        mass_shift_upper = roundf(mass_shift_upper*10000)/10000;
        mass_shift_lower = roundf(mass_shift_lower*10000)/10000;

        
		//================================//

		string output_file_path = dir_path + "results/" + sample + "/" + peptide + "_output.txt";

		// int break_scenario=0;

		//----------------------------------//
		// string peptide = "VSAMEDEMNEMK";
		// float mass_shift = -0.0011;
		//----------------------------------//

		//cout << "\n|------- " << (x / 3) + 1 << " of " << lines / 3 << " -------| \n\n";

		// counters
		//    int count;
		int combination = 0;
		int more_than_3 = 0;



		//=======================================================================================================//

		// start timer
		//auto start = chrono::steady_clock::now();

		// null vector of probability positions to populate

		vector<int> highest_prob_pos (5, list_end);
		//	vector<int> *highest_prob_pos_1 = new vector<int>();

		// number of loops to attempt counter
		int number_of_possible_ptms = 1;
		//cout << "Searching combinations of:" << "\n";

		//--------------------------------------------------------------------------------------------------//

		//––––––––––––––––––––––––––//
		//––––––––| 1 PTM |–––––––––//
		//––––––––––––––––––––––––––//

		//    cout << "1 PTM...";

		vector<int> highest_prob_pos_1 (1,list_end);


		int combination_1 = 0;

		for (int i = 0; i < length; i++) {

			if ((roundf(masses[i]*10000)/10000 <= mass_shift_upper) &&
                (roundf(masses[i]*10000)/10000 >= mass_shift_lower)) {
				combination_1++;
				
				if (prob[i] >= prob[highest_prob_pos_1[0]]) {
					highest_prob_pos_1[0] = i;

				}
			}
		}
		//cout << highest_prob_pos_1[0] <<endl;

		float prob_average_1 = 0.0;

		if (combination_1 > 0) {

			prob_average_1 = prob[highest_prob_pos_1[0]];
			//cout << "\n" <<  "The mass shift can be explained by 1 PTM  - ";
			//cout << "with score: " << prob_average_1 << "\n";
			highest_prob_pos[0] = highest_prob_pos_1[0];
			combination = combination_1;
			more_than_3 = 1;

		} else { number_of_possible_ptms++; } 




		//--------------------------------------------------------------------------------------------------//

		//––––––––––––––––––––––––––//
		//––––––––| 2 PTMs |––––––––//
		//––––––––––––––––––––––––––//

//    	cout << "2 PTMs...";
		vector<int> highest_prob_pos_2 (2,list_end);
		int combination_2 = 0;


		for (int j = 0; j < length; j++) {

			int k = j;

			while (masses[j] + masses[k] <= mass_shift_upper) {

				float sum = roundf((masses[j] + masses[k])*10000)/10000;
                
                    if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

					// masses that cancel each other
					if (masses[j] + masses[k] == 0) {
						break;
					}

//                 cout << masses[j] << " " << masses[k] << "\n";
					if (prob[j] + prob[k] >= prob[highest_prob_pos_2[0]] + prob[highest_prob_pos_2[1]]) {
						highest_prob_pos_2[0] = j; highest_prob_pos_2[1] = k;
						//cout << highest_prob_pos_2[0] << " and " << highest_prob_pos_2[1] << endl;
					}
					combination_2++;
				}
				k++;
			}
		}

		float prob_average_2 = 0.0;

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

//     cout << "3 PTMs...";

		vector<int> highest_prob_pos_3 (3,list_end);
		int combination_3 = 0;


		// increasing loop
		if (mass_shift < 160) {

			for (int i = 0; i < length; i++) {

				for (int j = i; j < length; j++) {

					int k = j;


					while (masses[i] + masses[j] + masses[k] <= mass_shift_upper && k <= length) {
						float sum = roundf((masses[i] + masses[j] + masses[k])*10000)/10000;

                        if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

							// masses that cancel each other
							if ((masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
								break;
							}


							if (prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_3[0]] + prob[highest_prob_pos_3[1]] + prob[highest_prob_pos_3[2]]) {
								highest_prob_pos_3[0] = i; highest_prob_pos_3[1] = j; highest_prob_pos_3[2] = k;
								//cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";
							}
							combination_3++;
						}
						k++;
					}
				}
			}

		} else {

			// decreasing loop
			for (int i = length; i >= 0; i--) {

				for (int j = length; j >= i; j--) {

					int k = length;
					while (masses[i] + masses[j] + masses[k] >= mass_shift_lower && k >= 0) {

						float sum = roundf((masses[i] + masses[j] + masses[k])*10000)/10000;
//                    cout << k << "\n";
//                    cout << sum << endl;
                        if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

							// masses that cancel each other
							if ((masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
								break;
							}
                        

							if (prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_3[0]] + prob[highest_prob_pos_3[1]] + prob[highest_prob_pos_3[2]] ) {
								highest_prob_pos_3[0] = i; highest_prob_pos_3[1] = j; highest_prob_pos_3[2] = k;
								//cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";
							}
							combination_3++;
						}
						k--;
					}
				}
			}
		}


		float prob_average_3 = 0.0;

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

			//		break;

		} else { number_of_possible_ptms++; }




		// highest values of all probabilities till 3
		vector<float> prob_vector_1 = { prob_average_1, prob_average_2, prob_average_3};
		float prob_vector_max_1 = *max_element(prob_vector_1.begin(), prob_vector_1.end());
		


		//--------------------------------------------------------------------------------------------------//

		//––––––––––––––––––––––––––//
		//––––––––| 4 PTMs |––––––––//
		//––––––––––––––––––––––––––//

//    cout << "4 PTMs...";
		// zero prob event
		vector<int> highest_prob_pos_4 (4, nat_end);
		int combination_4 = 0;


		//    cout <<nat_prob_length<< endl;

		for (int x = 0; x < nat_prob_length; x++) {

			for (int i = x; i < nat_prob_length; i++) {

				for (int j = i; j < nat_prob_length; j++) {

					int k = j;

					while (nat_prob[x] + nat_prob[i] + nat_prob[j] + nat_prob[k] > (prob_vector_max_1*4) && k <= nat_list) {
                        
                        float sum = roundf((nat_masses[x] + nat_masses[i] + nat_masses[j] + nat_masses[k])*10000)/10000;
                        
                        if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

							// masses that cancel each other
							if ((nat_masses[x] + nat_masses[i] == 0) || (nat_masses[x] + nat_masses[j] == 0) || (nat_masses[x] + nat_masses[k] == 0) ||
							    (nat_masses[i] + nat_masses[j] == 0) || (nat_masses[i] + nat_masses[k] == 0) || (nat_masses[j] + nat_masses[k] == 0)) {
								break;
							}


							float prob_sum = nat_prob[x] + nat_prob[i] + nat_prob[j] + nat_prob[k];

							if (prob_sum >= nat_prob[highest_prob_pos_4[0]] + nat_prob[highest_prob_pos_4[1]] + nat_prob[highest_prob_pos_4[2]] + nat_prob[highest_prob_pos_4[3]]) {
								highest_prob_pos_4[0] = x; highest_prob_pos_4[1] = i; highest_prob_pos_4[2] = j; highest_prob_pos_4[3] = k;
							//cout << x << " " << i << " " << j << " " << k << endl;
							}
							combination_4++;
						}
						k++;
					}
				}
			}
		}

		//--------------------------------------------------------------------------------------------------/n";



		float prob_average_4 = 0.0;

		if (combination_4 > 0) {

			prob_average_4 = (nat_prob[highest_prob_pos_4[0]] + nat_prob[highest_prob_pos_4[1]] + nat_prob[highest_prob_pos_4[2]] + nat_prob[highest_prob_pos_4[3]]) / 4;

			//cout << "\n" << "The mass shift can be explained by 4 PTMs - ";
			//cout << "with score: " << prob_average_4 << "\n";

			if ( prob_average_4 > prob_average_3 && prob_average_4 > prob_average_2 && prob_average_4 > prob_average_1) {
				//cout << "  >>> 4 PTMs have a higher weighting  <<< " << "\n" << "\n";
				highest_prob_pos[0] = highest_prob_pos_4[0]; highest_prob_pos[1] = highest_prob_pos_4[1];
				highest_prob_pos[2] = highest_prob_pos_4[2]; highest_prob_pos[3] = highest_prob_pos_4[3];
				combination = combination_4;
				more_than_3 = 4;
			}

		} else { number_of_possible_ptms++; }






		// highest values of all probabilities
		vector<float> prob_vector_2 = {prob_average_1, prob_average_2, prob_average_3, prob_average_4};
		float prob_vector_max_2 = *max_element(prob_vector_2.begin(), prob_vector_2.end());


		//--------------------------------------------------------------------------------------------------//

		//––––––––––––––––––––––––––//
		//––––––––| 5 PTMs |––––––––//
		//––––––––––––––––––––––––––//

//    cout << "5 PTMs...";
		// zero prob event
		vector<int> highest_prob_pos_5 = {nat_end, nat_end, nat_end, nat_end, nat_end};
		int combination_5 = 0;


		for (int y = 0; y < nat_prob_length; y++) {

			for (int x = y; x < nat_prob_length; x++) {

				for (int i = x; i < nat_prob_length; i++) {

					for (int j = i; j < nat_prob_length; j++) {

						int k = j;

						while (nat_prob[y] + nat_prob[x] + nat_prob[i] + nat_prob[j] + nat_prob[k] > (prob_vector_max_2*5) && k <= nat_list) {

                            float sum = roundf((nat_masses[y] + nat_masses[x] + nat_masses[i] + nat_masses[j] + nat_masses[k])*10000)/10000;
                            
                            if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

								// masses that cancel each other
								if ((nat_masses[y] + nat_masses[x] == 0) || (nat_masses[y] + nat_masses[i] == 0) || (nat_masses[y] + nat_masses[j] == 0) ||
								    (nat_masses[y] + nat_masses[k] == 0) || (nat_masses[x] + nat_masses[i] == 0) || (nat_masses[x] + nat_masses[j] == 0) ||
								    (nat_masses[x] + nat_masses[k] == 0) || (nat_masses[i] + nat_masses[j] == 0) || (nat_masses[i] + nat_masses[k] == 0) ||
								    (nat_masses[j] + nat_masses[k] == 0))  { break; }


								float prob_sum = nat_prob[y] + nat_prob[x] + nat_prob[i] + nat_prob[j] + nat_prob[k];

								if (prob_sum >= nat_prob[highest_prob_pos_5[0]] + nat_prob[highest_prob_pos_5[1]] + nat_prob[highest_prob_pos_5[2]] + nat_prob[highest_prob_pos_5[3]] + nat_prob[highest_prob_pos_5[4]]) {
									highest_prob_pos_5[0] = y; highest_prob_pos_5[1] = x; highest_prob_pos_5[2] = i; highest_prob_pos_5[3] = j; highest_prob_pos_5[4] = k;
									//cout << y << " " << x << " " << i << " " << j << " " << k << endl;
								}
								combination_5++;
							}
							k++;
						}
					}
				}
			}
		}

		//--------------------------------------------------------------------------------------------------//

		float prob_average_5 = 0.0;

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

		//------------------------------------------------------------------------------//


		//––––––––––––––––––––––––––//
		//––––––| 4 PTMs alt |––––––//
		//––––––––––––––––––––––––––//

		if (combination == 0 && number_of_possible_ptms == 6) {

//        cout << "4 PTMs (Alternative)...";

			number_of_possible_ptms = 4;


			// cout << "4 PTMs...";
			vector<int> highest_prob_pos_4 (4,list_end);
			int combination_4 = 0;


			// increasing loop
			if (mass_shift < 160) {

				for (int y = 0; y < length; y++) {

					for (int i = y; i < length; i++) {

						for (int j = i; j < length; j++) {

							int k = j;


							while (masses[y] + masses[i] + masses[j] + masses[k] <= mass_shift_upper && k <= length) {
								float sum = roundf((masses[y] + masses[i] + masses[j] + masses[k])*10000)/10000;

                                if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

									// masses that cancel each other
									if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
									    (masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
										break;
									}


									if (prob[y] + prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_4[0]] + prob[highest_prob_pos_4[1]] + prob[highest_prob_pos_4[2]] + prob[highest_prob_pos_4[3]]) {
										highest_prob_pos_4[0] = y; highest_prob_pos_4[1] = i; highest_prob_pos_4[2] = j; highest_prob_pos_4[3] = k;
										//                            cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";
									}
									combination_4++;
								}
								k++;
							}
						}
					}
				}

			} else {

				// decreasing loop
				for (int y = length; y > 0; y--) {

					for (int i = length; i >= y; i--) {

						for (int j = length; j >= i; j--) {

							int k = length;
							while (masses[y] + masses[i] + masses[j] + masses[k] >= mass_shift_lower && k >= 0) {

								float sum = roundf((masses[y] + masses[i] + masses[j] + masses[k])*10000)/10000;
								//                    cout << k << "\n";
								//                    cout << sum << endl;
                                if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

									// masses that cancel each other
									if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
									    (masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
										break;
									}
									//                        cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";

									if (prob[y] + prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_4[0]] + prob[highest_prob_pos_4[1]] + prob[highest_prob_pos_4[2]] + prob[highest_prob_pos_4[3]]) {
										highest_prob_pos_4[0] = y; highest_prob_pos_4[1] = i; highest_prob_pos_4[2] = j; highest_prob_pos_4[3] = k;
									}
									combination_4++;
								}
								k--;
							}
						}
					}
				}
			}


			float prob_average_4 = 0.0;

			if (combination_4 > 0) {

				prob_average_4 = (prob[highest_prob_pos_4[0]] + prob[highest_prob_pos_4[1]] + prob[highest_prob_pos_4[2]] + prob[highest_prob_pos_4[3]]) / 4;

				cout << "\n" << "The mass shift can be explained by 4 PTMs - ";
				cout << "with score: " << prob_average_4 << "\n";

				if ( prob_average_4 > prob_average_3 && prob_average_4 > prob_average_2 && prob_average_4 > prob_average_1) {
					cout << "  >>> 4 PTMs have a higher weighting  <<< " << "\n" << "\n";
					highest_prob_pos[0] = highest_prob_pos_4[0]; highest_prob_pos[1] = highest_prob_pos_4[1];
					highest_prob_pos[2] = highest_prob_pos_4[2]; highest_prob_pos[3] = highest_prob_pos_4[3];
					combination = combination_4;
					more_than_3 = 6;
				}

				//		break;

			} else {
				number_of_possible_ptms = 6;
			}
		}


	vector<int> highest_prob_pos_5 (5,list_end);
			int combination_5 = 0;


			// increasing loop
			if (mass_shift < 160) {

				for (int y = 0; y < length; y++) {

					for (int i = y; i < length; i++) {

						
						for (int j = i; j < length; j++) {

							int k = j;


							while (masses[y] + masses[i] + masses[j] + masses[k] <= mass_shift_upper && k <= length) {
								float sum = roundf((masses[y] + masses[i] + masses[j] + masses[k])*10000)/10000;

                                if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

									// masses that cancel each other
									if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
									    (masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
										break;
									}


									if (prob[y] + prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_4[0]] + prob[highest_prob_pos_4[1]] + prob[highest_prob_pos_4[2]] + prob[highest_prob_pos_4[3]]) {
										highest_prob_pos_4[0] = y; highest_prob_pos_4[1] = i; highest_prob_pos_4[2] = j; highest_prob_pos_4[3] = k;
										//                            cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";
									}
									combination_4++;
								}
								k++;
							}
						}
					}
				}

			} else {

				// decreasing loop
				for (int y = length; y > 0; y--) {

					for (int i = length; i >= y; i--) {

						for (int j = length; j >= i; j--) {

							int k = length;
							while (masses[y] + masses[i] + masses[j] + masses[k] >= mass_shift_lower && k >= 0) {

								float sum = roundf((masses[y] + masses[i] + masses[j] + masses[k])*10000)/10000;
								//                    cout << k << "\n";
								//                    cout << sum << endl;
                                if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

									// masses that cancel each other
									if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
									    (masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
										break;
									}
									//                        cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";

									if (prob[y] + prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_4[0]] + prob[highest_prob_pos_4[1]] + prob[highest_prob_pos_4[2]] + prob[highest_prob_pos_4[3]]) {
										highest_prob_pos_4[0] = y; highest_prob_pos_4[1] = i; highest_prob_pos_4[2] = j; highest_prob_pos_4[3] = k;
									}
									combination_4++;
								}
								k--;
							}
						}
					}
				}
			}


			float prob_average_4 = 0.0;

			if (combination_4 > 0) {

				prob_average_4 = (prob[highest_prob_pos_4[0]] + prob[highest_prob_pos_4[1]] + prob[highest_prob_pos_4[2]] + prob[highest_prob_pos_4[3]]) / 4;

				cout << "\n" << "The mass shift can be explained by 4 PTMs - ";
				cout << "with score: " << prob_average_4 << "\n";

				if ( prob_average_4 > prob_average_3 && prob_average_4 > prob_average_2 && prob_average_4 > prob_average_1) {
					cout << "  >>> 4 PTMs have a higher weighting  <<< " << "\n" << "\n";
					highest_prob_pos[0] = highest_prob_pos_4[0]; highest_prob_pos[1] = highest_prob_pos_4[1];
					highest_prob_pos[2] = highest_prob_pos_4[2]; highest_prob_pos[3] = highest_prob_pos_4[3];
					combination = combination_4;
					more_than_3 = 6;
				}

				//		break;

			} else {
				number_of_possible_ptms = 6;
			}
		}
/*
		//=======================================================================================================//

		// No results
		if (combination == 0 && number_of_possible_ptms == 6) {
			cout << "\n" << "The mass shift of " <<  peptide << " can NOT be explained till 5 PTMs" << "\n" << "\n";

			ofstream output_file;
			output_file.open(output_file_path, ios::out | ios::app);

			output_file << "peptide" << "\t" << "mass_shift" << "\t" << "1" << "\n";
			output_file << peptide << "\t" << mass_shift << "\t" << "NA" << "\n";
			output_file.close();


			//// skip to next peptide
			continue;
		}

		//=======================================================================================================//


		//================================| Output |================================//

		// timer
		auto end = chrono::steady_clock::now();
		auto diff = end - start;
		string elapsed_time;

		if (chrono::duration <double, milli> (diff).count() < 1) {
			elapsed_time = to_string(chrono::duration <double, milli> (diff).count()) + " ms";
		} else {
			elapsed_time = to_string(chrono::duration <double, milli> (diff).count() / 1000) + " s";
		}

		//// print output combination stats
		cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << "\n\n";
		cout << "peptide:  " << peptide << "\n\n";
		cout << "PTM list:   " << length << ";  number of PTMs: " << number_of_possible_ptms << "\n";
		cout << "mass shift: " << mass_shift << " ± " << mass_error << ";  combinations: " << combination << "\n" << "\n";
		cout << "=======================================================================" << "\n" << "\n";


//    cout << prob_average_3 << " " <<prob_average_4 << "\n";
		vector<string> out_name;

		//// 1-3 PTMs
		if (more_than_3 <= 3) {
			//// print PTM names from file
			fstream inputPTMnames;
			inputPTMnames.open(dir_path + "new_names.txt");
			string sline;


			cout << "PTM combination: ";
			for (int i = 0; i < 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					GotoLine(inputPTMnames, 1 + highest_prob_pos[i]);
					string line;
					inputPTMnames >> line;
					cout << line << "\t";
					out_name.push_back(line);
				}
			}

			cout << "\n\nmass shift sum: ";
			for (int i = 0; i < 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					cout << masses[highest_prob_pos[i]] << "\t";
				}
			}


			// print positions
			cout << "\n\npositions: ";

			for (int i = 0; i < 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					cout << 1 + highest_prob_pos[i] << " ";
				}
			}
		}

		//// 4-5 PTMs
		else if (more_than_3 <= 5) {

			cout << "\nPTM combination: ";
			for (int i = 0; i <= 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					cout << nat_names[highest_prob_pos[i]] << "\t";
				}
			}


			cout << "\n\nmass shift sum:  ";
			for (int i = 0; i <= 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					cout << nat_masses[highest_prob_pos[i]] << "\t";
				}
			}


			cout << "\npositions: ";
			for (int i = 0; i <= 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					cout << 1 + highest_prob_pos[i] << " ";
				}
			}
		}

		else if (more_than_3 == 6) {

			//// print PTM names from file
			fstream inputPTMnames;
			inputPTMnames.open(dir_path + "new_names.txt");
			string sline;


			cout << "PTM combination: ";
			for (int i = 0; i <= 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					GotoLine(inputPTMnames, 1 + highest_prob_pos[i]);
					string line;
					inputPTMnames >> line;
					cout << line << "\t";
					out_name.push_back(line);
				}
			}

			cout << "\n\nmass shift sum: ";
			for (int i = 0; i <= 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					cout << masses[highest_prob_pos[i]] << "\t";
				}
			}


			// print positions
			cout << "\n\npositions: ";

			for (int i = 0; i <= 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					cout << 1 + highest_prob_pos[i] << " ";
				}
			}
		}


		// print run-time
		cout << "\n\n" << "execution time = " << elapsed_time << "\n\n";
		cout << "=======================================================================" << "\n";



		//--------------------------------------------------------------------------------------------------/


		//================================| Write Output |================================//

		//// write output into file
		ofstream output_file;
		output_file.open(output_file_path, ios::out | ios::app);


		//// 1-3 PTMs
		if ((more_than_3 <= 3) || (more_than_3 == 6)) {


			output_file << "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t";
			for (int i = 1; i <= out_name.size(); i++) {
				output_file << i << "\t";
			}


			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
			for (int i = 0; i < out_name.size(); i++) {
				output_file << out_name[i] << "\t";
			}


			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
			for (int i = 0; i < out_name.size(); i++) {
				if (highest_prob_pos[i] != list_end) {
					output_file << masses[highest_prob_pos[i]] << "\t";
				}
			}

			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
			for (int i = 0; i < out_name.size(); i++) {
				if (highest_prob_pos[i] != list_end) {
					output_file << prob[highest_prob_pos[i]] << "\t";
				}
			}

			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
			for (int i = 0; i < out_name.size(); i++) {
				if (highest_prob_pos[i] != list_end) {
					output_file << 1 + highest_prob_pos[i] << "\t";
				}
			}



			//// 4-5 PTMs
		} else if (more_than_3 <= 5) {

			if (more_than_3 == 4) {
				output_file <<  "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t" << "1" << "\t" << "2" << "\t" << "3" << "\t" << "4" << "\t";
			} else if (more_than_3 == 5) { output_file << "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t" << "1" << "\t" << "2" << "\t" << "3" << "\t" << "4" << "\t" << "5" << "\t"; }


			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
			for (int i = 0; i <= 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					output_file << nat_names[highest_prob_pos[i]] << "\t";
				}
			}

			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
			for (int i = 0; i <= 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					output_file << nat_masses[highest_prob_pos[i]] << "\t";
				}
			}

			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
			for (int i = 0; i <= 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					output_file << nat_prob[highest_prob_pos[i]] << "\t";
				}
			}

			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
			for (int i = 0; i <= 4; i++) {
				if (highest_prob_pos[i] != list_end) {
					output_file << 1 + highest_prob_pos[i] << "\t";
				}
			}
		}

		output_file << "\n";
		output_file.close();

		cout << "\n\n";

		//--------------------------------------------------------------------------------------------------/


		//================================| Alt Output |================================//

		//// write alternative outputs into files

		string alt_output_file_path = dir_path + "results/" + sample + "/" + peptide + "_alternative.txt";

		ofstream alt_output_file;
		alt_output_file.open(alt_output_file_path, ios::out | ios::app);




		if (combination_1 > 0 && prob_average_1 > 0) {

			alt_output_file << "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t";

			for (int i = 1; i < 1; i++) {
				alt_output_file << i << "\t";
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
			for (int i = 0; i < 1; i++) {
//                alt_output_file << names[highest_prob_pos_1[i]] << "\t";
				alt_output_file << "\t";
			}


			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
			for (int i = 0; i < 1; i++) {
				if (highest_prob_pos_1[i] != list_end) {
					alt_output_file << masses[highest_prob_pos_1[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
			for (int i = 0; i < 1; i++) {
				if (highest_prob_pos_1[i] != list_end) {
					alt_output_file << prob[highest_prob_pos_1[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
			for (int i = 0; i < 1; i++) {
				if (highest_prob_pos_1[i] != list_end) {
					alt_output_file << 1 + highest_prob_pos_1[i] << "\t";
				}
			}


		} else if (combination_2 > 0 && prob_average_2 > 0) {

			alt_output_file << "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t";

			for (int i = 1; i < 2; i++) {
				alt_output_file << i << "\t";
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
			for (int i = 0; i < 2; i++) {
//                alt_output_file << names[highest_prob_pos_1[i]] << "\t";
				alt_output_file << "\t";
			}


			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
			for (int i = 0; i < 2; i++) {
				if (highest_prob_pos_2[i] != list_end) {
					alt_output_file << masses[highest_prob_pos_2[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
			for (int i = 0; i < 2; i++) {
				if (highest_prob_pos_2[i] != list_end) {
					alt_output_file << prob[highest_prob_pos_2[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
			for (int i = 0; i < 2; i++) {
				if (highest_prob_pos_2[i] != list_end) {
					alt_output_file << 1 + highest_prob_pos_2[i] << "\t";
				}
			}

		} else if (combination_3 > 0 && prob_average_3 > 0) {

			alt_output_file << "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t";

			for (int i = 1; i < 3; i++) {
				alt_output_file << i << "\t";
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
			for (int i = 0; i < 3; i++) {
				// alt_output_file << names[highest_prob_pos_1[i]] << "\t";
				alt_output_file << "\t";
			}


			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
			for (int i = 0; i < 3; i++) {
				if (highest_prob_pos_3[i] != list_end) {
					alt_output_file << masses[highest_prob_pos_3[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
			for (int i = 0; i < 3; i++) {
				if (highest_prob_pos_3[i] != list_end) {
					alt_output_file << prob[highest_prob_pos_3[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
			for (int i = 0; i < 3; i++) {
				if (highest_prob_pos_3[i] != list_end) {
					alt_output_file << 1 + highest_prob_pos_3[i] << "\t";
				}
			}
		}

		else if ( more_than_3 == 6 && combination_4 > 0 && prob_average_4 > 0) {

			alt_output_file << "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t";

			for (int i = 1; i < 4; i++) {
				alt_output_file << i << "\t";
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
			for (int i = 0; i < 4; i++) {
				// alt_output_file << names[highest_prob_pos_1[i]] << "\t";
				alt_output_file << "\t";
			}


			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
			for (int i = 0; i < 4; i++) {
				if (highest_prob_pos_4[i] != list_end) {
					alt_output_file << masses[highest_prob_pos_4[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
			for (int i = 0; i < 4; i++) {
				if (highest_prob_pos_4[i] != list_end) {
					alt_output_file << prob[highest_prob_pos_4[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
			for (int i = 0; i < 4; i++) {
				if (highest_prob_pos_4[i] != list_end) {
					alt_output_file << 1 + highest_prob_pos_4[i] << "\t";
				}
			}
		}





		if (more_than_3 == 4 && combination_4 > 0) {
			alt_output_file <<  "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t" << "1" << "\t" << "2" << "\t" << "3" << "\t" << "4" << "\t";

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
			for (int i = 0; i < 4; i++) {
				if (highest_prob_pos_4[i] != list_end) {
					alt_output_file << nat_names[highest_prob_pos_4[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
			for (int i = 0; i < 4; i++) {
				if (highest_prob_pos_4[i] != list_end) {
					alt_output_file << nat_masses[highest_prob_pos_4[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
			for (int i = 0; i < 4; i++) {
				if (highest_prob_pos_4[i] != list_end) {
					alt_output_file << nat_prob[highest_prob_pos_4[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
			for (int i = 0; i < 4; i++) {
				if (highest_prob_pos_4[i] != list_end) {
					alt_output_file << 1 + highest_prob_pos_4[i] << "\t";
				}
			}

		} else if (more_than_3 == 5 && combination_5 > 0) {
			alt_output_file << "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t" << "1" << "\t" << "2" << "\t" << "3" << "\t" << "4" << "\t" << "5" << "\t";

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";

			for (int i = 0; i < 5; i++) {
				if (highest_prob_pos_5[i] != list_end) {
					alt_output_file << nat_names[highest_prob_pos_5[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
			for (int i = 0; i < 5; i++) {
				if (highest_prob_pos_5[i] != list_end) {
					alt_output_file << nat_masses[highest_prob_pos_5[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
			for (int i = 0; i < 5; i++) {
				if (highest_prob_pos_5[i] != list_end) {
					alt_output_file << nat_prob[highest_prob_pos_5[i]] << "\t";
				}
			}

			alt_output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
			for (int i = 0; i < 5; i++) {
				if (highest_prob_pos_5[i] != list_end) {
					alt_output_file << 1 + highest_prob_pos_5[i] << "\t";
				}
			}
		}


		alt_output_file << "\n";
		alt_output_file.close();


*/
	}
		auto end = chrono::steady_clock::now();
		auto diff = end - start;
		string elapsed_time;

		if (chrono::duration <double, milli> (diff).count() < 1) {
			elapsed_time = to_string(chrono::duration <double, milli> (diff).count()) + " ms";
		} else {
			elapsed_time = to_string(chrono::duration <double, milli> (diff).count() / 1000) + " s";
		}
	cout << elapsed_time << endl;
	return 0;

}



//a->push_back(3);
//cout << "here " << a->at(0);


