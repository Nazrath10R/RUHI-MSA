
//  main.cpp
//  ptm
//
//  Created by Nazrath Nawaz on 28/08/2019.
//  Developed by Pedro Cardoso on 09/06/2020.
//  Copyright © 2019 Nazrath Nawaz. All rights reserved.
//
//uiuiuiu
#include <iostream>
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


//=======================================================================================================//


using namespace std;


// function to extract line of file based on position specified
void GotoLine(fstream& file, unsigned int num) {
	file.seekg(ios::beg);
	for (int i = 0; i < num - 1; ++i) {
		file.ignore(numeric_limits<streamsize>::max(), '\n');
	}
}

//fuction that finds all elements of a list that are within a certain range
vector<int> binarySearch( float masses[], const int &array_Size, const float &mass_shift_lower, const float &mass_shift_upper) {

	size_t first = 0;                 //first array element
	size_t last = array_Size - 1;            //last array element
	size_t middle {};                        //mid point of search
	vector<int> position {};                 //position of search value

bool found = false;            //flag

	while ((first<=last) && (found==false)) {
		middle = (first + last) / 2; //this finds the mid point
		if ((roundf(masses[middle]*10000)/10000 <= mass_shift_upper) &&
                (roundf(masses[middle]*10000)/10000 >= mass_shift_lower)) {
			position.push_back(middle); //we might have more than one match! need to change code for that 
			int set_middle = middle; // to use on the way down the array
			middle = middle + 1; //moving up the array
			while (roundf(masses[middle]*10000)/10000 <= mass_shift_upper) {
				position.push_back(middle);
				middle++;
			}
			set_middle = set_middle - 1; //moving down the array
			while (roundf(masses[set_middle]*10000)/10000 >= mass_shift_lower) {
				position.push_back(set_middle);
				set_middle--;
			}
			
			found = true;
		}
		else if (masses[middle] > mass_shift_upper) { // if it's in the lower half
			last = middle - 1;
		}
		else {
			first = middle + 1;                 //if it's in the upper half
		}
	}
	return position;
}

//using Binary Search principle to obtain specific value
int binarySearchSmallerElement_decreasing(float masses[], const float &mass_shift, const int &size_of_array){
    int l = 0;
	int r = size_of_array -1;
    int resultantIndex = size_of_array;
    while(l<=r) {
        int mid = l + (r-l) / 2;
        if(masses[mid] >= mass_shift){
            l = mid +1;
        }
        else{
            r = mid -1;
            resultantIndex = mid;
        }
    }
	if (resultantIndex == size_of_array)
		return (resultantIndex -1);
	else
		return (resultantIndex);
}

//using Binary Search principle to obtain specific value
int SpecificBinarySearch( float masses[], const int &array_Size, const float &number) {

	size_t first = 0;                 //first array element
	size_t last = array_Size - 1;            //last array element
	size_t middle {};                        //mid point of search

	while (first<=last)  {
		middle = (first + last) / 2; //this finds the mid point
		if ((roundf(masses[middle]*10000)/10000 == number)) {
			return middle;
			
		}
		else if ((roundf(masses[middle]*10000)/10000 > number)) { // if it's in the lower half
			last = middle - 1;
		}
		else {
			first = middle + 1;                 //if it's in the upper half
		}
	}
	return -1; //not found
}

//using Binary Search principle to obtain fisrt value of the list that is bigger than a specific value (mass_shift)
int binarySearchLargerElement( float masses[], const float &mass_shift, const int &size_of_array) {
    int l = 0;
	int r = size_of_array -1;
    int resultantIndex = size_of_array;
    while(l<=r) {
        int mid = l + (r-l) / 2;
        if(masses[mid] <= mass_shift){
            l = mid +1;
        }
        else{
            r = mid -1;
            resultantIndex = mid;
        }
    }
    if (resultantIndex == size_of_array)
		return (resultantIndex -1);
	else
		return (resultantIndex);
}

std::tuple<vector<int>*, vector<int>*, int*> finding1PTM ( float masses [],float prob [], const int &PTM_list, const float &mass_shift_lower, const float &mass_shift_upper) {
		
		vector<int> highest_prob_pos_1 {PTM_list-1};
		
		
		//Storing all possible PTMs on position_combination1
		vector<int> position_combination1 =  binarySearch(masses, PTM_list, mass_shift_lower, mass_shift_upper); //all PTMs that are withing our tolerance
		int pos_length = position_combination1.size(); 		
		if (pos_length > 1) //checking which of the PTMs is more probable and storing
			for ( size_t i {0}; i < pos_length; i++) { //
					//cout << position_combination1[i] << " is the position! The number is " << masses[position_combination1[i]] 
					//cout << " that is bigger than " << mass_shift_lower << " and lower than " << mass_shift_upper << endl;
					if (prob[position_combination1[i]] >= prob[highest_prob_pos_1[0]])
						highest_prob_pos_1[0] = position_combination1[i];
			}
		else if (pos_length ==1) {
			highest_prob_pos_1[0] = position_combination1[0];
		}
		else {
			//cout << "\n" <<  "The mass shift can not be explained by 1 PTM " << endl;
		}
		return {&highest_prob_pos_1, &position_combination1, &pos_length};
}


std::tuple<vector<int>, vector<int>, int> findingCombinationOf2PTM ( float masses [],float prob [], const int &PTM_list, const float &mass_shift_lower, const float &mass_shift_upper) {
	vector<int> position_combination2 {}; // here is were we store all matched positions
	int list_end = PTM_list-1;
	vector<int> highest_prob_pos_2 {list_end, list_end};
	int combination_2 = 0;

	if (mass_shift_upper > 0) {
		float maximum = mass_shift_upper + abs(masses[0]); //find the maximum mass shift that we need to look.
		int limit {};
		limit = binarySearchLargerElement(masses, maximum, PTM_list); //find position of first larger number of the maximum
		
		float minimum = mass_shift_upper/2; 
		int minimum_limit {};
		minimum_limit = binarySearchLargerElement( masses, minimum, PTM_list) - 1;
		//cout << minimum << " " << minimum_limit << "  " << masses[minimum_limit] << endl;
				
		for (int j = 0; j < limit; j++) {

			int k = minimum_limit;

			while (masses[j] + masses[k] <= mass_shift_upper) {

				float sum = roundf((masses[j] + masses[k])*10000)/10000;
                
					if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

					// masses that cancel each other
					if (masses[j] + masses[k] == 0) {
						k++;
						continue;//does break statement make sense here?? think not because we can have more matches
					}
					else {
						position_combination2.push_back(j);
						position_combination2.push_back(k);
//                 cout << masses[j] << " " << masses[k] << "\n";
						if (prob[j] + prob[k] >= prob[highest_prob_pos_2[0]] + prob[highest_prob_pos_2[1]]) {
							highest_prob_pos_2[0] = j; highest_prob_pos_2[1] = k;
							//cout << " \n" << highest_prob_pos_2[0] << " and " << highest_prob_pos_2[1] << endl;
						}
					combination_2++;
					}
				}
				if (k==j)
					k = j + 1; // dont repeat combinations
				else
					k++; 
			}
		}
	}
	else {
		int maximum = abs(masses[0]);//the maximum mass shift is the absolute value of the first value 
		int limit {};
		limit = binarySearchLargerElement(masses, maximum, PTM_list); //find position of first larger number of the maximum
		
		float minimum {0}; // because adding two positive number will never be a negative number
		int minimum_limit {};
		minimum_limit = binarySearchLargerElement( masses, minimum, PTM_list); // no need of -1 becasue here we are looking for the first postivie number. mass_shift will not be 0.
		//cout << minimum << " " << minimum_limit << "  " << masses[minimum_limit] << endl;
		//hcek this values above...		
		for (int j = 0; j < limit; j++) {

			int k = minimum_limit;

			while (masses[j] + masses[k] <= mass_shift_upper) {

				float sum = roundf((masses[j] + masses[k])*10000)/10000;
                
					if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {
					
					// masses that cancel each other
						if (masses[j] + masses[k] == 0) {
							k++;
							continue; //does break statement make sense here?? think not because we can have more matches
						}
						else {
							position_combination2.push_back(j); // can use insert below to keep the most probable in 0 and 1;
							position_combination2.push_back(k);
//                 cout << masses[j] << " " << masses[k] << "\n";
							if (prob[j] + prob[k] >= prob[highest_prob_pos_2[0]] + prob[highest_prob_pos_2[1]]) {
								highest_prob_pos_2[0] = j; highest_prob_pos_2[1] = k;
								//cout << " \n" << highest_prob_pos_2[0] << " and " << highest_prob_pos_2[1] << endl;
							}
							combination_2++;
						}
					}
					if (k==j)
						k = j+ 1;
					else
						k++;
			}
		}
	}
	return {highest_prob_pos_2, position_combination2, combination_2};
	}

std::tuple<vector<int>*, int*> findingCombinationOf3PTM ( float masses [],float prob [], const int &PTM_list, const float &mass_shift_lower, const float &mass_shift_upper) {
		
		int list_end {PTM_list-1};
		vector<int> highest_prob_pos_3 {list_end, list_end, list_end};
		int combination_3 = 0;

		float maximum = mass_shift_upper + (abs(masses[0])*2);//the maximum mass shift we need to check in the first loop 
		int limit {};
		if (maximum < masses[list_end]) 
			limit = binarySearchLargerElement(masses, maximum, PTM_list);
		else
			limit = list_end;
		
		// increasing loop
		if (mass_shift_upper < 160) {
				for (int i = 0; i < limit; i++) {
					float maximum_1 = -(masses[i]) + abs(masses[0]) + mass_shift_upper;
					int limit_1 = binarySearchLargerElement( masses, maximum_1, PTM_list);
					
					
					for (int j = i; j < limit_1; j++) {
						
						float minimum_1 {};
						minimum_1 = -(masses[i]) - (masses[j]) + mass_shift_lower;
						int minimum_limit_1 = binarySearchLargerElement (masses, minimum_1, PTM_list) -1;
						int k {};
						if (minimum_limit_1 < j)
							k = j;
						else
							k = minimum_limit_1;

						while (masses[i] + masses[j] + masses[k] <= mass_shift_upper && k < PTM_list) {
								float sum = roundf((masses[i] + masses[j] + masses[k])*10000)/10000;

								if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

								// masses that cancel each other
								if ((masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
									k++;
									continue;
								}


								if ((prob[i] + prob[j] + prob[k]) >= (prob[highest_prob_pos_3[0]] + prob[highest_prob_pos_3[1]] + prob[highest_prob_pos_3[2]])) {
									highest_prob_pos_3[0] = i; highest_prob_pos_3[1] = j; highest_prob_pos_3[2] = k;
                            		//cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";
								}
								combination_3++;
								}
								if (k==j)
									k = j + 1;
								k++;
						}
					}
				}

		} else {
			// decreasing loop
			for (int i = limit; i >= 0; i--) {
				float maximum_1 = -(masses[i]) + abs(masses[0]) + mass_shift_upper;
				int limit_1 = binarySearchLargerElement( masses, maximum_1, PTM_list);

				for (int j = limit_1; j >= i; j--) {
					float minimum_1 {};
					minimum_1 = -(masses[i]) - (masses[j]) + mass_shift_lower;
					int minimum_limit_1 = binarySearchLargerElement (masses, minimum_1, PTM_list); //not -1 because in this case we are decreasing in the loop
					int k = minimum_limit_1;

					while (masses[i] + masses[j] + masses[k] >= mass_shift_lower && k >= 0) {

						float sum = roundf((masses[i] + masses[j] + masses[k])*10000)/10000;
//                    cout << k << "\n";
//                    cout << sum << endl;
						if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

							// masses that cancel each other
							if ((masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
								k--;
								continue;
							}
							if (prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_3[0]] + prob[highest_prob_pos_3[1]] + prob[highest_prob_pos_3[2]] ) {
								highest_prob_pos_3[0] = i; highest_prob_pos_3[1] = j; highest_prob_pos_3[2] = k;
								//cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";
							}
							combination_3++;
						}
						if (k==j)
							k = j - 1;
						else
							k--;
					}
				}
			}
		}
		return {&highest_prob_pos_3, &combination_3};
}

std::tuple<vector<int>,vector<int>, int, vector<int>, int> findingCombinationOf2and3PTM ( float masses [],float prob [], const int &PTM_list, const float &mass_shift_lower, const float &mass_shift_upper) {
	
	vector<int> position_combination2 {}; // here is were we store all matched positions
	int list_end = PTM_list-1;

	vector<int> highest_prob_pos_2 {list_end, list_end};
	int combination_2 = 0;
	vector<int> highest_prob_pos_3 {list_end, list_end, list_end};
	int combination_3 = 0;
	
	float list_top_value = masses[list_end]; //maximum value of the list
	float list_lowest_value = masses[0]; //lowest value of the list
	float dif_upper_lower = mass_shift_upper - mass_shift_lower;

	float maximum = mass_shift_upper + (abs(list_lowest_value)*2);//the maximum mass shift we need to check in the first loop 
	int limit {};
	if (maximum < list_top_value) 
		limit = binarySearchLargerElement(masses, maximum, PTM_list);
	else
		limit = list_end;
		
	// increasing loop
		if (mass_shift_upper < 160) {
			for (int i = 0; i < limit; i++) {
				float maximum_1 = -(masses[i]) - (list_lowest_value) + mass_shift_upper;
				int limit_1 = binarySearchLargerElement( masses, maximum_1, PTM_list);
				int j {i};
				float maximum_lowest = -(masses[i]) - (list_top_value) + mass_shift_upper;
				if (maximum_1 < list_lowest_value || maximum_lowest > list_top_value) 
					break;
//				else
//					j = binarySearchLargerElement( masses, maximum_lowest, PTM_list);
//				if (j<i)
//					j=i;
				
				while (masses[i] + masses[j] <= mass_shift_upper && j < limit_1) {
						
					float sum = roundf((masses[i] + masses[j])*10000)/10000;
                
					if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

					// masses that cancel each other
					if (masses[i] + masses[j] == 0) {
						j++;
						continue;//does break statement make sense here?? think not because we can have more matches
					}
					else {
						//position_combination2.push_back(i);
						//position_combination2.push_back(j);
//                 cout << masses[j] << " " << masses[k] << "\n";
						if (prob[i] + prob[j] >= prob[highest_prob_pos_2[0]] + prob[highest_prob_pos_2[1]]) {
							highest_prob_pos_2[0] = i; highest_prob_pos_2[1] = j;
							//cout << masses[i] << " " << masses[j] << "\n";
						}
					combination_2++;
					}
					}
					
					float minimum_1 {};
					minimum_1 = -(masses[i]) - (masses[j]) + mass_shift_lower;
					float minimum_1_highest = minimum_1 + dif_upper_lower;
					int k{}; 
					if (minimum_1_highest < list_lowest_value && minimum_1 > list_top_value) {
						j++; continue; }
					else
						k = binarySearchLargerElement (masses, minimum_1, PTM_list) -1;
					if (k<j)
						k = j;

					while (masses[i] + masses[j] + masses[k] <= mass_shift_upper && k < PTM_list) {
							float sum = roundf((masses[i] + masses[j] + masses[k])*10000)/10000;

							if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

							// masses that cancel each other
							if ((masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
								k++;
								continue;
							}
							if ((prob[i] + prob[j] + prob[k]) >= (prob[highest_prob_pos_3[0]] + prob[highest_prob_pos_3[1]] + prob[highest_prob_pos_3[2]])) {
								highest_prob_pos_3[0] = i; highest_prob_pos_3[1] = j; highest_prob_pos_3[2] = k;
								//cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";
							}
							combination_3++;
							}
								k++;
					}
						j++;
				}
			}
	} else {
		// decreasing loop
		for (int i = limit; i >= 0; i--) {
			float maximum_1 = -(masses[i]) - (list_lowest_value) + mass_shift_upper;
				int limit_1 = binarySearchLargerElement( masses, maximum_1, PTM_list);
				int j {};
				float maximum_lowest = -(masses[i]) - (list_top_value) + mass_shift_upper;
				if (maximum_1 < list_lowest_value || maximum_lowest > list_top_value) {
				i--; continue; }
				else
					j = limit_1;
				if (j > i)
					j = i;
				

			while (masses[i] + masses[j] >= mass_shift_lower && j >= 0) {
				float sum = roundf((masses[i] + masses[j])*10000)/10000;
                
					if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {
					
					// masses that cancel each other
						if (masses[i] + masses[j] == 0) {
							j--;
							continue; //does break statement make sense here?? think not because we can have more matches
						} 
						else {
							//position_combination2.push_back(i); // can use insert below to keep the most probable in 0 and 1;
							//position_combination2.push_back(j);
							//cout << masses[j] << " " << masses[k] << "\n";
							if (prob[i] + prob[j] >= prob[highest_prob_pos_2[0]] + prob[highest_prob_pos_2[1]]) {
								highest_prob_pos_2[0] = i; highest_prob_pos_2[1] = j;
								//cout << masses[i] << " " << masses[j] << "\n";
							}
							combination_2++;
						}
					}

				float maximum_2 = -(masses[i]) - (masses[j]) + mass_shift_upper;
//				int limit_2 = binarySearchLargerElement( masses, maximum_2, PTM_list);
				int k {};
				float maximum2_lowest = -(masses[i]) - (list_top_value) + mass_shift_upper;
				if (maximum_2 < list_lowest_value || maximum2_lowest > list_top_value) {
				j--; continue; }
				else
					k = limit_1;
				if (k > j)
					k = j;

				while (masses[i] + masses[j] + masses[k] >= mass_shift_lower && k >= 0) {

					float sum = roundf((masses[i] + masses[j] + masses[k])*10000)/10000;

					if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

					// masses that cancel each other
						if ((masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
							k--;
							continue;
						}
						if (prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_3[0]] + prob[highest_prob_pos_3[1]] + prob[highest_prob_pos_3[2]] ) {
							highest_prob_pos_3[0] = i; highest_prob_pos_3[1] = j; highest_prob_pos_3[2] = k;
							//cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";
						}
						combination_3++;
					}
					if (k>j)
						k = j;
					else
						k--;
				}
				if (j>i)
					j = i;
				else
					j--;
			}
			
		}
	}
		return {highest_prob_pos_2, position_combination2, combination_2, highest_prob_pos_3, combination_3};
}

std::tuple<vector<int>, int> findingCombinationOf4NaturalPTM (float nat_prob [], float nat_masses [],float prob [], const int &nat_prob_length, const float &mass_shift_lower, const float &mass_shift_upper, float &sum_of_probs_4) {
	int nat_end = nat_prob_length -1;
	// zero prob event
	vector<int> highest_prob_pos_4 = {nat_end, nat_end, nat_end, nat_end};
	int combination_4 = 0;
	
	if (sum_of_probs_4 < (4*nat_prob[0])) { //check if is possible to obtain a better score 
		for (int x = 0; x < nat_prob_length; x++) {
			float maximum_prob = -(nat_prob[x]) - (nat_prob[x]*2) + sum_of_probs_4;
			int limit_prob {};
			if (maximum_prob > nat_prob[x])
				break; // if maximum is more than nat_prob[0] then it means that sum_of_probs will never be reached
			else if ( maximum_prob < 0)
				limit_prob = nat_prob_length; // if is 
			else {
				limit_prob = binarySearchSmallerElement_decreasing (nat_prob, maximum_prob, nat_end);
					}
			for (int i = x; i < limit_prob; i++) {
					float maximum_1 = -(nat_prob[x]) - (nat_prob[i]) - (nat_prob[i]) + sum_of_probs_4;
					int limit1_prob {};
					if (maximum_1 > nat_prob[i])
						break; // if maximum is more than nat_prob[0] then it means that sum_of_probs will never be reached
					else if ( maximum_1 < 0)
						limit1_prob = nat_prob_length;
					else {
						limit1_prob = binarySearchSmallerElement_decreasing (nat_prob, maximum_1, nat_end);
					}
						
				for (int j = i; j < limit1_prob; j++) {
					float maximum_2 = -(nat_prob[x]) - (nat_prob[i]) - (nat_prob[j]) + sum_of_probs_4;
					int k {j};
					int limit2 {};
					
					if (maximum_2 > nat_prob[j])
						break; // if maximum is more than nat_prob[0] then it means that sum_of_probs will never be reached
					else if ( maximum_2 < 0)
						limit2 = nat_prob_length;
					else {
						limit2 = binarySearchSmallerElement_decreasing (nat_prob, maximum_2, nat_end);
					}
					
					while (nat_prob[x] + nat_prob[i] + nat_prob[j] + nat_prob[k] > sum_of_probs_4 && k < limit2) {
                        
                        float sum = roundf((nat_masses[x] + nat_masses[i] + nat_masses[j] + nat_masses[k])*10000)/10000;

                        if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

							// masses that cancel each other
							if ((nat_masses[x] + nat_masses[i] == 0) || (nat_masses[x] + nat_masses[j] == 0) || (nat_masses[x] + nat_masses[k] == 0) ||
							    (nat_masses[i] + nat_masses[j] == 0) || (nat_masses[i] + nat_masses[k] == 0) || (nat_masses[j] + nat_masses[k] == 0)) {
								k++;
								continue;
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
	}
	return {highest_prob_pos_4, combination_4};
}

std::tuple<vector<int>, int> findingCombinationOf4PTM ( float masses [],float prob [], const int &PTM_list, const float &mass_shift_lower, const float &mass_shift_upper) {


	// cout << "4 PTMs...";
	int list_end = PTM_list -1;
	vector<int> highest_prob_pos_4 = {list_end, list_end, list_end, list_end};
	int combination_4 = 0;


	float maximum = mass_shift_upper + (abs(masses[0])*3);//the maximum mass shift we need to check in the first loop 
	int limit {};
	if (maximum < masses[list_end]) 
		limit = binarySearchLargerElement(masses, maximum, PTM_list);
	else
		limit = list_end;


	// increasing loop
	if (mass_shift_upper < 160) {

		for (int y = 0; y < limit; y++) {
			float maximum_2 = -(masses[y]) + (abs(masses[0])*2) + mass_shift_upper;
			int limit_2 = binarySearchLargerElement( masses, maximum_2, PTM_list);


			for (int i = y; i < limit_2; i++) {
				float maximum_1 = -(masses[y]) -(masses[i]) + abs(masses[0]) + mass_shift_upper;
				int limit_1 = binarySearchLargerElement( masses, maximum_1, PTM_list);


				for (int j = i; j < limit_1; j++) {
					float minimum_1 {};
					minimum_1 = -(masses[y]) -(masses[i]) - (masses[j]) + mass_shift_lower;
					int minimum_limit_1 = binarySearchLargerElement (masses, minimum_1, PTM_list) -1;
					int k {};
					if (minimum_limit_1 < j)
						k = j;
					else
						k = minimum_limit_1;


					while (masses[y] + masses[i] + masses[j] + masses[k] <= mass_shift_upper && k < PTM_list) {
						float sum = roundf((masses[y] + masses[i] + masses[j] + masses[k])*10000)/10000;

						if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {
							// masses that cancel each other
							if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
								(masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
								k++;
								continue;
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
		for (int y = limit; y >= 0; y--) {
			float maximum_2 = -(masses[y]) + (abs(masses[0])*3) + mass_shift_upper;
			int limit_2 = binarySearchLargerElement( masses, maximum_2, PTM_list);


			for (int i = limit_2; i >= y; i--) {
				float maximum_1 = -(masses[y]) -(masses[i]) + (abs(masses[0])*2) + mass_shift_upper;
				int limit_1 = binarySearchLargerElement( masses, maximum_1, PTM_list);


				for (int j = limit_1; j >= i; j--) {
					float minimum_1 {};
					minimum_1 = -(masses[y]) -(masses[i]) - (masses[j]) + abs(masses[0]) + mass_shift_upper;
					int minimum_limit_1 = binarySearchLargerElement (masses, minimum_1, PTM_list) -1;
					int k {};
					if (minimum_limit_1 > j)
						k = j;
					else
						k = minimum_limit_1;


					while (masses[y] + masses[i] + masses[j] + masses[k] >= mass_shift_lower && k >= 0) {

						float sum = roundf((masses[y] + masses[i] + masses[j] + masses[k])*10000)/10000;
						//                    cout << k << "\n";
						//                    cout << sum << endl;
						if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

							// masses that cancel each other
							if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
								(masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
								k++;
								continue;
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
	return {highest_prob_pos_4, combination_4};
}

std::tuple<vector<int>, int, vector<int>, int> findingCombinationOf4and5PTM_A (float masses [],float prob [], const int &PTM_list, const float &mass_shift_lower, const float &mass_shift_upper) {

	int list_end = PTM_list -1;
	vector<int> highest_prob_pos_4 {list_end, list_end, list_end, list_end};
	int combination_4 = 0;
	vector<int> highest_prob_pos_5 {list_end, list_end, list_end, list_end, list_end};
	int combination_5 = 0;


	float maximum = mass_shift_upper + (abs(masses[0])*3);//the maximum mass shift we need to check in the first loop 
	int limit {};
	if (maximum < masses[list_end]) 
		limit = binarySearchLargerElement(masses, maximum, PTM_list);
	else
		limit = list_end;


	// increasing loop
	if (mass_shift_upper < 160) {
		for (int y = 0; y <= limit; y++) {
			float maximum_2 = -(masses[y]) + (abs(masses[0])*2) + mass_shift_upper;
			int limit_2 = binarySearchLargerElement( masses, maximum_2, PTM_list);


			for (int i = y; i <= limit_2; i++) {
				float maximum_1 = -(masses[y]) -(masses[i]) + (abs(masses[0])) + mass_shift_upper;
				int limit_1 = binarySearchLargerElement( masses, maximum_1, PTM_list);


				for (int j = i; j <= limit_1; j++) {
					float minimum_1 {};
					minimum_1 = -(masses[y]) -(masses[i]) - (masses[j]) + mass_shift_upper;
					int minimum_limit_1 = binarySearchLargerElement (masses, minimum_1, PTM_list) -1;
					int k {};
					if (minimum_limit_1 < j)
						k = j;
					else
						k = minimum_limit_1;


					while (masses[y] + masses[i] + masses[j] + masses[k] <= mass_shift_upper && k < PTM_list) {
						float sum = roundf((masses[y] + masses[i] + masses[j] + masses[k])*10000)/10000;
						if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {
							// masses that cancel each other
							if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
								(masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
								k++;
								continue;
							}
							if (prob[y] + prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_4[0]] + prob[highest_prob_pos_4[1]] + prob[highest_prob_pos_4[2]] + prob[highest_prob_pos_4[3]]) {
								highest_prob_pos_4[0] = y; highest_prob_pos_4[1] = i; highest_prob_pos_4[2] = j; highest_prob_pos_4[3] = k;
								//                            cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";
							}
							combination_4++;
						}
						float minimum_2 {};
						minimum_2 = -(masses[y]) -(masses[i]) - (masses[j]) - (masses[k]) + mass_shift_upper;
						int minimum_limit_2 = binarySearchLargerElement (masses, minimum_2, PTM_list) -1;
						int x {};
						if (minimum_limit_2 < k)
							x = k;
						else
							x = minimum_limit_2; 

						while (masses[y] + masses[i] + masses[j] + masses[k] + masses[x] >= mass_shift_lower && k < PTM_list) {
							float sum1 = roundf((masses[y] + masses[i] + masses[j] + masses[k] + masses[x])*10000)/10000;
							//                    cout << k << "\n";
							//                    cout << sum << endl;
							if ((sum1 <= mass_shift_upper) && (sum1 >= mass_shift_lower)) {

								// masses that cancel each other
								if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
									(masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0) ||
									(masses[x] + masses[y] == 0) || (masses[x] + masses[i] == 0) || (masses[x] + masses[j] == 0) ||
									(masses[x] + masses[k] == 0)) {
									x++;
									continue;
								}
							//                        cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";

							if (prob[y] + prob[i] + prob[j] + prob[k] + prob[x] >= prob[highest_prob_pos_5[0]] + prob[highest_prob_pos_5[1]] + prob[highest_prob_pos_5[2]] + prob[highest_prob_pos_5[3]] + prob[highest_prob_pos_5[4]]) {
								highest_prob_pos_5[0] = y; highest_prob_pos_5[1] = i; highest_prob_pos_5[2] = j; highest_prob_pos_5[3] = k; highest_prob_pos_5[4] = x;
							}
							combination_5++;
							}
							x++;
						}
						k++;
					}
				}
			}
		}

	} else {

		// decreasing loop
		for (int y = limit; y >= 0; y--) {
			float maximum_2 = -(masses[y]) + (abs(masses[0])*2) + mass_shift_upper;
			int limit_2 = binarySearchLargerElement( masses, maximum_2, PTM_list);


			for (int i = limit_2; i >= y; i--) {
				float maximum_1 = -(masses[y]) -(masses[i]) + abs(masses[0]) + mass_shift_upper;
				int limit_1 = binarySearchLargerElement( masses, maximum_1, PTM_list);


				for (int j = limit_1; j >= i; j--) {
					float minimum_1 {};
					minimum_1 = -(masses[y]) -(masses[i]) - (masses[j]) + mass_shift_lower;
					int minimum_limit_1 = binarySearchLargerElement (masses, minimum_1, PTM_list) -1;
					int k {};
					if (minimum_limit_1 > j)
						k = j;
					else
						k = minimum_limit_1;


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
						float minimum_2 {};
						minimum_2 = -(masses[y]) -(masses[i]) - (masses[j]) + abs(masses[k]) + mass_shift_upper;
						int minimum_limit_2 = binarySearchLargerElement (masses, minimum_2, PTM_list) -1;
						int x {};
						if (minimum_limit_2 > k)
							x = k;
						else
							x = minimum_limit_2; 

						while (masses[y] + masses[i] + masses[j] + masses[k] + masses[x] >= mass_shift_lower && k >= 0) {
							float sum = roundf((masses[y] + masses[i] + masses[j] + masses[k] + masses[x])*10000)/10000;
							//                    cout << k << "\n";
							//                    cout << sum << endl;
							if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

							// masses that cancel each other
								if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
									(masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0) ||
									(masses[x] + masses[y] == 0) || (masses[x] + masses[i] == 0) || (masses[x] + masses[j] == 0) ||
									(masses[x] + masses[k] == 0)) {
									x++;
									continue;
								}
							//                        cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";

							if (prob[y] + prob[i] + prob[j] + prob[k] + prob[x] >= prob[highest_prob_pos_5[0]] + prob[highest_prob_pos_5[1]] + prob[highest_prob_pos_5[2]] + prob[highest_prob_pos_5[3]] + prob[highest_prob_pos_5[4]]) {
								highest_prob_pos_5[0] = y; highest_prob_pos_5[1] = i; highest_prob_pos_5[2] = j; highest_prob_pos_5[3] = k; highest_prob_pos_5[4] = x;
							}
							combination_5++;
							}
							x--;
						}
						k--;
					}
				}
			}
		}
	}
	return {highest_prob_pos_4, combination_4, highest_prob_pos_5, combination_5};
}

std::tuple<vector<int>, int, vector<int>, int, float> findingCombinationOf4and5NaturalPTM (float nat_prob [], float nat_masses [],float prob [], const int &nat_prob_length, const float &mass_shift_lower, const float &mass_shift_upper,
						float &sum_of_probs_4, float &prob_average_1, float &prob_average_2, float &prob_average_3, float &prob_average_4, float &prob_vector_max_1) {

	int nat_end = nat_prob_length-1;
	// zero prob event
	vector<int> highest_prob_pos_4 = {nat_end, nat_end, nat_end, nat_end};
	vector<int> highest_prob_pos_5 = {nat_end, nat_end, nat_end, nat_end, nat_end};
	int combination_4 = 0;
	int combination_5 = 0;
	float sum_of_probs_5 = prob_vector_max_1 * 5;
		//if (sum_of_probs_4 < (4*nat_prob[0]) || sum_of_probs_5 < (5*nat_prob[0]))
		for (int y = 0; y < nat_prob_length; y++) {
//			float maximum_prob_5 = -(nat_prob[y]) - (nat_prob[y]*3) + sum_of_probs_5;
//			int limit_prob_5 = {};
//			if (maximum_prob_5 > nat_prob[y])
//					break; // if maximum is more than nat_prob[0] then it means that sum_of_probs will never be reached
//			else if ( maximum_prob_5 < 0)
//					limit_prob_5 = nat_prob_length; // if is 
//			else {
//					limit_prob_5 = binarySearchSmallerElement_decreasing (nat_prob, maximum_prob_5, nat_end);
//				}

			for (int x = y; x < nat_prob_length; x++) {
//				float maximum_prob = -(nat_prob[y]) -(nat_prob[x]) - (nat_prob[x]*2) + sum_of_probs_4;
//				int limit_prob {};
//				if (maximum_prob > nat_prob[x])
//					break; // if maximum is more than nat_prob[0] then it means that sum_of_probs will never be reached
//				else if ( maximum_prob < 0)
//					limit_prob = nat_prob_length; // if is 
//				else {
//					limit_prob = binarySearchSmallerElement_decreasing (nat_prob, maximum_prob, nat_end);
//				}

				for (int i = x; i < nat_prob_length; i++) {
//					float maximum_1 = -(nat_prob[y]) -(nat_prob[x]) - (nat_prob[i]*2) + sum_of_probs_4;
//					int limit1_prob {};
//					if (maximum_1 > nat_prob[i])
//						break; // if maximum is more than nat_prob[0] then it means that sum_of_probs will never be reached
//					else if ( maximum_1 < 0)
//						limit1_prob = nat_prob_length;
//					else {
//						limit1_prob = binarySearchSmallerElement_decreasing (nat_prob, maximum_1, nat_end);
//					}
					int j {i};

					while (nat_prob[y] + nat_prob[x] + nat_prob[i] + nat_prob[j] > sum_of_probs_4 && j < nat_prob_length) {
						
						float sum = roundf((nat_masses[x] + nat_masses[i] + nat_masses[j] + nat_masses[y])*10000)/10000;

                        if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

							// masses that cancel each other
							if ((nat_masses[x] + nat_masses[i] == 0) || (nat_masses[x] + nat_masses[j] == 0) || (nat_masses[x] + nat_masses[y] == 0) ||
							    (nat_masses[i] + nat_masses[j] == 0) || (nat_masses[i] + nat_masses[y] == 0) || (nat_masses[j] + nat_masses[y] == 0)) {
								j++;
								continue;
							}
							float prob_sum = nat_prob[x] + nat_prob[i] + nat_prob[j] + nat_prob[y];

							if (prob_sum >= nat_prob[highest_prob_pos_4[0]] + nat_prob[highest_prob_pos_4[1]] + nat_prob[highest_prob_pos_4[2]] + nat_prob[highest_prob_pos_4[3]]) {
								highest_prob_pos_4[0] = y; highest_prob_pos_4[1] = x; highest_prob_pos_4[2] = i; highest_prob_pos_4[3] = j;
								//cout << y << " " << x << " " << i << " " << j << endl;
								prob_average_4 = (nat_prob[highest_prob_pos_4[0]] + nat_prob[highest_prob_pos_4[1]] + nat_prob[highest_prob_pos_4[2]] + nat_prob[highest_prob_pos_4[3]]) / 4;
								vector<float> prob_vector_2 = {prob_average_1, prob_average_2, prob_average_3, prob_average_4};
								float prob_vector_max_2 = *max_element(prob_vector_2.begin(), prob_vector_2.end());
								sum_of_probs_5 = prob_vector_max_2*5;
							}
							combination_4++;
						}
						
						
						//float maximum_2 = -(nat_prob[y]) -(nat_prob[x]) - (nat_prob[i]) - (nat_prob[j]) + sum_of_probs_5;
						int k {j};
						//int limit2 {};
					
//						if (maximum_2 > nat_prob[j]) {
//							break; // if maximum is more than nat_prob[0] then it means that sum_of_probs will never be reached
//						}
//						else if ( maximum_2 < 0)
//							limit2 = nat_prob_length;
//						else {
//							limit2 = binarySearchSmallerElement_decreasing (nat_prob, maximum_2, nat_end);
//						}

						while (nat_prob[y] + nat_prob[x] + nat_prob[i] + nat_prob[j] + nat_prob[k] > sum_of_probs_5 && k < nat_prob_length) {

                            float sum1 = roundf((nat_masses[y] + nat_masses[x] + nat_masses[i] + nat_masses[j] + nat_masses[k])*10000)/10000;
                            
                            if ((sum1 <= mass_shift_upper) && (sum1 >= mass_shift_lower)) {

								// masses that cancel each other
								if ((nat_masses[x] + nat_masses[i]) == 0 || (nat_masses[x] + nat_masses[j] == 0) || (nat_masses[x] + nat_masses[y] == 0) ||
							    (nat_masses[i] + nat_masses[j] == 0) || (nat_masses[i] + nat_masses[y] == 0) || (nat_masses[j] + nat_masses[y] == 0) ||
								(nat_masses[x] + nat_masses[k] == 0) || (nat_masses[y] + nat_masses[k] == 0) || (nat_masses[i] + nat_masses[k] == 0) ||
								    (nat_masses[j] + nat_masses[k] == 0)) { k++; continue; }


								float prob_sum1 = nat_prob[y] + nat_prob[x] + nat_prob[i] + nat_prob[j] + nat_prob[k];

								if (prob_sum1 >=(nat_prob[highest_prob_pos_5[0]] + nat_prob[highest_prob_pos_5[1]] + nat_prob[highest_prob_pos_5[2]] + nat_prob[highest_prob_pos_5[3]] + nat_prob[highest_prob_pos_5[4]])) {
									highest_prob_pos_5[0] = y; highest_prob_pos_5[1] = x; highest_prob_pos_5[2] = i; highest_prob_pos_5[3] = j; highest_prob_pos_5[4] = k;
									//cout << y << " " << x << " " << i << " " << j << " " << k << endl;
								}
								combination_5++;
							}
							k++;
						}
					j++;
					} 
				}
			}
		}
		return {highest_prob_pos_4, combination_4, highest_prob_pos_5, combination_5, prob_average_4};
}

std::tuple<vector<int>, int, vector<int>, int> findingCombinationOf4and5PTM_Best (float masses [],float prob [], const int &PTM_list, const float &mass_shift_lower, const float &mass_shift_upper) {

	int list_end = PTM_list -1;
	vector<int> highest_prob_pos_4 {list_end, list_end, list_end, list_end};
	int combination_4 = 0;
	vector<int> highest_prob_pos_5 {list_end, list_end, list_end, list_end, list_end};
	int combination_5 = 0;


	float list_top_value = masses[list_end]; //maximum value of the list
	float list_lowest_value = masses[0]; //lowest value of the list
	float dif_upper_lower = mass_shift_upper - mass_shift_lower;


	float maximum = mass_shift_upper + (abs(list_lowest_value)*3);//the maximum mass shift we need to check in the first loop 
	int limit {};
	if (maximum < list_top_value) 
		limit = binarySearchLargerElement(masses, maximum, PTM_list);
	else
		limit = list_end;
	
	// increasing loop
	if (mass_shift_upper < 160) {
		for (int y = 0; y <= limit; y++) {
			float maximum_2 = -(masses[y]) + (abs(list_lowest_value)*2) + mass_shift_lower;
			int limit_2 {};
			int maximum_2_highest = maximum_2 + dif_upper_lower;
			if (maximum_2 > list_top_value || maximum_2_highest < list_lowest_value)
				break; //if its bigger then there is now way to achieve the mass_shift wanted. if the 
			else
				limit_2 = binarySearchLargerElement( masses, maximum_2, PTM_list);


			for (int i = y; i <= limit_2; i++) {
				float maximum_1 = -(masses[y]) -(masses[i]) + (abs(list_lowest_value)) + mass_shift_lower;
				int limit_1 {};
				int maximum_1_high = maximum_1 + dif_upper_lower;
				if (maximum_1 > list_top_value || maximum_1_high < list_lowest_value)
					break;
				else
					limit_1 = binarySearchLargerElement( masses, maximum_1, PTM_list);


				for (int j = i; j <= limit_1; j++) {
					float minimum_1 {};
					minimum_1 = -(masses[y]) -(masses[i]) - (masses[j]) + mass_shift_lower;
					int minimum_1_highest = minimum_1 + dif_upper_lower;
					int k {}; 
					
					if (minimum_1 > list_top_value || minimum_1_highest < list_lowest_value)
						break;
					else 
						k = binarySearchLargerElement (masses, minimum_1, PTM_list) -1;


					while (masses[y] + masses[i] + masses[j] + masses[k] <= mass_shift_upper && k < PTM_list) {
						float sum = roundf((masses[y] + masses[i] + masses[j] + masses[k])*10000)/10000;
						if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {
							// masses that cancel each other
							if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
								(masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
								k++;
								continue;
							}
							if (prob[y] + prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_4[0]] + prob[highest_prob_pos_4[1]] + prob[highest_prob_pos_4[2]] + prob[highest_prob_pos_4[3]]) {
								highest_prob_pos_4[0] = y; highest_prob_pos_4[1] = i; highest_prob_pos_4[2] = j; highest_prob_pos_4[3] = k;
								//                            cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";
							}
							combination_4++;
						}
						float minimum_2 {};
						minimum_2 = -(masses[y]) -(masses[i]) - (masses[j]) - (masses[k]) + mass_shift_lower;
						int minimum_2_highest = minimum_2 + dif_upper_lower;
						int x {};
						if (minimum_2_highest > list_top_value || minimum_2 < list_lowest_value)
							break;
						else
							x = binarySearchLargerElement (masses, minimum_2, PTM_list) -1; 


						while (masses[y] + masses[i] + masses[j] + masses[k] + masses[x] <= mass_shift_upper && k < PTM_list) {
							float sum1 = roundf((masses[y] + masses[i] + masses[j] + masses[k] + masses[x])*10000)/10000;
							//                    cout << k << "\n";
							//                    cout << sum << endl;
							if ((sum1 <= mass_shift_upper) && (sum1 >= mass_shift_lower)) {

								// masses that cancel each other
								if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
									(masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0) ||
									(masses[x] + masses[y] == 0) || (masses[x] + masses[i] == 0) || (masses[x] + masses[j] == 0) ||
									(masses[x] + masses[k] == 0)) {
									x++;
									continue;
								}
							//                        cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";

							if (prob[y] + prob[i] + prob[j] + prob[k] + prob[x] >= prob[highest_prob_pos_5[0]] + prob[highest_prob_pos_5[1]] + prob[highest_prob_pos_5[2]] + prob[highest_prob_pos_5[3]] + prob[highest_prob_pos_5[4]]) {
								highest_prob_pos_5[0] = y; highest_prob_pos_5[1] = i; highest_prob_pos_5[2] = j; highest_prob_pos_5[3] = k; highest_prob_pos_5[4] = x;
							}
							combination_5++;
							}
							if (x <= k)
								x = k + 1;
							else
								x++;
						}
						if ( k <= j)
							k = j + 1;
						else
							k++;
					}
				}
			}
		}

	} else {

		// decreasing loop
		for (int y = limit; y >= 0; y--) {
			float maximum_2 = -(masses[y]) + (abs(list_lowest_value)*2) + mass_shift_upper;
			int limit_2 {};
			int maximum_2_lower = maximum_2 - dif_upper_lower;
			if (maximum_2_lower > list_top_value || maximum_2 < list_lowest_value)
				break;
			else
				limit_2 = binarySearchLargerElement( masses, maximum_2, PTM_list);


			for (int i = limit_2; i >= y; i--) {
				float maximum_1 = -(masses[y]) -(masses[i]) + abs(list_lowest_value) + mass_shift_upper;
				int maximum_1_lowest = maximum_1 - dif_upper_lower;
				int limit_1 {};
				if (maximum_1_lowest > list_top_value || maximum_1 < list_lowest_value)
					break;
				else
					limit_1= binarySearchLargerElement( masses, maximum_1, PTM_list);


				for (int j = limit_1; j >= i; j--) {
					float minimum_1 {};
					minimum_1 = -(masses[y]) -(masses[i]) - (masses[j]) + mass_shift_upper;
					int minimum_1_lowest = minimum_1 - dif_upper_lower;
					int k {};
					if (minimum_1_lowest > list_top_value || minimum_1 < list_lowest_value)
						break;
					else
						k = binarySearchLargerElement (masses, minimum_1, PTM_list) -1;


					while (masses[y] + masses[i] + masses[j] + masses[k] >= mass_shift_lower && y >= 0) {
						float sum = roundf((masses[y] + masses[i] + masses[j] + masses[k])*10000)/10000;
						//                    cout << k << "\n";
						//                    cout << sum << endl;
						if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

							// masses that cancel each other
							if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
								(masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
								y--;
								continue;
							}
							//                        cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";

							if (prob[y] + prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_4[0]] + prob[highest_prob_pos_4[1]] + prob[highest_prob_pos_4[2]] + prob[highest_prob_pos_4[3]]) {
								highest_prob_pos_4[0] = y; highest_prob_pos_4[1] = i; highest_prob_pos_4[2] = j; highest_prob_pos_4[3] = k;
							}
							combination_4++;
						}
						float minimum_2 {};
						minimum_2 = -(masses[y]) -(masses[i]) - (masses[j]) - (masses[k]) + mass_shift_upper;
						
						int minimum_2_lowest = minimum_2 - dif_upper_lower;
  
						int x {};
						if (minimum_2_lowest > list_top_value || minimum_2 < list_lowest_value)
							break;
						else
							x = binarySearchLargerElement (masses, minimum_2, PTM_list) -1; 

						while (masses[y] + masses[i] + masses[j] + masses[k] + masses[x] >= mass_shift_lower && x >= 0) {
							float sum = roundf((masses[y] + masses[i] + masses[j] + masses[k] + masses[x])*10000)/10000;
							//                    cout << k << "\n";
							//                    cout << sum << endl;
							if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

							// masses that cancel each other
								if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
									(masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0) ||
									(masses[x] + masses[y] == 0) || (masses[x] + masses[i] == 0) || (masses[x] + masses[j] == 0) ||
									(masses[x] + masses[k] == 0)) {
									x--;
									continue;
								}
							//                        cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";

							if (prob[y] + prob[i] + prob[j] + prob[k] + prob[x] >= prob[highest_prob_pos_5[0]] + prob[highest_prob_pos_5[1]] + prob[highest_prob_pos_5[2]] + prob[highest_prob_pos_5[3]] + prob[highest_prob_pos_5[4]]) {
								highest_prob_pos_5[0] = y; highest_prob_pos_5[1] = i; highest_prob_pos_5[2] = j; highest_prob_pos_5[3] = k; highest_prob_pos_5[4] = x;
							}
							combination_5++;
							}
							if (x==k)
								x = k - 1;
							else
								x--;
						}
						if (k ==j)
							k = j -1;
						else
							k--;
					}
				}
			}
		}
	}
	return {highest_prob_pos_4, combination_4, highest_prob_pos_5, combination_5};
}

std::tuple<vector<int>, int, vector<int>, int, vector<int>, int, vector<int>, int> findingCombinationOf2_3_4_5_PTM ( float masses [],float prob [], const int &PTM_list, const float &mass_shift_lower, const float &mass_shift_upper) {
	
	vector<int> position_combination2 {}; // here is were we store all matched positions
	int list_end = PTM_list-1;

	vector<int> highest_prob_pos_2 {list_end, list_end};
	int combination_2 = 0;
	vector<int> highest_prob_pos_3 {list_end, list_end, list_end};
	int combination_3 = 0;
	vector<int> highest_prob_pos_4 {list_end, list_end,list_end, list_end,list_end};
	int combination_4 = 0;
	vector<int> highest_prob_pos_5 {list_end, list_end, list_end, list_end, list_end};
	int combination_5 = 0;
	
	float list_top_value = masses[list_end]; //maximum value of the list
	float list_lowest_value = masses[0]; //lowest value of the list
	

	float maximum = mass_shift_upper + (abs(list_lowest_value));//the maximum mass shift we need to check in the first loop 
	int limit {};
	if (maximum < list_top_value) 
		limit = binarySearchLargerElement(masses, maximum, PTM_list);
	else
		limit = list_end;
		
	// increasing loop
	if (mass_shift_upper < 160) {
			for (int i = 0; i < limit; i++) {
				int j {i};
				
				
				while (masses[i] + masses[j] <= mass_shift_upper && j < PTM_list) {
						
					float sum = roundf((masses[i] + masses[j])*10000)/10000;
                
					if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {
					//cout << masses[i] << " " << masses[j] << "\n";
					// masses that cancel each other
					if (masses[i] + masses[j] == 0) {
						j++;
						continue;//does break statement make sense here?? think not because we can have more matches
					}
					else {
						//position_combination2.push_back(i);
						//position_combination2.push_back(j);
//                 cout << masses[j] << " " << masses[k] << "\n";
						if (prob[i] + prob[j] >= prob[highest_prob_pos_2[0]] + prob[highest_prob_pos_2[1]]) {
							highest_prob_pos_2[0] = i; highest_prob_pos_2[1] = j;
							//cout << masses[i] << " " << masses[j] << "\n";
						}
					combination_2++;
					}
					}
					int k {j};


					while (masses[i] + masses[j] + masses[k] <= mass_shift_upper && k < PTM_list) {
							float sum1 = roundf((masses[i] + masses[j] + masses[k])*10000)/10000;

							if ((sum1 <= mass_shift_upper) && (sum1 >= mass_shift_lower)) {

							// masses that cancel each other
							if ((masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
								k++;
								continue;
							}
							if ((prob[i] + prob[j] + prob[k]) >= (prob[highest_prob_pos_3[0]] + prob[highest_prob_pos_3[1]] + prob[highest_prob_pos_3[2]])) {
								highest_prob_pos_3[0] = i; highest_prob_pos_3[1] = j; highest_prob_pos_3[2] = k;
								//cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";
							}
							combination_3++;
							}
							int y {k};


							while (masses[y] + masses[i] + masses[j] + masses[k] <= mass_shift_upper && k < PTM_list) {
								float sum2 = roundf((masses[y] + masses[i] + masses[j] + masses[k])*10000)/10000;
								if ((sum2 <= mass_shift_upper) && (sum2 >= mass_shift_lower)) {
									// masses that cancel each other
									if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
										(masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
										y++;
										continue;
									}
									if (prob[y] + prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_4[0]] + prob[highest_prob_pos_4[1]] + prob[highest_prob_pos_4[2]] + prob[highest_prob_pos_4[3]]) {
										highest_prob_pos_4[0] = y; highest_prob_pos_4[1] = i; highest_prob_pos_4[2] = j; highest_prob_pos_4[3] = k;
										//cout << masses[i] << " " << masses[j] << " " << masses[k] << " "<< masses[y] << "\n";
									}
									combination_4++;
								}
								int x {y};
								


								while (masses[y] + masses[i] + masses[j] + masses[k] + masses[x] <= mass_shift_upper && k < PTM_list) {
									float sum3 = roundf((masses[y] + masses[i] + masses[j] + masses[k] + masses[x])*10000)/10000;
									//                    cout << k << "\n";
									//                    cout << sum << endl;
									if ((sum3 <= mass_shift_upper) && (sum3 >= mass_shift_lower)) {

										// masses that cancel each other
										if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
											(masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0) ||
											(masses[x] + masses[y] == 0) || (masses[x] + masses[i] == 0) || (masses[x] + masses[j] == 0) ||
											(masses[x] + masses[k] == 0)) {
											x++;
											continue;
										}
									//                        cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";

									if (prob[y] + prob[i] + prob[j] + prob[k] + prob[x] >= prob[highest_prob_pos_5[0]] + prob[highest_prob_pos_5[1]] + prob[highest_prob_pos_5[2]] + prob[highest_prob_pos_5[3]] + prob[highest_prob_pos_5[4]]) {
										highest_prob_pos_5[0] = y; highest_prob_pos_5[1] = i; highest_prob_pos_5[2] = j; highest_prob_pos_5[3] = k; highest_prob_pos_5[4] = x;
										//cout << masses[i] << " " << masses[j] << " " << masses[k] << " "<< masses[y] << " " << masses[x] << "\n";
									}
									combination_5++;
									}
								
										x++;
								}
								
									y++;
							}
							
								k++;
					}
					
						j++;
				}
			}
	} else {
		// decreasing loop
		for (int i = limit; i >= 0; i--) {
			int j {i};

			while (masses[i] + masses[j] >= mass_shift_lower && j >= 0) {
				float sum = roundf((masses[i] + masses[j])*10000)/10000;
                
					if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {
					
					// masses that cancel each other
						if (masses[i] + masses[j] == 0) {
							j++;
							continue; //does break statement make sense here?? think not because we can have more matches
						} 
						else {
							//position_combination2.push_back(i); // can use insert below to keep the most probable in 0 and 1;
							//position_combination2.push_back(j);
							//cout << masses[j] << " " << masses[k] << "\n";
							if (prob[i] + prob[j] >= prob[highest_prob_pos_2[0]] + prob[highest_prob_pos_2[1]]) {
								highest_prob_pos_2[0] = i; highest_prob_pos_2[1] = j;
								//cout << masses[i] << " " << masses[j] << "\n";
							}
							combination_2++;
						}
					}
				int k {j};

				while (masses[i] + masses[j] + masses[k] >= mass_shift_lower && k >= 0) {

					float sum = roundf((masses[i] + masses[j] + masses[k])*10000)/10000;

					if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

					// masses that cancel each other
						if ((masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
							k--;
							continue;
						}
						if (prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_3[0]] + prob[highest_prob_pos_3[1]] + prob[highest_prob_pos_3[2]] ) {
							highest_prob_pos_3[0] = i; highest_prob_pos_3[1] = j; highest_prob_pos_3[2] = k;
							//cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";
						}
						combination_3++;
					}
				
					int y {k};
				 
				

					while (masses[y] + masses[i] + masses[j] + masses[k] >= mass_shift_lower && y >= 0) {
						float sum = roundf((masses[y] + masses[i] + masses[j] + masses[k])*10000)/10000;
						//                    cout << k << "\n";
						//                    cout << sum << endl;
						if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

							// masses that cancel each other
							if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
								(masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0)) {
								y--;
								continue;
							}
							//                        cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";

							if (prob[y] + prob[i] + prob[j] + prob[k] >= prob[highest_prob_pos_4[0]] + prob[highest_prob_pos_4[1]] + prob[highest_prob_pos_4[2]] + prob[highest_prob_pos_4[3]]) {
								highest_prob_pos_4[0] = y; highest_prob_pos_4[1] = i; highest_prob_pos_4[2] = j; highest_prob_pos_4[3] = k;
								//cout << masses[i] << " " << masses[j] << " " << masses[k] << " "<< masses[y] << " " << "\n";
							}
							combination_4++;
						}
	
						int x {y};
//						float minimum_3 {};
//						minimum_3 = -(masses[y]) -(masses[i]) - (masses[j]) - (masses[k]) + mass_shift_upper;
//						float minimum_3_highest = minimum_3 + dif_upper_lower;
//						if (minimum_3_highest < list_lowest_value || minimum_3 > list_top_value)
//							break;
//						else
//							x = binarySearchLargerElement (masses, minimum_3_highest, PTM_list); 

						while (masses[y] + masses[i] + masses[j] + masses[k] + masses[x] >= mass_shift_lower && x >= 0) {
							float sum = roundf((masses[y] + masses[i] + masses[j] + masses[k] + masses[x])*10000)/10000;
							//                    cout << k << "\n";
							//                    cout << sum << endl;
							if ((sum <= mass_shift_upper) && (sum >= mass_shift_lower)) {

							// masses that cancel each other
								if ((masses[y] + masses[i] == 0) || (masses[y] + masses[j] == 0) || (masses[y] + masses[k] == 0) ||
									(masses[i] + masses[j] == 0) || (masses[i] + masses[k] == 0) || (masses[j] + masses[k] == 0) ||
									(masses[x] + masses[y] == 0) || (masses[x] + masses[i] == 0) || (masses[x] + masses[j] == 0) ||
									(masses[x] + masses[k] == 0)) {
									x--;
									continue;
								}
							//                        cout << masses[i] << " " << masses[j] << " " << masses[k] << "\n";

							if (prob[y] + prob[i] + prob[j] + prob[k] + prob[x] >= prob[highest_prob_pos_5[0]] + prob[highest_prob_pos_5[1]] + prob[highest_prob_pos_5[2]] + prob[highest_prob_pos_5[3]] + prob[highest_prob_pos_5[4]]) {
								highest_prob_pos_5[0] = y; highest_prob_pos_5[1] = i; highest_prob_pos_5[2] = j; highest_prob_pos_5[3] = k; highest_prob_pos_5[4] = x;
								//cout << masses[i] << " " << masses[j] << " " << masses[k] << " "<< masses[y] << " " << masses[x] << "\n";
							}
							combination_5++;
							}
						
								x--;
						}
					
							y--;
					}
				
						k--;
				}
				
					j--;
			}
		}
	}
	
	//cout << "---------------------------------------------------------------------------------------------------------------------------------";
		return {highest_prob_pos_2, combination_2, highest_prob_pos_3, combination_3, highest_prob_pos_4, combination_4, highest_prob_pos_5, combination_5};
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



int main(int argc, char** argv) {

	/* Main arguments:
	 * first = sample name
	 * second = tolerance
	 */ 
	auto start = chrono::steady_clock::now();
	


	
	
	string sample = "JOHNPP_l1";

	string dir_path = "/Users/pedrocardoso/Documents/RUHI-MSA/";
	string input_path = dir_path + "input/" + sample + ".txt";
	
	
	
	

	
	//-------------------------------- Read Input files ------------------------------------------//
/* this is probably not needed since we are reading the file now with a while loop
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

//    cout <<  "No of Peptides: " << lines << "\n";
*/
	//--------------------------------//



	//================================| Input |================================//


	//// PTM masses and probabilities

	// mass shift array
	const int PTM_list = 1498;  //1497 PTMs
	// list_end is -1 because on an array with 10 element position 9 is the last element!!!!
	const int list_end = PTM_list; //last element to create null prob for no matches
	const int nat_list = 90; //89 natural PTMs
	const int nat_end = nat_list;

	//// UniMod masses and probabilities

	// mass shift array // PTMs array
	float* masses = openFile(dir_path + "libs/new_masses.txt", PTM_list);

	// Probabilities array
	float* prob = openFile(dir_path + "libs/new_prob.txt", (PTM_list+1));

	// null probability for no matches (last new vector position)
	prob[list_end] = 0;


	//// Natural frequency masses and probabilities

	// Nat probabilities array
	float* nat_prob = openFile(dir_path + "libs/new_nat_prob.txt", (nat_list+1));

	nat_prob[nat_end] = 0;

	// Nat masses array
	float* nat_masses = openFile(dir_path + "libs/new_nat_masses.txt", nat_list);


	// Nat names array
	string nat_names[nat_list];
	ifstream fakeNames;
	fakeNames.open(dir_path + "libs/new_nat_names.txt");

	int names_count = 0;
	while (!fakeNames.eof()) {
		fakeNames >> nat_names[names_count];
		names_count++;
	}

	// sizes
	
	int nat_prob_length = nat_list;

	//-------------------------------------------------------//

	//-----------OUTPUT file--------------------------------//
	string output_file_path = dir_path + "results/" + sample + "_output.txt";
	
	ofstream output_file;
	output_file.open(output_file_path, ios::out);

	output_file << "Peptide" << "\t" << "Mass_shift" << "\t" << "Peptide Mass" << "\t" << "Mass of PTMs" << "\t" << "Score" << "\t" <<  "PTM 1" << "\t"<< "PTM 2" << "\t"<< "PTM 3" << "\t"<< "PTM4" << "\t" << "PTM5";
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
		
		//outfile << peptide << " " << mass_shift << " " << peptide_mass << endl;		
//		auto end1 = chrono::steady_clock::now();
        // mass error
//        float peptide_mass = 1000; 
        float tolerance = 10; // this has been defined on the  .main arguments
        float mass_error = peptide_mass * tolerance / 1000000;
        //cout << mass_error << endl;
        float mass_shift_upper = mass_shift + mass_error;
        float mass_shift_lower = mass_shift - mass_error;
        //cout << mass_shift_lower << "-" << mass_shift_upper << endl;
        mass_shift_upper = roundf(mass_shift_upper*10000)/10000; //check if this works without rounding. If we round we are losing precision
        mass_shift_lower = roundf(mass_shift_lower*10000)/10000;
		//cout << mass_shift_lower << "-" << mass_shift_upper << endl;
        

		//================================//

		//string output_file_path = dir_path + "results/" + sample + "/" + peptide + "_output.txt";

		// int break_scenario=0;

		//----------------------------------//
		// string peptide = "VSAMEDEMNEMK";
		// float mass_shift = -0.0011;
		//----------------------------------//

		// counters
		//    int count;
		int combination = 0;
		int more_than_3 = 0;
		//=======================================================================================================//

		// start timer
		//auto start = chrono::steady_clock::now();

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
		//––––––––| 4 and/or 5 NAtural PTMs |––––––––//
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

			output_file << peptide << "\t" << mass_shift << "\t" << peptide_mass << "NA" << "NA" << "NA" << "NA" << "NA" << "\n";
			output_file.close();


			//// skip to next peptide
			continue;
		}

		//=======================================================================================================//


		//================================| Output |================================//

		// timers
		
//		auto end = chrono::steady_clock::now();
//		auto diff = end - start;
//		string elapsed_time;
//
//		if (chrono::duration <double, milli> (diff).count() < 1) {
//			elapsed_time = to_string(chrono::duration <double, milli> (diff).count()) + " ms";
//		} else {
//			elapsed_time = to_string(chrono::duration <double, milli> (diff).count() / 1000) + " s";
//		}
		

		//// print output combination stats
//		cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << "\n\n";
//		cout << "peptide:  " << peptide << "\n\n";
//		cout << "PTM list:   " << length << ";  number of PTMs: " << number_of_possible_ptms << "\n";
//		cout << "mass shift: " << mass_shift << " ± " << mass_error << ";  combinations: " << combination << "\n" << "\n";
//		cout << "=======================================================================" << "\n" << "\n";


//    cout << prob_average_3 << " " <<prob_average_4 << "\n";
		vector<string> out_PTM_name;
		float sum_of_masses {}; // mass of all PTMS
		float prob_PTMs {}; // score of PTMs
		//// 1-3 PTMs
		if (more_than_3 <= 3) {
			//// print PTM names from file
			fstream inputPTMnames;
			inputPTMnames.open(dir_path + "new_names.txt");
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

			//cout << "\n\nmass shift sum: ";
//			for (int i = 0; i < 4; i++) {
//				if (highest_prob_pos[i] != list_end) {
//					cout << masses[highest_prob_pos[i]] << "\t";
//				}
			


			// print positions
//			cout << "\n\npositions: ";
//
//			for (int i = 0; i < 4; i++) {
//				if (highest_prob_pos[i] != list_end) {
//					cout << 1 + highest_prob_pos[i] << " ";
//				}
//			}


		//// 4-5 PTMs
		} else if (more_than_3 <= 5) {
			
			fstream Natural_PTM_Names;
			Natural_PTM_Names.open(dir_path + "new_nat_names.txt");
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


//			cout << "\n\nmass shift sum:  ";
//			for (int i = 0; i <= 4; i++) {
//				if (highest_prob_pos[i] != list_end) {
//					cout << nat_masses[highest_prob_pos[i]] << "\t";
//				}
//			}


//			cout << "\npositions: ";
//			for (int i = 0; i <= 4; i++) {
//				if (highest_prob_pos[i] != list_end) {
//					cout << 1 + highest_prob_pos[i] << " ";
//				}
//			}
		} else if (more_than_3 == 6) {

			//// print PTM names from file
			fstream inputPTMnames;
			inputPTMnames.open(dir_path + "new_names.txt");
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

//			cout << "\n\nmass shift sum: ";
//			for (int i = 0; i <= 4; i++) {
//				if (highest_prob_pos[i] != list_end) {
//					cout << masses[highest_prob_pos[i]] << "\t";
//				}
//			}
//
//
//			// print positions
//			cout << "\n\npositions: ";
//
//			for (int i = 0; i <= 4; i++) {
//				if (highest_prob_pos[i] != list_end) {
//					cout << 1 + highest_prob_pos[i] << " ";
//				}
//			}
		}


		// print run-time
		//cout << "\n\n" << "execution time = " << elapsed_time << "\n\n";
		//cout << "=======================================================================" << "\n";


		//cout <<  prob_average_1 << " " << prob_average_2 << " " << prob_average_3 << " " << prob_average_4 << " " << prob_average_5 << "\n";
		//--------------------------------------------------------------------------------------------------/


		//================================| Write Output |================================//

		//// write output into file
		//ofstream output_file;
		output_file.open(output_file_path, ios::out | ios::app); //ios::app to open file in append mode

			//output_file << "Peptide" << "\t" << "Mass_shift" << "\t" << "field" << "\t"; //DONE at the beginning
//			for (int i = 1; i < out_PTM_name.size(); i++) {
//				output_file << "PTM " << i << "\t";
//			}
			

//"Peptide" << "\t" << "Mass_shift" << "\t" << "Peptide Mass" << "\t" << "Mass of PTMs" << "\t" << "Score" << "\t" << "PTM 1" << "PTM 2" etc...;
			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << peptide_mass << "\t" << sum_of_masses << "\t" << prob_PTMs << "\t";
			for (int i = 0; i < out_PTM_name.size(); i++) {
				output_file << out_PTM_name[i] << "\t";
			}
		//output_file << "\n";
		output_file.close();


//			for (int i = 0; i < out_name.size(); i++) {
//				if (highest_prob_pos[i] != list_end) {
//					output_file << masses[highest_prob_pos[i]] << "\t";
//				}
//			}
//
//			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
//			for (int i = 0; i < out_name.size(); i++) {
//				if (highest_prob_pos[i] != list_end) {
//					output_file << prob[highest_prob_pos[i]] << "\t";
//				}
//			}

//			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
//			for (int i = 0; i < out_name.size(); i++) {
//				if (highest_prob_pos[i] != list_end) {
//					output_file << 1 + highest_prob_pos[i] << "\t";
//				}
//			}



			//// 4-5 PTMs
//		 else if (more_than_3 <= 5) {
//
//			if (more_than_3 == 4) {
//				output_file <<  "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t" << "1" << "\t" << "2" << "\t" << "3" << "\t" << "4" << "\t";
//			} else if (more_than_3 == 5) { output_file << "peptide" << "\t" << "mass_shift" << "\t" << "field" << "\t" << "1" << "\t" << "2" << "\t" << "3" << "\t" << "4" << "\t" << "5" << "\t"; }
//
//
//			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "name" << "\t";
//			for (int i = 0; i <= 4; i++) {
//				if (highest_prob_pos[i] != list_end) {
//					output_file << nat_names[highest_prob_pos[i]] << "\t";
//				}
//			}
//
//			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "mass" << "\t";
//			for (int i = 0; i <= 4; i++) {
//				if (highest_prob_pos[i] != list_end) {
//					output_file << nat_masses[highest_prob_pos[i]] << "\t";
//				}
//			}
//
//			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "score" << "\t";
//			for (int i = 0; i <= 4; i++) {
//				if (highest_prob_pos[i] != list_end) {
//					output_file << nat_prob[highest_prob_pos[i]] << "\t";
//				}
//			}
//
//			output_file << "\n" << peptide << "\t" << mass_shift << "\t" << "position" << "\t";
//			for (int i = 0; i <= 4; i++) {
//				if (highest_prob_pos[i] != list_end) {
//					output_file << 1 + highest_prob_pos[i] << "\t";
//				}
//			}
//		}

		

		//cout << "\n\n";

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


