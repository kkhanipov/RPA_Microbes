#ifndef PCR_Profile_ARRAY_H
#define PCR_Profile_ARRAY_H

#include <iostream>
#include "PCR_Profile.h"

class PCR_Profile_Array
{
	PCR_Profile ** pcr_profiles;
	int * to_be_considered; //0=added, 1=pareto, -1=marked deletion
	unsigned int number_of_profiles;
	unsigned int array_size;

public:
	PCR_Profile_Array();
	~PCR_Profile_Array();
	bool add_pcr_profile(PCR_Profile * pcr_profile);
	bool delete_pcr_profile(unsigned int position); 
	bool mark_profile_pareto(unsigned int position);
	bool get_stats();
	bool compress_left();
	PCR_Profile * get_pcr_profile(unsigned int position);
	int get_number_of_profiles() { return number_of_profiles;}
	int * get_to_be_considered() { return to_be_considered; }
};
bool PCR_Profile_Array::get_stats()
{
	cout << "Number of Profiles " << number_of_profiles << endl;
	return true;
}
PCR_Profile_Array::PCR_Profile_Array()
{
	array_size = 4096;
	pcr_profiles = new PCR_Profile *[array_size]; assert(pcr_profiles);
	to_be_considered = new int [array_size]; assert(to_be_considered);
	number_of_profiles = 0;
}
PCR_Profile_Array::~PCR_Profile_Array()
{
	for (int i = 0; i < number_of_profiles; i++)
	{
		if (pcr_profiles[i] != NULL) delete pcr_profiles[i];
	}

	delete[]pcr_profiles;
	delete[]to_be_considered;
}

bool PCR_Profile_Array::add_pcr_profile(PCR_Profile * _pcr_profile)
{
	if (number_of_profiles + 1 >= array_size) //array size too small, needs to be reallocated
	{
		array_size *= 1.5;
		to_be_considered = (int*)realloc(to_be_considered, sizeof(int)*array_size);
		pcr_profiles = (PCR_Profile **)realloc(pcr_profiles, sizeof(PCR_Profile *)*array_size);
	}
	pcr_profiles[number_of_profiles] = new PCR_Profile(_pcr_profile);
	to_be_considered[number_of_profiles] = 0;
	number_of_profiles++;
	return true;
}
bool PCR_Profile_Array::delete_pcr_profile(unsigned int position)
{
	if (position > number_of_profiles)
	{
		cout << "PCR_Profile_Array::delete_pcr_profile, that position does not exist" << endl;
		assert (position < number_of_profiles);
	}
	if (pcr_profiles[position]==NULL)
	{
		cout << "PCR_Profile_Array::delete_pcr_profile, that position has been deleted" << endl;
		assert(pcr_profiles[position] == NULL);
	}
	delete pcr_profiles[position];
	pcr_profiles[position] = NULL;
	to_be_considered[position] = -1;
	return true;
}

bool PCR_Profile_Array::mark_profile_pareto(unsigned int position)
{
	if (position > number_of_profiles)
	{
		cout << "PCR_Profile_Array::mark_profile_pareto position does not exist" << endl;
		assert(position < number_of_profiles);
	}
	if (pcr_profiles[position] == NULL)
	{
		cout << "PCR_Profile_Array::mark_profile_pareto, that position has been deleted" << endl;
		assert(pcr_profiles[position] == NULL);
	}
	to_be_considered[position] = 1;
	return true;
}
PCR_Profile* PCR_Profile_Array::get_pcr_profile(unsigned int position)
{
	if (position > number_of_profiles)
	{
		cout << "PCR_Profile* PCR_Profile_Array::get_pcr_profile position does not exist" << endl;
		assert(position < number_of_profiles);
	}
	if (pcr_profiles[position] == NULL)
	{
		cout << "PCR_Profile* PCR_Profile_Array::get_pcr_profile, that position has been deleted" << endl;
		assert(pcr_profiles[position] == NULL);
	}
	return pcr_profiles[position];
}

bool PCR_Profile_Array::compress_left()
{
	PCR_Profile ** temp_pcr_profiles;
	temp_pcr_profiles= new PCR_Profile *[array_size]; assert(pcr_profiles);
	int * temp_to_be_considered;
	temp_to_be_considered = new int[array_size]; assert(to_be_considered);
	unsigned int temp_num_of_profiles = 0;

	for (int i = 0; i < number_of_profiles; i++)
	{
		if (to_be_considered[i] == 1)
		{
			temp_to_be_considered[temp_num_of_profiles] = to_be_considered[i];
			temp_pcr_profiles[temp_num_of_profiles] = pcr_profiles[i];
			temp_num_of_profiles++;
		}
	}

	delete pcr_profiles;
	pcr_profiles=temp_pcr_profiles;

	delete to_be_considered;
	to_be_considered = temp_to_be_considered;

	number_of_profiles = temp_num_of_profiles;
	return true;
}

#endif

