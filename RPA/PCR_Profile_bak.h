#ifndef PCR_Profile_H
#define PCR_Profile_H


#include "Primer_Set.h"
#include "Sequence.h"
#include <iostream>

using namespace std;
// Object PCR profile will be used to store the primer sets and information about how they are mapping to the sequence of interest

class PCR_Profile
{
	static unsigned int min_primer_distance;
	static unsigned int max_primer_distance;
	Primer_Set * p_set; // this will be the primer set used to create the primer locations profile

	unsigned int profile_length; //length of the primer_locations array, should be the size of the sequence of interest

	int * location_of_primer;
	int * type_of_primer;
	unsigned int primer_profile_array_size;
	unsigned int number_of_primers;
	// an array the size of the seqeunce of interest to denote where the primers from the p_set are mapping. +1 will denote forward mapping primer location, 
	//0 will denote a non mapped location, -1 will denote a reverse mapping primer location, -5 will denote positions which are not available due to unknown nucletodies 
	//being present there, 5 will denote if there are primers in both directions at the same location.
	

	struct Stats
	{
		unsigned int number_forward_primers; // calculated by calculate_statistics(), contains the number of forward primers in the primer locations profile
		unsigned int number_reverse_primers;// calculated by calculate_statistics(), contains the number of reverse primers in the primer locations profile
		unsigned int number_short_amplicons; // calculated by calculate_statistics(), contains the number of amplicons shorter than designed by a global variable in the primer locations profile
		unsigned int number_long_amplicons;// calculated by calculate_statistics(), contains the number of amplicons of a good size than designed by a global variable in the primer locations profile
		unsigned int total_lenght_short_amplicons; // calculated by calculate_statistics(),  total length of the primer locations profile covered by short amplicons
		unsigned int total_lenght_long_amplicons;// calculated by calculate_statistics(),  total length of the primer locations profile covered by long amplicons
		unsigned int total_lenght_too_long_amplicons;
		unsigned int total_length_uncovered;
	} stats;

public:
	PCR_Profile(PCR_Profile * _pcr_profile_a, PCR_Profile * _pcr_profile_b, ostream &err_msg = cout);
	PCR_Profile(PCR_Profile * _pcr_profile, ostream &err_msg = cout);
	PCR_Profile(Primer_Set * _primer_set, Sequence * _sequence, ostream &err_msg = cout);
	
	~PCR_Profile();
	bool PCR_profile_calculation(Sequence * _sequence, ostream &err_msg = cout);
	bool calculate_statistics(ostream & out = cout, ostream &err_msg = cout); //calculates statistics on the primer_locations array

	bool show_statistics(ostream & out = cout, ostream &err_msg = cout); // gives the values of the calculate_statistics();
	bool show_All(ostream & out = cout, ostream &err_msg = cout); // first shows the show_statistics() and then prints the primer_location profile
	bool add_primer_location_to_profile(unsigned int position, int primer_type, ostream &err_msg = cout);
	unsigned int get_number_forward_primers() { return stats.number_forward_primers;}
	unsigned int get_number_reverse_primers() { return stats.number_reverse_primers; }
	unsigned int get_number_short_amplicons() { return stats.number_short_amplicons; }
	unsigned int get_number_long_amplicons() { return stats.number_long_amplicons; }
	unsigned int get_total_lenght_short_amplicons() { return stats.total_lenght_short_amplicons; }
	unsigned int get_total_lenght_long_amplicons() { return stats.total_lenght_long_amplicons; }
	unsigned int get_profile_length() { return profile_length; }
	unsigned int get_total_lenght_too_long_amplicons() { return stats.total_lenght_too_long_amplicons; }
	unsigned int get_total_length_uncovered() {	return stats.total_length_uncovered;}
	Primer_Set * get_pointer_to_primer_set() { return p_set; }
	unsigned int get_number_of_primers_primer_set() { return p_set->get_number_of_primers(); }
	unsigned int get_primer_profile_array_size() { return primer_profile_array_size; }
	unsigned int get_number_of_primers_location_profile() { return number_of_primers; }
	int * get_location_of_primer() {return location_of_primer;}
	int * get_type_of_primer() { return type_of_primer; }


	Stats get_Stats() { return stats; }
};

unsigned int PCR_Profile::min_primer_distance = 100;
unsigned int PCR_Profile::max_primer_distance = 500;

PCR_Profile::~PCR_Profile()
{
	delete p_set;
	delete[] location_of_primer;
	delete[] type_of_primer;
}
PCR_Profile::PCR_Profile(PCR_Profile * _pcr_profile_a, PCR_Profile * _pcr_profile_b, ostream &err_msg)
{
	stats = {};
	
	if (_pcr_profile_a->get_profile_length() != _pcr_profile_b->get_profile_length())
	{
		err_msg << "ERROR: PCR_Profile::PCR_Profile(...) Profile Lenght of Profiles A and B is not the same" << endl;
		err_msg << "PCR Profiles must have same profile lenght" << endl;
		assert(NULL);
	}
	profile_length = _pcr_profile_a->get_profile_length();
	unsigned int max_number_of_primers = _pcr_profile_a->get_number_of_primers_primer_set() + _pcr_profile_b->get_number_of_primers_primer_set();
	p_set = new Primer_Set(max_number_of_primers, _pcr_profile_a->get_pointer_to_primer_set()->get_primer_length());
	assert(p_set);


	for (int i = 0; i < _pcr_profile_a->get_number_of_primers_primer_set(); i++)
	{
		p_set->add_primer(_pcr_profile_a->get_pointer_to_primer_set()->get_primer_as_value(i));
	}

	for (int i = 0; i < _pcr_profile_b->get_number_of_primers_primer_set(); i++)
	{
		p_set->add_primer(_pcr_profile_b->get_pointer_to_primer_set()->get_primer_as_value(i));
	}

	primer_profile_array_size = _pcr_profile_a->get_number_of_primers_location_profile() + _pcr_profile_b->get_number_of_primers_location_profile();
	number_of_primers = primer_profile_array_size;
	location_of_primer = new int[primer_profile_array_size];
	type_of_primer = new int[primer_profile_array_size];

	unsigned int counter_profile_a = 0;
	unsigned int counter_profile_b = 0;

	for (int i = 0; i < number_of_primers; i++)
	{
		if (_pcr_profile_a->get_location_of_primer()[counter_profile_a] < _pcr_profile_b->get_location_of_primer()[counter_profile_b])
		{
			location_of_primer[i] = _pcr_profile_a->get_location_of_primer()[counter_profile_a];
			type_of_primer[i] = _pcr_profile_a->get_type_of_primer()[counter_profile_a];
			counter_profile_a++;
			continue;
		}
		if (_pcr_profile_a->get_location_of_primer()[counter_profile_a] > _pcr_profile_b->get_location_of_primer()[counter_profile_b])
		{
			location_of_primer[i] = _pcr_profile_b->get_location_of_primer()[counter_profile_b];
			type_of_primer[i] = _pcr_profile_b->get_type_of_primer()[counter_profile_b];
			counter_profile_b++;
			continue;
		}
		if (_pcr_profile_a->get_location_of_primer()[counter_profile_a] == _pcr_profile_b->get_location_of_primer()[counter_profile_b])
		{
			location_of_primer[i] = _pcr_profile_a->get_location_of_primer()[counter_profile_b];
			
			if (_pcr_profile_a->get_type_of_primer()[counter_profile_a] == 5 || _pcr_profile_b->get_type_of_primer()[counter_profile_b] == 5)
			{
				type_of_primer[i] = 5;

			}
			else if (_pcr_profile_a->get_type_of_primer()[counter_profile_a] == _pcr_profile_b->get_type_of_primer()[counter_profile_b])
			{
				type_of_primer[i] = _pcr_profile_b->get_type_of_primer()[counter_profile_b];
			}
			else if (_pcr_profile_a->get_type_of_primer()[counter_profile_a] != _pcr_profile_b->get_type_of_primer()[counter_profile_b])
			{
				type_of_primer[i] = 5;
			}

			counter_profile_a++;
			counter_profile_b++;
			number_of_primers--;
		}
	}


	calculate_statistics();

}

PCR_Profile::PCR_Profile(PCR_Profile * _pcr_profile, ostream &err_msg)
{
	stats = _pcr_profile->get_Stats();

	profile_length = _pcr_profile->get_profile_length();
	unsigned int max_number_of_primers = _pcr_profile->get_number_of_primers_primer_set();
	p_set = new Primer_Set(max_number_of_primers, _pcr_profile->get_pointer_to_primer_set()->get_primer_length());

	for (int i = 0; i < _pcr_profile->get_number_of_primers_primer_set(); i++)
	{
		p_set->add_primer(_pcr_profile->get_pointer_to_primer_set()->get_primer_as_value(i));
	}

	primer_profile_array_size = _pcr_profile->get_number_of_primers_location_profile();
	number_of_primers = _pcr_profile->get_number_of_primers_location_profile();
	location_of_primer = new int[primer_profile_array_size];
	type_of_primer = new int[primer_profile_array_size];

	for (int i = 0; i < number_of_primers; i++)
	{
		location_of_primer[i] = _pcr_profile->get_location_of_primer()[i];
		type_of_primer[i] = _pcr_profile->get_type_of_primer()[i];
	}


}

PCR_Profile::PCR_Profile(Primer_Set * _primer_set, Sequence * _sequence, ostream &err_msg)
{
	stats = {};


	p_set = new Primer_Set(_primer_set, err_msg); assert(p_set);
	assert(_sequence->get_sequence_length() >= p_set->get_primer_length());
	profile_length = _sequence->get_sequence_length()-p_set->get_primer_length();

	primer_profile_array_size = 10000;
	number_of_primers = 0;
	location_of_primer = new int[primer_profile_array_size];
	type_of_primer = new int[primer_profile_array_size];
	for (int i = 0; i < primer_profile_array_size; i++)
	{
		location_of_primer[i] = -1;
		type_of_primer[i] = -1;
	}
	
	PCR_profile_calculation(_sequence);
	
	calculate_statistics();
}

bool PCR_Profile::add_primer_location_to_profile(unsigned int position, int primer_type, ostream &err_msg)
{
	for (int i = 0; i < number_of_primers; i++)
	{
		if (location_of_primer[i] == position)
		{
			if (type_of_primer[i] != primer_type)
			{
				type_of_primer[i] = 5;
			}
			return true;
		}
	}
	location_of_primer[number_of_primers] = position;
	type_of_primer[number_of_primers] = primer_type;
	number_of_primers++;
	return true;

}
bool PCR_Profile::show_statistics(ostream & out, ostream &err_msg)
{
	out << "number_forward_primers " << get_number_forward_primers() << endl;
	out << "number_reverse_primers " << get_number_reverse_primers() << endl;
	out << "number_short_amplicons " << get_number_short_amplicons() << endl;
	out << "number_long_amplicons " << get_number_long_amplicons() << endl;
	out << "total_lenght_short_amplicons " << get_total_lenght_short_amplicons() << endl;
	out << "total_lenght_long_amplicons " << get_total_lenght_long_amplicons() << endl;
	out << "total_lenght_too_long_amplicons " << get_total_lenght_too_long_amplicons() << endl;
	out << "total_length_uncovered" << get_total_length_uncovered() << endl;
	out << "total_length" << get_total_length_uncovered()+ get_total_lenght_too_long_amplicons() + get_total_lenght_long_amplicons() + get_total_lenght_short_amplicons() << endl;
	cout << "distance" << endl;
	for (int z = 0; z < number_of_primers; z++)cout << "Position " << location_of_primer[z] << "type " << type_of_primer[z] << endl;
	cout << "distance" << endl;
	return true;
}
bool PCR_Profile::show_All(ostream & out, ostream &err_msg)
{
	out << "Number of Primers " << p_set->get_number_of_primers() << endl;
	for (int i = 0; i < p_set->get_number_of_primers(); i++)
	{
		out << "Primers " << p_set->get_pointer_to_primer_array()[i] << endl;
	}
	show_statistics(out, err_msg);
	return true;
}

bool PCR_Profile::PCR_profile_calculation(Sequence * _sequence,ostream &err_msg)
{
	
	for (int i = 0; i < p_set->get_number_of_primers(); i++)
	{
		for (int j = 0; j < profile_length; j++)
		{
			if (number_of_primers >= primer_profile_array_size)
			{
				cout << "Realloc is being attempted made" << endl;
				primer_profile_array_size *= 1.5;
				location_of_primer = (int*)realloc(location_of_primer, sizeof(int)*primer_profile_array_size);
				type_of_primer = (int*)realloc(type_of_primer, sizeof(int)*primer_profile_array_size);
			}

			if (_sequence->get_pointer_to_sequence_int()[j] == p_set->get_pointer_to_primer_array()[i])
			{		
				add_primer_location_to_profile(j, 1); 
			}
			if (_sequence->get_pointer_to_sequence_int()[j] == p_set->get_pointer_to_reverse_complement_primer_array()[i])
			{
				add_primer_location_to_profile(j, -1); 
			}
		}
	}
	return true;
}


bool PCR_Profile::calculate_statistics(ostream & out, ostream &err_msg)
{
	int start_position = 0;
	int end_position = 0;
	unsigned int distance = 0;
	bool start_found = false;
	bool end_found = false;
	for (int i = 0; i < number_of_primers; i++)
	{
		if (start_found == true && end_found == true)
		{
			start_found = false;
			end_found = false;
		}
		if (!start_found && (type_of_primer[i] == 1 || type_of_primer[i] == 5))
		{
			start_found = true;
			end_found = false;
			start_position = location_of_primer[i];
			stats.number_forward_primers++;
		}
		if (start_found && !end_found && (type_of_primer[i] == -1 || type_of_primer[i] == 5))
		{
			start_found = false;
			end_found = true;

			end_position = location_of_primer[i];
			stats.number_reverse_primers++;
			distance = end_position - start_position + p_set->get_primer_length();
			if (distance > 10000000)
			{
				cout << "distance broke" << endl;
				for (int z = 0; z < number_of_primers; z++)cout << "Position " << location_of_primer[z] << "type " << type_of_primer[z] << endl;
				cout << "distance broke" << endl;
			}
			if (distance >= min_primer_distance && distance <= max_primer_distance)
			{
				stats.total_lenght_long_amplicons += distance;
				stats.number_long_amplicons++;
			}
			if (distance < min_primer_distance)
			{
				stats.total_lenght_short_amplicons += distance;
				stats.number_short_amplicons++;
			}
			if (distance>max_primer_distance) stats.total_lenght_too_long_amplicons += distance;

		}

	}
	return true;
}


#endif