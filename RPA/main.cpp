#include <iostream>
#include "Array_Sequences.h"
#include "Primer_Set.h"
#include "Optimization_Toolbox.h"
#include <time.h>
#include "PCR_Profile_Array.h"
#include <fstream>
#include <cstdlib>

using namespace std;

bool prepare_pareto(PCR_Profile_Array * pcr_profile_array_background, PCR_Profile_Array * pcr_profile_array_target, unsigned int & optimal_solution, ostream &err_msg = cout)
{
	double *target_genome_coverage, *background_genome_coverage;
	bool *pareto_set;
	unsigned int number_of_profiles;

	number_of_profiles = pcr_profile_array_background->get_number_of_profiles();

	target_genome_coverage = new double[number_of_profiles];
	background_genome_coverage = new double[number_of_profiles];
	pareto_set = new bool[number_of_profiles];

	for (int i = 0; i < number_of_profiles; i++)
	{
		background_genome_coverage[i] = pcr_profile_array_background->get_pcr_profile(i)->get_total_lenght_long_amplicons();
		target_genome_coverage[i] = pcr_profile_array_target->get_pcr_profile(i)->get_total_lenght_long_amplicons();
		pareto_set[i] = false;
	}

	Optimization_Toolbox::calculate_pareto_frontier(target_genome_coverage, background_genome_coverage, pareto_set, number_of_profiles, true, false, err_msg);
	int *tmp_pareto = new int[number_of_profiles];
	memset(tmp_pareto, 0, sizeof(int)*number_of_profiles);
	int n_pareto = 0;

	int position = -1;
	bool found_new_pareto = false;
	int count = 0;
	for (int i = 0; i < number_of_profiles; i++)
	{

		if (pareto_set[i] == true && pcr_profile_array_target->get_to_be_considered()[i] == 0 && pcr_profile_array_background->get_pcr_profile(i)->get_total_lenght_long_amplicons() <2335499)
		{

			tmp_pareto[n_pareto] = i;
			n_pareto++;
			found_new_pareto = true;
			pcr_profile_array_target->mark_profile_pareto(i);
			pcr_profile_array_background->mark_profile_pareto(i);
		}

		if (pareto_set[i] == false)
		{
			pcr_profile_array_target->delete_pcr_profile(i);
			pcr_profile_array_background->delete_pcr_profile(i);

			count++;
		}
	}
	//cout << count << endl;
	if (!found_new_pareto) return false;
	assert(n_pareto);

	position = tmp_pareto[rand() % n_pareto];
	optimal_solution = position;
	delete[]tmp_pareto;
	delete[] target_genome_coverage;
	delete[] background_genome_coverage;
	delete[] pareto_set;
	return true;
	return true;
}
bool prepare_pareto(PCR_Profile_Array * pcr_profile_array_4_comparison, unsigned int & optimal_solution, ostream &err_msg = cout)
{
	double *genome_coverage, *genome_overrepresentation;
	bool *pareto_set;
	unsigned int number_of_profiles;

	number_of_profiles=pcr_profile_array_4_comparison->get_number_of_profiles();

	genome_coverage = new double[number_of_profiles];
	genome_overrepresentation = new double[number_of_profiles];
	pareto_set = new bool[number_of_profiles];

	for (int i = 0; i < number_of_profiles; i++)
	{
		genome_coverage[i]= pcr_profile_array_4_comparison->get_pcr_profile(i)->get_total_lenght_long_amplicons();
		genome_overrepresentation[i] = pcr_profile_array_4_comparison->get_pcr_profile(i)->get_total_lenght_short_amplicons();
		pareto_set[i] = false;
	}

	Optimization_Toolbox::calculate_pareto_frontier(genome_coverage, genome_overrepresentation, pareto_set, number_of_profiles, true, false, err_msg);
	int *tmp_pareto = new int[number_of_profiles];
	memset(tmp_pareto, 0, sizeof(int)*number_of_profiles);
	int n_pareto = 0;

	int position = -1;
	bool found_new_pareto = false;
	int count = 0;
	for (int i = 0; i < number_of_profiles; i++)
	{

		if (pareto_set[i] == true && pcr_profile_array_4_comparison->get_to_be_considered()[i] == 0)
		{

			tmp_pareto[n_pareto] = i;
			n_pareto++;
			found_new_pareto = true;
			pcr_profile_array_4_comparison->mark_profile_pareto(i);
		}

		if (pareto_set[i] == false)
		{
			pcr_profile_array_4_comparison->delete_pcr_profile(i);
			count++;
		}
	}
	//cout << count << endl;
	if (!found_new_pareto) return false;
	assert(n_pareto);

	position = tmp_pareto[rand() % n_pareto];
	optimal_solution = position;
	delete[]tmp_pareto;
	delete[] genome_coverage;
	delete[] genome_overrepresentation;
	delete[] pareto_set;
	return true;
}

int main(int argc, char **argv)
{
	srand(time(0));
	if (argc != 4)
	{
		cout << "<Primers> <Background> <Target>" << endl;
		return -1;
	}
	Array_Sequences * background;
	Array_Sequences * target;
	Primer_Set * primers = new Primer_Set(argv[1]); assert(primers);
	background = new Array_Sequences(argv[2], primers->get_primer_length());	assert(background);
	target = new Array_Sequences(argv[3], primers->get_primer_length()); assert(target);
	//as = new Array_Sequences("synthetic.fasta",primers->get_primer_length()); 


	background->show_Statistics();
	target->show_Statistics();
	
//	as->show_All();

	unsigned int number_of_individual_primers = primers->get_number_of_primers();
	Primer_Set ** individual_primers;
	individual_primers = new Primer_Set *[number_of_individual_primers]; assert(individual_primers);


	PCR_Profile_Array * PCR_Profile_array_background;
	PCR_Profile_array_background = new PCR_Profile_Array(); assert(PCR_Profile_array_background);

	PCR_Profile_Array * PCR_Profile_array_target;
	PCR_Profile_array_target = new PCR_Profile_Array(); assert(PCR_Profile_array_target);

	PCR_Profile ** individual_PCR_profiles_background;
	individual_PCR_profiles_background = new PCR_Profile *[number_of_individual_primers]; assert(individual_PCR_profiles_background);

	PCR_Profile ** individual_PCR_profiles_target;
	individual_PCR_profiles_target = new PCR_Profile *[number_of_individual_primers]; assert(individual_PCR_profiles_target);

	unsigned int count = 0;
	#pragma omp parallel for
	for (int i = 0; i < number_of_individual_primers; i++)
	{
		individual_primers[i] = new Primer_Set(1, 6); assert(individual_primers[i]);
		individual_primers[i]->add_primer(primers->get_primer_as_value(i));

		individual_PCR_profiles_background[i] = new PCR_Profile(individual_primers[i], background->get_pointer_to_sequence_object(0));
		assert(individual_PCR_profiles_background[i]);

		individual_PCR_profiles_target[i] = new PCR_Profile(individual_primers[i], target->get_pointer_to_sequence_object(0));
		assert(individual_PCR_profiles_target[i]);
		count++;
	}

	for (int i = 0; i < number_of_individual_primers; i++)
	{
		PCR_Profile_array_background->add_pcr_profile(individual_PCR_profiles_background[i]);
		PCR_Profile_array_target->add_pcr_profile(individual_PCR_profiles_target[i]);
	}
	
	unsigned int pareto_index;
	prepare_pareto(PCR_Profile_array_background, PCR_Profile_array_target, pareto_index);
	

	//prepare_pareto(individual_PCR_profiles, number_of_individual_primers, pareto_index);
	PCR_Profile * pareto_PCR_profile_background = new PCR_Profile(PCR_Profile_array_background->get_pcr_profile(pareto_index)); assert(pareto_PCR_profile_background);
	PCR_Profile * pareto_PCR_profile_target = new PCR_Profile(PCR_Profile_array_target->get_pcr_profile(pareto_index)); assert(pareto_PCR_profile_target);
	
	cout << "Background % Amplified" << "\t" << "Target % Amplified" << endl;
	cout << (double)pareto_PCR_profile_background->get_total_lenght_long_amplicons() / (double)background->get_pointer_to_sequence_object(0)->get_usable_length() << "\t" << (double)pareto_PCR_profile_target->get_total_lenght_long_amplicons() / (double)target->get_pointer_to_sequence_object(0)->get_usable_length() << endl;
	PCR_Profile_array_background->compress_left();
	PCR_Profile_array_target->compress_left();


	count = 0;
	while (true)
	{

		PCR_Profile ** temp_pareto_PCR_profile_background;
		temp_pareto_PCR_profile_background = new PCR_Profile *[number_of_individual_primers]; assert(temp_pareto_PCR_profile_background);


		PCR_Profile ** temp_pareto_PCR_profile_target;
		temp_pareto_PCR_profile_target = new PCR_Profile *[number_of_individual_primers]; assert(temp_pareto_PCR_profile_target);

#pragma omp parallel for
		for (int i = 0; i < number_of_individual_primers; i++)
		{
			temp_pareto_PCR_profile_background[i] = new PCR_Profile(pareto_PCR_profile_background, individual_PCR_profiles_background[i]); assert(temp_pareto_PCR_profile_background[i]);
			temp_pareto_PCR_profile_target[i] = new PCR_Profile(pareto_PCR_profile_target, individual_PCR_profiles_target[i]); assert(temp_pareto_PCR_profile_target[i]);

		}
		count++;

		for (int i = 0; i < number_of_individual_primers; i++)
		{
			PCR_Profile_array_background->add_pcr_profile(temp_pareto_PCR_profile_background[i]);
			PCR_Profile_array_target->add_pcr_profile(temp_pareto_PCR_profile_target[i]);

		}
		if (!prepare_pareto(PCR_Profile_array_background, PCR_Profile_array_target, pareto_index)) break;



		if (
			PCR_Profile_array_background->get_pcr_profile(pareto_index)->get_total_lenght_long_amplicons() == pareto_PCR_profile_background->get_total_lenght_long_amplicons()
			&&
			PCR_Profile_array_background->get_pcr_profile(pareto_index)->get_total_lenght_short_amplicons() == pareto_PCR_profile_background->get_total_lenght_short_amplicons()
			&&
			PCR_Profile_array_target->get_pcr_profile(pareto_index)->get_total_lenght_long_amplicons() == pareto_PCR_profile_target->get_total_lenght_long_amplicons()
			&&
			PCR_Profile_array_target->get_pcr_profile(pareto_index)->get_total_lenght_short_amplicons() == pareto_PCR_profile_target->get_total_lenght_short_amplicons()

			) break;
			delete pareto_PCR_profile_background; delete pareto_PCR_profile_target;
			pareto_PCR_profile_background = new PCR_Profile(PCR_Profile_array_background->get_pcr_profile(pareto_index)); assert(pareto_PCR_profile_background);
			pareto_PCR_profile_target = new PCR_Profile(PCR_Profile_array_target->get_pcr_profile(pareto_index)); assert(pareto_PCR_profile_target);

			cout << (double)pareto_PCR_profile_background->get_total_lenght_long_amplicons() / (double)background->get_pointer_to_sequence_object(0)->get_usable_length() << "\t" << (double)pareto_PCR_profile_target->get_total_lenght_long_amplicons() / (double)target->get_pointer_to_sequence_object(0)->get_usable_length() << endl;

			PCR_Profile_array_background->compress_left();
			PCR_Profile_array_target->compress_left();
			//PCR_Profile_array->get_stats();
		//}
		//else break;
		
			for (int i = 0; i < number_of_individual_primers; i++)
			{
				delete temp_pareto_PCR_profile_background[i];
				delete temp_pareto_PCR_profile_target[i];
			}
		delete[]temp_pareto_PCR_profile_background;
		delete[]temp_pareto_PCR_profile_target;

		
	}

	return 1;
}