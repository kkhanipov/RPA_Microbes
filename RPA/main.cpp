#include <iostream>
#include "Array_Sequences.h"
#include "Primer_Set.h"
#include "Optimization_Toolbox.h"
#include <time.h>
#include "PCR_Profile_Array.h"
#include <fstream>
#include <cstdlib>

using namespace std;

void test_1()
{
	double * x = new double[100];
	double * y = new double[100];
	bool * pareto = new bool[100];
	srand(time(0));
	for (int i = 0; i < 100; i++)
	{
		x[i] = (double)rand();
		y[i] = (double)rand();
		pareto[i] = 0;
	}
	
	Optimization_Toolbox::calculate_pareto_frontier(x, y, pareto, 100, false, true);
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
bool prepare_pareto(PCR_Profile ** pcr_profiles_4_comparison, unsigned int num_pcr_profiles, unsigned int & optimal_solution, ostream &err_msg=cout)
{
	double *genome_coverage, *genome_overrepresentation;
	bool *pareto_set;

	genome_coverage = new double[num_pcr_profiles];
	genome_overrepresentation = new double[num_pcr_profiles];
	pareto_set = new bool[num_pcr_profiles];

	for (int i = 0; i < num_pcr_profiles; i++)
	{
		genome_coverage[i] = pcr_profiles_4_comparison[i]->get_total_lenght_long_amplicons();
		genome_overrepresentation[i] = pcr_profiles_4_comparison[i]->get_total_lenght_short_amplicons();
		pareto_set[i] = false;
	}

	Optimization_Toolbox::calculate_pareto_frontier(genome_coverage, genome_overrepresentation, pareto_set, num_pcr_profiles, true, false, err_msg);

	int position = -1;
	
	int *tmp_pareto = new int[num_pcr_profiles];
	memset(tmp_pareto, 0, sizeof(int)*num_pcr_profiles);
	int n_pareto = 0;
	for (int i = 0; i < num_pcr_profiles; i++) if (pareto_set[i] == true) tmp_pareto[n_pareto++] = i;
	assert(n_pareto);
	position = tmp_pareto[rand() % n_pareto];
	


	optimal_solution = position;
	delete[]tmp_pareto;
	delete[] genome_coverage;
	delete[] genome_overrepresentation;
	delete[] pareto_set;
	return true;
}
int main()
{
	srand(time(0));
	char * file = "test";
	Array_Sequences * as;
	
	Primer_Set * primers = new Primer_Set("primers.fasta"); assert(primers);
	as = new Array_Sequences("sequence.fasta", primers->get_primer_length());
	//as = new Array_Sequences("synthetic.fasta",primers->get_primer_length()); 
	assert(as);

	as->show_Statistics();
	
//	as->show_All();

	unsigned int number_of_individual_primers = 4096;
	Primer_Set ** individual_primers;
	individual_primers = new Primer_Set *[number_of_individual_primers]; assert(individual_primers);

	PCR_Profile_Array * PCR_Profile_array;
	PCR_Profile_array = new PCR_Profile_Array(); assert(PCR_Profile_array);

	PCR_Profile ** individual_PCR_profiles;
	individual_PCR_profiles = new PCR_Profile *[number_of_individual_primers]; assert(individual_PCR_profiles);

	unsigned int count = 0;
	#pragma omp parallel for
	for (int i = 0; i < number_of_individual_primers; i++)
	{
		individual_primers[i] = new Primer_Set(1, 6); assert(individual_primers[i]);
		individual_primers[i]->add_primer(primers->get_primer_as_value(i));
		individual_PCR_profiles[i] = new PCR_Profile(individual_primers[i], as->get_pointer_to_sequence_object(0));
		assert(individual_PCR_profiles[i]);
		count++;
	}

	for (int i = 0; i < number_of_individual_primers; i++)
	{
		PCR_Profile_array->add_pcr_profile(individual_PCR_profiles[i]);
	}
	
	unsigned int pareto_index;
	prepare_pareto(PCR_Profile_array, pareto_index);
	

	//prepare_pareto(individual_PCR_profiles, number_of_individual_primers, pareto_index);
	PCR_Profile * pareto_PCR_profile = new PCR_Profile(PCR_Profile_array->get_pcr_profile(pareto_index)); assert(pareto_PCR_profile);
	pareto_PCR_profile->show_statistics();

	PCR_Profile_array->compress_left();
	//PCR_Profile_array->get_stats();

	count = 0;
	while (true)
	{
		system("PAUSE");
		PCR_Profile ** temp_pareto_PCR_profile;
		temp_pareto_PCR_profile = new PCR_Profile *[number_of_individual_primers]; assert(temp_pareto_PCR_profile);

#pragma omp parallel for
		for (int i = 0; i < number_of_individual_primers; i++)
		{
			temp_pareto_PCR_profile[i] = new PCR_Profile(pareto_PCR_profile, individual_PCR_profiles[i]); assert(temp_pareto_PCR_profile[i]);
		}
		count++;

		for (int i = 0; i < number_of_individual_primers; i++)
		{
			PCR_Profile_array->add_pcr_profile(temp_pareto_PCR_profile[i]);
		}
		if (!prepare_pareto(PCR_Profile_array, pareto_index)) break;

		//prepare_pareto(temp_pareto_PCR_profile, number_of_individual_primers, pareto_index);

		//if (temp_pareto_PCR_profile[pareto_index]->get_total_lenght_long_amplicons() > pareto_PCR_profile->get_total_lenght_long_amplicons() ||
		//	temp_pareto_PCR_profile[pareto_index]->get_total_lenght_short_amplicons() < pareto_PCR_profile->get_total_lenght_short_amplicons())
		//{
			delete pareto_PCR_profile;
			pareto_PCR_profile = new PCR_Profile(PCR_Profile_array->get_pcr_profile(pareto_index)); assert(pareto_PCR_profile);
			pareto_PCR_profile->show_statistics();
			PCR_Profile_array->compress_left();
			//PCR_Profile_array->get_stats();
		//}
		//else break;
		
		for (int i = 0; i<number_of_individual_primers; i++) delete temp_pareto_PCR_profile[i];
		delete[]temp_pareto_PCR_profile;
		
	}

	system("PAUSE");
	return 1;
}