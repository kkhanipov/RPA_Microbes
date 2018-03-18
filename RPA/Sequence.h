#ifndef Sequence_H
#define Sequence_H

#include <iostream>
#include <assert.h>
#include "Optimization_Toolbox.h"

using namespace std;

// this class is meant to store sequences read from a FASTA file
class Sequence
{
	char * dna_sequence;
	unsigned int * int_dna_sequence;
	unsigned int seq_length;
	unsigned int usable_length;

public:
	
	Sequence(char * sequence, unsigned int sequence_length, unsigned int primer_length=0,ostream & err_msg=cout);
	Sequence(Sequence *_sequence, unsigned int primer_length = 0, ostream & err_msg = cout);
	~Sequence();
	bool show_statistics(ostream & out=cout, ostream &err_msg=cout); //prints out the nucleotide contribution of the sequence
	bool show_All(ostream & out = cout, ostream &err_msg = cout); //prints out the show_statistics() and then the actual sequence

	char * get_pointer_to_sequence() {return dna_sequence;}
	unsigned int * get_pointer_to_sequence_int() { return int_dna_sequence; }
	int get_sequence_length() { return seq_length; }
	unsigned int get_usable_length() { return usable_length; }
};
Sequence::~Sequence()
{
	if (dna_sequence != NULL)
	{
		delete[] dna_sequence;
	}
	if (int_dna_sequence != NULL)
	{
		delete[] int_dna_sequence;
	}
}
Sequence::Sequence(Sequence *_sequence, unsigned int _primer_length, ostream & err_msg)
{
	usable_length = 0;
	if (_sequence->get_sequence_length() == 0)
	{
		err_msg << "ERROR: Sequence::Sequence ==> _sequence_length == 0" << endl; 
		err_msg << "Sequence could not be allocated b/c length is 0" << endl;
		assert(NULL);
	}
	dna_sequence = new char[_sequence->get_sequence_length()+1];
	seq_length = _sequence->get_sequence_length();
	if (dna_sequence == NULL)
	{
		err_msg << "ERROR: Sequence::Sequence ==> sequence == NULL" << endl; 
		err_msg << "Sequence could not be allocated" << endl;
		assert(NULL);
	}
	for (int i = 0; i < _sequence->get_sequence_length(); i++)
	{
		dna_sequence[i] = _sequence->get_pointer_to_sequence()[i];
	}
	dna_sequence[seq_length] = '\0';
	if (_primer_length == 0)int_dna_sequence = NULL;
	else
	{
		int_dna_sequence = new unsigned int[seq_length - _primer_length];
		unsigned int primer_length = _primer_length;
		unsigned int primer_value = 555555;
		for (int i = 0; i < seq_length; i++)
		{
			primer_value = 555555;
			Optimization_Toolbox::convert_primer_txt_to_int(&dna_sequence[i], primer_length, primer_value);
			int_dna_sequence[i] = primer_value;
			if (primer_value != 555555)usable_length++;
		}
	}
	
}

Sequence::Sequence(char * _sequence, unsigned int _sequence_length, unsigned int _primer_length, ostream & err_msg)
{
	usable_length = 0;
	if (_sequence_length == 0)
	{
		err_msg << "ERROR: Sequence::Sequence ==> _sequence_length == 0" << endl; 
		err_msg << "Sequence could not be allocated b/c length is 0" << endl;
		assert(NULL);
	}
	seq_length = _sequence_length;
	dna_sequence = new char[_sequence_length + 1];
	if (dna_sequence == NULL)
	{
		err_msg << "ERROR: Sequence::Sequence ==> sequence == NULL" << endl; 
		err_msg << "Sequence could not be allocated" << endl;
		assert(NULL);
	}
	for (int i = 0; i < _sequence_length; i++)
	{
		dna_sequence[i] = _sequence[i];
	}
	dna_sequence[seq_length] = '\0';
	if (_primer_length == 0)int_dna_sequence = NULL;
	else
	{
		int_dna_sequence = new unsigned int[seq_length - _primer_length];
		unsigned int primer_length = _primer_length;
		unsigned int primer_value = 555555;
		for (int i = 0; i < seq_length; i++)
		{
			primer_value = 555555;
			Optimization_Toolbox::convert_primer_txt_to_int(&dna_sequence[i], primer_length, primer_value);
			int_dna_sequence[i] = primer_value;
			if (primer_value != 555555)usable_length++;
		}
	}
}

bool Sequence::show_statistics(ostream & out, ostream &err_msg)
{
	int count_A = 0;
	int count_C = 0;
	int count_T = 0;
	int count_G = 0;
	int count_N = 0;
	for (int i = 0; i < seq_length; i++)
	{
		switch (dna_sequence[i])
		{
		case 'A':case 'a': count_A++; break;
		case 'T':case 't': count_T++; break;
		case 'C':case 'c': count_C++; break;
		case 'G':case 'g': count_G++; break;
		case 'N':case 'n': count_N++; break;
		default:
			err_msg << "ERROR: Sequence::Sequence ==> unexpected character: sequence[" << i << "]=" << dna_sequence[i] << endl;
		}
	}

	out << "Nucleotide Contribution is:" << endl;
	out << "A :"<< count_A << endl;
	out << "C :" << count_C << endl;
	out << "T :" << count_T << endl;
	out << "G :" << count_G << endl;
	out << "N :" << count_N << endl;

	return true;
}

bool Sequence::show_All(ostream & out, ostream &err_msg)
{
	show_statistics();

	for (int i = 0; i < seq_length; i++)
	{
		out << dna_sequence[i] << "\t";
		out << int_dna_sequence[i] << "\t";
	}
	out << endl;

	return true;
}


#endif