#ifndef Array_Sequences_H
#define Array_Sequences_H

#include "Sequence.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <assert.h>


using namespace std;

//This class is meant to store a list of sequences read in from FASTA files. 
//The class uses class Sequences for each individual sequence, and this is an array that ties them all together.

class Array_Sequences
{
	Sequence ** sequence; //an array of pointers to sequence objects
	unsigned int number_of_sequences; // number of sequences located inside the array

public:

	static unsigned int max_sequence_length;
	static unsigned int max_number_of_sequences;

	Array_Sequences(char * filename, unsigned int primer_value=0, ostream &err_msg = cout); //constructor to be able to read a FASTA file and create the individual sequences and array of sequences
	Array_Sequences(ostream &err_msg = cout); //constructor to be able to read a FASTA file and create the individual sequences and array of sequences

	//Array_Sequences(Array_Sequences * _arr_seq, unsigned int primer_value = 0, ostream &err_msg = cout);
	~Array_Sequences();
	bool show_Statistics(ostream & out = cout, ostream &err_msg = cout); //show basic statistics about the array of sequences, number of sequences present,  and their nucleotide contribution

	bool show_All(ostream & out = cout, ostream &err_msg = cout); // shows the basic statistics about the array of sequences and then prints to screen the actual sequences

	bool add_sequence(Sequence * _sequence, unsigned int primer_value = 0, ostream &err_msg = cout); //ability to add a sequence to the array of sequences by sequence object
	bool add_sequence(char * _sequence, unsigned int _sequence_length, unsigned int primer_value = 0, ostream &err_msg = cout); //ability to add a sequence to the array of sequences by giving the actual sequence
	
	unsigned int get_number_of_sequences() { return number_of_sequences; }
	Sequence * get_pointer_to_sequence_object(unsigned int position) { return sequence[position]; }
};

unsigned int Array_Sequences::max_sequence_length = 1000000000;//TODO:make 4
unsigned int Array_Sequences::max_number_of_sequences = 5000;

Array_Sequences::~Array_Sequences()
{
	if (sequence != NULL)
	{
		for (int i = 0; i < number_of_sequences; i++) if(sequence[i]!=NULL) delete sequence[i];
		delete[] sequence;
	}
}
Array_Sequences::Array_Sequences(ostream &err_msg)
{
	number_of_sequences = 0;
	sequence = new Sequence*[max_number_of_sequences]; assert(sequence);
}

//NOTE: will not read NCBI fasta with multiple lines
Array_Sequences::Array_Sequences(char * filename, unsigned int _primer_value, ostream &err_msg)
{
	number_of_sequences = 0;
	// Check if input file is opened successfully
	if (filename == NULL)
	{
		err_msg << "ERROR: Array_Sequences::Array_Sequences(...) ==> filename == NULL" << endl; 
		err_msg << "No File Name Provided" << endl;
		assert(NULL);
	}

	sequence = new Sequence*[max_number_of_sequences]; 
	if (sequence == NULL)
	{
		err_msg << "ERROR: Array_Sequences::Array_Sequences(...) ==> Not enough memory (sequence == NULL)" << endl;
		assert(NULL);
	}

	ifstream in;
	in.open(filename);

	if (!in.is_open())
	{
		err_msg << "ERROR: Array_Sequences::Array_Sequences(...) ==> Cannot open FASTA file" << endl; 
		err_msg << "Failed to open file" << endl;
		assert(NULL);
	}

	char * _header_text = new char[10000 + 1];
	if (_header_text == NULL)
	{
		err_msg << "ERROR: Array_Sequences::Array_Sequences(...) ==> Not enough memory (_header_text == NULL)" << endl;
		err_msg << "Failed to read header" << endl;
		assert(NULL);
	}

	char * _sequence = new char[max_sequence_length + 1]; 
	if (_sequence == NULL)
	{
		err_msg << "ERROR: Array_Sequences::Array_Sequences(...) ==>  Not enough memory (_sequence == NULL)" << endl;
		assert(NULL);
	}

	// Read sequences in file
	while (!in.eof())
	{
		in.getline(_header_text, max_sequence_length);
		if (_header_text[0] == '\0' || _header_text[0] == EOF) break;

		in.getline(_sequence, max_sequence_length + 1);
		
		unsigned int tmp_sequence_length = strlen(_sequence);
		sequence[number_of_sequences] = new Sequence(_sequence, tmp_sequence_length, _primer_value); assert(sequence[number_of_sequences]);

		number_of_sequences++;
		//TODO: implement if necessary 
		/*Increase array size if needed
		if (number_of_sequences >= max_number_of_sequences)
		{
			//realloc stuff
		}
		*/
	}

	delete[] _header_text;
	delete[] _sequence;
}

bool Array_Sequences::add_sequence(char * _sequence, unsigned int _sequence_length, unsigned int _primer_value, ostream &err_msg)
{
	if (number_of_sequences >= max_number_of_sequences)
	{
		err_msg << "ERROR: Array_Sequences::add_sequence(...) ==> Not array size needs to be increased == NULL)" << endl;
		return false;
	}
	sequence[number_of_sequences] = new Sequence(_sequence, _sequence_length, _primer_value);
	assert(sequence[number_of_sequences]);
	number_of_sequences++;
	return true;
}

bool Array_Sequences::add_sequence(Sequence * _sequence, unsigned int _primer_value,ostream &err_msg)
{
	return add_sequence(_sequence->get_pointer_to_sequence(), _sequence->get_sequence_length(), _primer_value);
}

bool Array_Sequences::show_Statistics(ostream & out, ostream &err_msg)
{
	for (int i = 0; i < number_of_sequences; i++)
	{
		sequence[i]->show_statistics(out,err_msg);
	}
	return true;
}
bool Array_Sequences::show_All(ostream & out, ostream &err_msg)
{
	for (int i = 0; i < number_of_sequences; i++)
	{
		sequence[i]->show_All(out,err_msg);
	}
	return true;
}

#endif