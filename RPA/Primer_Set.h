#ifndef Primer_Set_H
#define Primer_Set_H



#include <iostream>
#include <cstring>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include "Array_Sequences.h"

using namespace std;
//this object will hold a list of primers to be used in the computation. It will be part of the PCR profile object. 

class Primer_Set
{
	static unsigned int max_number_of_primers;
	unsigned int primer_length; // the length of the primers

	unsigned int * primer; // contains a list of primers in their integer value
	unsigned int * reverse_complement;

	unsigned int number_of_primers; // the number of primers in the list

public:

	Primer_Set(unsigned int array_size_primers, unsigned int _primer_length, ostream & err_msg = cout);
	Primer_Set(char * filename, ostream & err_msg = cout); // ability to read a list or primers from a FASTA file
	Primer_Set(Primer_Set * primer_set, ostream & err_msg = cout);
	Primer_Set(unsigned int * primers, unsigned int number_primers, unsigned int _max_number_of_primers, ostream & err_msg = cout); //constructor to make a primer set from a array of primers converted to integers
	Primer_Set(char ** primers, unsigned int number_primers, unsigned int _max_number_of_primers, ostream & err_msg = cout);//constructor to make a primer set from a array of primers in char array
	~Primer_Set();
	bool write_to_file(char * filename, ostream & err_msg = cout); //writes the list of primers in FASTA format to file

	bool show_statistics(ostream & out = cout, ostream &err_msg = cout); //prints out the nucleotide contribution of the sequence

	bool show_All(ostream & out = cout, ostream &err_msg = cout); //prints out the show_statistics() and then the actual sequences

	unsigned int get_primer_as_value(unsigned int position) { return primer[position]; }

	//	bool get_primer_as_txt(unsigned int position, char *& primer, ostream & err_msg = cout); //copies the primers to the char ** primers object as text

	bool delete_primer(int *position, unsigned int number_of_primers_to_delete, ostream & err_msg = cout); //ability to mark a primer for deletion at a specific position in the array. 
//	bool delete_primer(int primer_value, ostream & err_msg = cout); //mark a primer for deletion based on its value
//	bool delete_primer(char * primer, ostream & err_msg = cout); //mark a primer for deletion based on its sequence

	bool add_primer(unsigned int primer_value, ostream & err_msg = cout); // add a primer to the array based on its integer value
	bool add_primer(char * primer, ostream & err_msg = cout); // add a primer to the array based on its sequence

	bool convert_primer_int_to_txt(unsigned int primer_value, char *& primer, ostream & err_msg = cout);
	bool convert_primer_txt_to_int(char * primer, unsigned int & primer_value, ostream & err_msg = cout);
	unsigned int convert_primer_to_reverse_complement(int position, ostream & err_msg = cout);
	unsigned int * get_pointer_to_primer_array() { return primer; };
	unsigned int * get_pointer_to_reverse_complement_primer_array() { return reverse_complement; };

	unsigned int get_number_of_primers() { return number_of_primers; };
	unsigned int get_primer_length() { return primer_length; };

};

unsigned int Primer_Set::max_number_of_primers = 5000;

Primer_Set::~Primer_Set()
{
	if (primer != NULL) delete[] primer;
	if (reverse_complement != NULL) delete[] reverse_complement;
}
Primer_Set::Primer_Set(unsigned int array_size_primers, unsigned int _primer_length, ostream & err_msg)
{
	assert(array_size_primers <= max_number_of_primers);

	primer_length = _primer_length;
	number_of_primers = 0;
	primer = new unsigned int[array_size_primers]; assert(primer);
	reverse_complement = new unsigned int[array_size_primers]; assert(reverse_complement);
}
Primer_Set::Primer_Set(Primer_Set * primer_set, ostream & err_msg)
{
	primer_length = primer_set->get_primer_length();
	number_of_primers = 0;
	primer = new unsigned int[max_number_of_primers]; assert(primer);
	reverse_complement = new unsigned int[max_number_of_primers]; assert(reverse_complement);
	for(int i=0;i<primer_set->get_number_of_primers();i++) add_primer(primer_set->get_primer_as_value(i));
}
bool Primer_Set::convert_primer_txt_to_int(char * _primer, unsigned int & primer_value, ostream & err_msg)
{
	primer_value = 0;
	char * _primer_value = new char[primer_length]; assert(_primer_value);
	for (int i = 0; i < primer_length; i++)
	{
		switch (_primer[i])
		{
		case 'A':_primer_value[i] = '1'; break;
		case 'T':_primer_value[i] = '2'; break;
		case 'C':_primer_value[i] = '3'; break;
		case 'G':_primer_value[i]=	'4'; break;
		default: strcpy(_primer_value, "555555"); primer_value = atoi(_primer_value);  return true;
		}
	}
	primer_value = atoi(_primer_value);
	return true;
}
Primer_Set::Primer_Set(char * filename, ostream & err_msg)
{
	Array_Sequences * as = new Array_Sequences(filename,0, err_msg); assert(as);
	if (as->get_number_of_sequences() == 0)
	{
		err_msg << "Primer_Set::Primer_Set(...) 0 primers were read from the file";
		err_msg << "could not read file with primers" << endl;
		assert(NULL);
	}

	number_of_primers = as->get_number_of_sequences();
	
	primer_length = as->get_pointer_to_sequence_object(0)->get_sequence_length();
	primer = new unsigned int[max_number_of_primers]; assert(primer);
	reverse_complement = new unsigned int[max_number_of_primers]; assert(reverse_complement);

	for (int i = 0; i < number_of_primers; i++)
	{
		unsigned int primer_val = 555555;
		convert_primer_txt_to_int(as->get_pointer_to_sequence_object(i)->get_pointer_to_sequence(), primer_val);
		primer[i] = primer_val;
		reverse_complement[i] = convert_primer_to_reverse_complement(i);
	}

	delete as;
}

bool Primer_Set::write_to_file(char * filename, ostream & err_msg)
{
	ofstream out;
	out.open(filename);
	if (!out.is_open())
	{
		err_msg << "ERROR: bool Primer_Set::write_to_file(), could not open " << filename << " file" << endl;
		assert(NULL);
	}
	if (number_of_primers == 0)
	{
		err_msg << "ERROR: bool Primer_Set::write_to_file() there are no primers to write" << endl;
		assert(NULL);
	}
	for (int i = 0; i < number_of_primers; i++)
	{
		out << ">" << i << endl;
		char * text_primer;
		if (!convert_primer_int_to_txt(primer[i], text_primer))
		{
			err_msg << "ERROR: bool Primer_Set::write_to_file() could not convert primer " << i << " to text" << endl;
			assert(NULL);
		}
		out << text_primer << endl;
	}
	out.close();
	return true;
}

bool Primer_Set::show_statistics(ostream & out, ostream &err_msg)
{
	int count_A = 0;
	int count_C = 0;
	int count_T = 0;
	int count_G = 0;
	int count_N = 0;

	out << "Primer length is " << primer_length << endl;
	out << "Number of primers is " << number_of_primers << endl;
	for (int i = 0; i < number_of_primers; i++)
	{
		char * primer_seq = NULL;
		
		if (!convert_primer_int_to_txt(primer[i], primer_seq))
		{
			err_msg << "ERROR: bool Primer_Set::show_statistics() could not convert primer " << i << " to text" << endl;
			assert(NULL);
		}
		for (int i = 0; i < primer_length; i++)
		{
			switch (primer_seq[i])
			{
			case 'A':case 'a': count_A++; break;
			case 'T':case 't': count_T++; break;
			case 'C':case 'c': count_C++; break;
			case 'G':case 'g': count_G++; break;
			case 'N':case 'n': count_N++; break;
			default:
				err_msg << "ERROR: Sequence::Sequence ==> unexpected character: sequence[" << i << "]=" << primer_seq[i] << endl;
				assert(NULL);
			}
		}

		delete[] primer_seq;
	}
	out << "Nucleotide Contribution is:" << endl;
	out << "A :" << count_A << endl;
	out << "C :" << count_C << endl;
	out << "T :" << count_T << endl;
	out << "G :" << count_G << endl;
	out << "N :" << count_N << endl;

	return true;
}

bool Primer_Set::convert_primer_int_to_txt(unsigned int primer_value, char *& _primer, ostream & err_msg)
{
	_primer = new char[primer_length + 1]; assert(_primer);
	char * conversion_primer = new char[primer_length]; assert(conversion_primer);
	sprintf(conversion_primer, "%d", primer_value);
	for (int i = 0; i < primer_length; i++)
	{
		switch (conversion_primer[i])
		{
		case '1':_primer[i] = 'A'; break;
		case '2':_primer[i] = 'T'; break;
		case '3':_primer[i] = 'C'; break;
		case '4':_primer[i] = 'G'; break;
		default: _primer[i] = 'N';
		}
	}
	_primer[primer_length] = '\0';
	
	return true;
}
unsigned int Primer_Set::convert_primer_to_reverse_complement(int position, ostream & err_msg)
{
	char * conversion_primer = new char[primer_length]; assert(conversion_primer);
	char * reverse_complement_primer = new char[primer_length]; assert(reverse_complement_primer);
	sprintf(conversion_primer, "%d", primer[position]);
	int ii = 0;
	for (int i = primer_length-1; i >= 0; i--)
	{
		
		switch (conversion_primer[i])
		{
		case '1':reverse_complement_primer[ii] = '2'; break;
		case '2':reverse_complement_primer[ii] = '1'; break;
		case '3':reverse_complement_primer[ii] = '4'; break;
		case '4':reverse_complement_primer[ii] = '3'; break;
		default: reverse_complement_primer[ii] = '5';
		}
		ii++;
	}
	return atoi(reverse_complement_primer);
}
bool Primer_Set::show_All(ostream & out, ostream &err_msg)
{
	show_statistics(out,err_msg);

	for (int i = 0; i < number_of_primers; i++)
	{
		char * primer_seq = NULL;
		if (!convert_primer_int_to_txt(primer[i], primer_seq))
		{
			err_msg << "ERROR: bool Primer_Set::show_statistics() could not convert primer " << i << " to text" << endl;
			assert(NULL);
		}
		out << "Primer id: " << i << " Primer Value :" << primer[i] << " text: ";
		for (int j = 0; j < primer_length; j++)out << primer_seq[j];
		out << endl;

		delete[] primer_seq;
	}
	return true;
}
bool Primer_Set::add_primer(unsigned int primer_value, ostream & err_msg)
{

	if (number_of_primers >= max_number_of_primers)
	{
		err_msg << "ERROR: bool Primer_Set::add_primer() number_of_primers >= max_number_of_primers" << endl;
		assert(NULL);
	}
	for (int i = 0; i < number_of_primers; i++)
	{
		if (primer[i] == primer_value)
		{
			//err_msg << "WARNING: bool Primer_Set::add_primer() primer already included in the list" << endl;
			return true;
		}
	}
	primer[number_of_primers] = primer_value;
	reverse_complement[number_of_primers] = convert_primer_to_reverse_complement(number_of_primers);//TODO: make argument unsigned int
	number_of_primers++;
	return true;
}
bool Primer_Set::add_primer(char * _primer, ostream & err_msg)
{
	unsigned int primer_value = 0;
	convert_primer_txt_to_int(_primer, primer_value);
	add_primer(primer_value, err_msg);
	return true;
}
bool Primer_Set::delete_primer(int *position, unsigned int number_of_primers_to_delete, ostream & err_msg)
{
	assert(number_of_primers >= 1);

	unsigned int * _primer = new unsigned[max_number_of_primers]; assert(_primer);
	unsigned int * _reverse_complement = new unsigned[max_number_of_primers]; assert(_reverse_complement);
	int ii = 0;
	bool found_primer;
	for (int i = 0; i < number_of_primers; i++)
	{
		found_primer = false;
		for (int j = 0; j < number_of_primers_to_delete; j++)
		{
			if (i == position[j])
			{
				found_primer = true;
				break;
			}
		}
		if (!found_primer)
		{
			_primer[ii] = primer[i];
			_reverse_complement[ii] = reverse_complement[i];
		}
		ii++;
	}
	number_of_primers--;
	delete[] primer;
	delete[] reverse_complement;
	primer = _primer;
	reverse_complement = _reverse_complement;
	return true;
}



#endif // !Primer_Set_H