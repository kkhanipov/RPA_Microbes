#ifndef Optimization_Toolbox_H
#define Optimization_Toolbox_H

#include <iostream>
#include <assert.h>
#include <time.h>
#include <limits>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

class Optimization_Toolbox
{
public:
	static bool calculate_pareto_frontier(double * x, double * y, bool * paretto_set, unsigned int number_of_values, bool maximize_x, bool maximize_y, ostream &err_msg = cout);
	
	//the function will calculate the parretto frontier of a dataset. It will take an array of x and y coordinates. Maximize_x and Maximize_y are to denote whether the values need to maximized or minimized
	//Use quicksort on them and pick the best point.
	//From there it will go left and right to determine all of the points which belong in the parretto frontier. The function will return the paretto_set marking positions in which points are part of the 
	//paretto frontier as "true".
	static bool convert_primer_txt_to_int(char * _primer, unsigned int primer_length, unsigned int & primer_value, ostream & err_msg = cout);

};
bool Optimization_Toolbox::convert_primer_txt_to_int(char * _primer, unsigned int primer_length, unsigned int & primer_value, ostream & err_msg)
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
		case 'G':_primer_value[i] = '4'; break;
		default: strcpy(_primer_value, "555555"); primer_value = atoi(_primer_value);  return true;
		}
	}
	primer_value = atoi(_primer_value);
	return true;
}
struct Sortable_Pareto
{
	unsigned int index;
	double x, y;
	bool pareto_set;
};

int compare_Sortable_Pareto(const void *_a,const void *_b)
{
	Sortable_Pareto a = *((Sortable_Pareto*)_a);
	Sortable_Pareto b = *((Sortable_Pareto*)_b);
	if (a.x > b.x) return -1;
	if (a.x == b.x && a.y<b.y) return -1;
	else return 1;
}
int compare_Sortable_Pareto_Max_X_Min_Y(const void *_a, const void *_b)
{
	Sortable_Pareto *a = (Sortable_Pareto*)_a;
	Sortable_Pareto *b = (Sortable_Pareto*)_b;
	if (a->x > b->x) return -1;
	if (a->x == b->x && a->y < b->y) return -1;
	return 1;
}
int compare_Sortable_Pareto_Max_X_Max_Y(const void *_a, const void *_b)
{
	Sortable_Pareto *a = (Sortable_Pareto*)_a;
	Sortable_Pareto *b = (Sortable_Pareto*)_b;
	if (a->x > b->x) return -1;
	if (a->x == b->x && a->y > b->y) return -1;
	return 1;
}
int compare_Sortable_Pareto_Min_X_Min_Y(const void *_a, const void *_b)
{
	Sortable_Pareto *a = (Sortable_Pareto*)_a;
	Sortable_Pareto *b = (Sortable_Pareto*)_b;
	if (a->x < b->x) return -1;
	if (a->x == b->x && a->y < b->y) return -1;
	return 1;
}
int compare_Sortable_Pareto_Min_X_Max_Y(const void *_a, const void *_b)
{
	Sortable_Pareto *a = (Sortable_Pareto*)_a;
	Sortable_Pareto *b = (Sortable_Pareto*)_b;
	if (a->y > b->y) return -1;
	if (a->y == b->y && a->x < b->x) return -1;
	return 1;
}

bool Optimization_Toolbox::calculate_pareto_frontier(double * x, double * y, bool * pareto_set, unsigned int number_of_values, bool maximize_x, bool maximize_y, ostream &err_msg)
{
	Sortable_Pareto *sorted_set = new Sortable_Pareto[number_of_values]; assert(sorted_set);
	double max_x = std::numeric_limits<double>::min();
	double max_y = std::numeric_limits<double>::min();
	double min_x = std::numeric_limits<double>::max();
	double min_y = std::numeric_limits<double>::max();

	for (unsigned int i = 0; i < number_of_values; i++)
	{
		sorted_set[i].index = i;
		sorted_set[i].x = x[i];
		sorted_set[i].y = y[i];
		sorted_set[i].pareto_set = false;
	}

	if (maximize_x && !maximize_y) qsort(sorted_set, number_of_values, sizeof(Sortable_Pareto), compare_Sortable_Pareto_Max_X_Min_Y);
	if (maximize_x && maximize_y) qsort(sorted_set, number_of_values, sizeof(Sortable_Pareto), compare_Sortable_Pareto_Max_X_Max_Y);
	if (!maximize_x && !maximize_y) qsort(sorted_set, number_of_values, sizeof(Sortable_Pareto), compare_Sortable_Pareto_Min_X_Min_Y);
	if (!maximize_x && maximize_y) qsort(sorted_set, number_of_values, sizeof(Sortable_Pareto), compare_Sortable_Pareto_Min_X_Max_Y);

	
	for (unsigned int i = 0; i < number_of_values; i++)
	{
		if (maximize_x && !maximize_y)
		{
			if (max_x < sorted_set[i].x)
			{
				max_x = sorted_set[i].x;
				sorted_set[i].pareto_set = true;
			}
			if (min_y > sorted_set[i].y)
			{
				min_y = sorted_set[i].y;
				sorted_set[i].pareto_set = true;
			}

		}

		//TODO: missing implementation of maxx maxy
		//if (maximize_x && maximize_y)
		/*
		if (sorted_set[i].pareto_set)
		{
			err_msg << "index :" << sorted_set[i].index << "\t";
			err_msg << "x :" << sorted_set[i].x << "\t";
			err_msg << "y :" << sorted_set[i].y << "\t";
			err_msg << "x/y :" << sorted_set[i].x / sorted_set[i].y << "\t";
			err_msg << "pareto set :" << sorted_set[i].pareto_set << endl;
		}
		*/
	}

	for (unsigned int i = 0; i < number_of_values; i++)
	{
		if (sorted_set[i].pareto_set == true) pareto_set[sorted_set[i].index] = true;
	}

	return true;
}

#endif