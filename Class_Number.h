#define pi 3.14159265359
#define e 2.7182818284590452353602874

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include <gsl/gsl_sf_expint.h>

int kronecker(int n, int p);
int get_class_number(int d, long double R);
long double L(int D, long double R);
bool isSquare(int n);
void fill_units(std::string filename, int length, std::vector<long double> &arr);
void load_correct_answers(std::string filename, int length, std::vector<int> &arr);