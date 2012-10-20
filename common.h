#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>

#include "structs.h"

using namespace std;

void header();

string capslock (string input);

const string int_to_string(int i);
const string char_to_string (char i);

double SIGNUM (double in);

double SIN  (double in);
double COS 	(double in);
double ASIN (double in);
double ACOS (double in);
double ATAN (double in);

int rounding (double in);
int mm_to_point (int i);

VCTR  crossproduct (VCTR in1, VCTR in2);
double dotproduct (VCTR in1, VCTR in2, bool normalisation);

double det_3 (vector <vector <double> > in);

vector <vector <double> >  declare_3x3_matrix (double a, double b, double c, double d, double e, double f, double g, double h, double i);

vector <vector <double> > init_matrix (const int dimension);
vector <vector <double> > init_matrix (const int i, const int j);

vector <vector <double> > identity_matrix (vector <vector <double> >);

vector <double>  init_vector (int dimension);
int search_max_off_diagonal_element_in_mtrx (vector <vector <double> > in);
int m_from_max_element (int max_element, size_t m, size_t n);
double teta (vector <vector <double> > in, int m, int n);
vector <vector <double> > init_rotation_mtrx (double teta, int m, int n, int dimension);
vector <vector <double> > outer_product (vector <double> in);
vector <vector <double> > transpose (vector <vector <double> > in);
vector <vector <double> > mult_mtrx (vector <vector <double> > in1, vector <vector <double> > in2);
vector <vector <double> > add_mtrx (vector <vector <double> > in1, vector <vector <double> > in2);

size_t return_second_eigenvalue (vector <vector< double > > in);
size_t return_first_eigenvalue (vector <vector< double > > in);

vector <vector <double> > jacobi (vector <vector <double> > in);

vector <vector <double> > gaussian_elimination (vector <vector <double> > in);

vector <vector <double> > LU_decomposition (vector <vector <double> > in);
vector <vector <double> > compute_Z (vector <vector <double> > L, vector <vector <double> > c);
vector <vector <double> > compute_X (vector <vector <double> > U, vector <vector <double> > Z);


vector <vector <double> > row_division_diagonal (vector <vector <double> > in,  double rownumber, double value);
vector <vector <double> > row_addition (vector <vector <double> > in,  double actual_row_number, double zero_row_number, double value);

vector <vector <double> > row_addition_LU (vector <vector <double> > in,  double actual_row_number, double zero_row_number, double value);

vector < vector < double > > generate_A (vector < vector < double > > EVEV);
vector < vector < double > > generate_D (vector < vector < double > > EVEV);

bool check_fry_matrix (size_t second_eigenvalue, vector <vector <double> > in_eigenvector);
bool check_correct_stressfield (STRESSFIELD sf);

VCTR generarte_stress_colors (double value);

VCTR unitvector (VCTR in);
CENTR_VECT unitvector (CENTR_VECT in);
vector < double > unitvector (vector < double > in);

VCTR declare_vector (double a, double b, double c);
CENTR_VECT declare_vector (double a, double b, double c, double d, double e, double f);

VCTR flip_D_vector (VCTR in);
VCTR flip_N_vector (VCTR in);
VCTR flip_ptn_vector (VCTR in);

VCTR compute_d_for_SC (GDB_ i);

VCTR DXDYDZ_from_dipdir_dip (DIPDIR_DIP i);
VCTR NXNYNZ_from_dipdir_dip (DIPDIR_DIP i);

VCTR inversion_DXDYDZ_from_DXDYDZ (VCTR i);
STRESSFIELD stressvector_to_DXDYDZ (STRESSFIELD in);

VCTR DXDYDZ_from_NXNYNZ (VCTR i);
VCTR NXNYNZ_from_DXDYDZ (VCTR i);

DIPDIR_DIP dipdir_dip_from_DXDYDZ (VCTR i);
DIPDIR_DIP dipdir_dip_from_NXNYNZ (VCTR i);

VCTR ROTATE (VCTR axis, VCTR torotate, double angle);

bool existence (string expression, vector<GDB_> inGDB);
bool existence_of_group (int expression, vector <int> whichgroup);
bool existence_of_group_GDB (string expression, vector <GDB_> inGDB);
bool existence_of_groupcodes (vector <GDB_> inGDB);

vector <double> cubic_solution (double A, double B, double C, double D);
vector <double>  quartic_solution (double A, double B, double C, double D, double E);

STRESSFIELD eigenvalue_eigenvector (STRESSTENSOR st);
STRESSTENSOR stresstensor_from_eigenvalue_eigenvector (STRESSFIELD sf);

STRESSFIELD computestressfield_DXDYDZ (STRESSFIELD in);
STRESSFIELD computestressfield_NXNYNZ (STRESSFIELD in);
STRESSFIELD stress_regime (STRESSFIELD in);

STRESSTENSOR invert_stress_tensor (STRESSTENSOR st);

VCTR return_stressvector (STRESSTENSOR st, GDB_ inGDB, bool compression_positive);
VCTR return_normalstress (STRESSTENSOR st, GDB_ inGDB, bool compression_positive);
VCTR return_shearstress (STRESSTENSOR st, GDB_ inGDB, bool compression_positive);
VCTR return_upsilon (STRESSTENSOR st, GDB_ inGDB, string method, bool compression_positive);

double return_ANG (STRESSTENSOR st, GDB_ inGDB, bool compression_positive);
double return_RUP (STRESSTENSOR st, GDB_ inGDB, string method, bool compression_positive);

double return_average_misfit (STRESSTENSOR st, vector <GDB_> inGDB, bool compression_positive);

vector <double> hyperplane_from_GDB (GDB_ inGDB);
vector < vector <double> > shan_matrix_from_GDB (GDB_ inGDB);

#endif