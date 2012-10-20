#ifndef CLUSTER_HPP_
#define CLUSTER_HPP_

#include <string>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>


#include "structs.h"


vector <CENTR_VECT>  init_centriod (int cluster_number, vector <GDB_> inGDB);
vector <vector <double> > init_distance_matrix (int cluster_number, vector <GDB_> inGDB);
vector <int> init_whichgroup (int cluster_number, vector <GDB_> inGDB);
vector <int> init_whichgroup_II (int cluster_number, vector <GDB_> inGDB);

double compute_distance (CENTR_VECT centroid, GDB_ inGDB);

vector <CENTR_VECT> compute_centroid_from_which_group (int cluster_number, vector <int> which_group, vector <GDB_> inGDB);
vector <vector <double> > compute_distance_matrix_from_centroid (vector <vector <double> > distance_matrix, vector <GDB_> inGDB, vector <CENTR_VECT> centroid);
vector <int> compute_whichgroup_from_distances (vector <vector <double> > distance_matrix, vector <GDB_> inGDB);

vector <GDB_> attach_group_codes (vector <int> which_group, vector <GDB_> inGDB);

double cumulative_distance (vector <vector <double> > distance_matrix, vector <int> which_group);

vector <vector <double> > clustering_cycle (int cluster_number, vector <GDB_> inGDB, INPSET_ i);
int autoclustering_cycle (int cluster_number, vector <GDB_> inGDB);

vector <GDB_> k_means_clustering (int cluster_number, vector <GDB_> inGDB,  INPSET_ i);
vector <GDB_> k_means (INPSET_ i, vector <GDB_> inGDB);

#endif
