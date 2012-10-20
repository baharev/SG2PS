#ifndef ANGELIER_HPP_
#define ANGELIER_HPP_

#include <string>
#include <vector>
#include <cmath>

#include "structs.h"
#include "common.h"


using namespace std;

vector <vector < double> > michael_parameters (vector <GDB_> inGDB); //OK
vector <vector < double> > stressvector_parameters (vector <GDB_> inGDB);

ANG_PRM angelier_parameters (vector <GDB_> inGDB);
STRESSTENSOR compute_angelier_stresstensor (ANG_PRM parameters, vector <GDB_> inGDB);

STRESSTENSOR ptn_P (vector <GDB_> inGDB);
STRESSTENSOR ptn_T (vector <GDB_> inGDB);
STRESSTENSOR ptn_N (vector <GDB_> inGDB);

vector <vector <double> > FRY (vector <GDB_> inGDB, INPSET_ inset);

vector <vector <double> > SHAN (vector <GDB_> inGDB, INPSET_ inset);

STRESSTENSOR ANGELIER (vector <GDB_> inGDB, INPSET_ inset);

STRESSTENSOR MICHAEL (vector <GDB_> inGDB, INPSET_ inset);
STRESSFIELD MICHAEL_PROCESS (vector <GDB_> inGDB, INPSET_ inset);

STRESSTENSOR NDA (vector <GDB_> inGDB, INPSET_ inset);
STRESSFIELD NDA_PROCESS (vector <GDB_> inGDB, INPSET_ inset);

STRESSTENSOR BINGHAM (vector <GDB_> inGDB);
STRESSFIELD BINGHAM_PROCESS (vector <GDB_> inGDB);

vector <GDB_> return_stressvector_estimators (STRESSTENSOR st, vector <GDB_> inGDB, string method, bool compression_positive);

vector <GDB_> generate_virtual_striae (vector <GDB_> inGDB);

vector <GDB_> inversion (string method, vector <GDB_> inGDB, ofstream& o, INPSET_ inset, CENTER center, CENTER mohr_center, PAPER P);

#endif
