#ifndef RGF_HPP_
#define RGF_HPP_

#include <string>
#include <vector>

#include "structs.h"
#include "data_io.h"
#include "common.h"
#include "cluster.h"
#include "angelier.h"

string GC_from_tempGC (string GCtemp);
double DIPDIR_from_DIPDIRtemp (string DIPDIRtemp, string DATAGROUP);
double corrDIPDIR_from_DIPDIR (double DIPDIR, INPSET_ inset, string DATAGROUP);

double DIP_from_DIPtemp (string DIPtemp, string DATAGROUP);
double corrDIP_from_DIP (double DIP, string DATAGROUP);

double LDIR_from_LDIRtemp (string LDIRtemp, string DATAGROUP);
double corrLDIR_from_LDIR (double LDIR, string DATAGROUP);

string produce_LINEATION (string LDIRtemp, string LDIPtemp);
string produce_OFFSET (string SENSEtemp);
vector <GDB_> compete_colorcode (vector <GDB_> inGDB);
vector <GDB_> black_colorcode (vector <GDB_> inGDB);
vector <GDB_> colorcode_grom_groupcode (vector <GDB_> inGDB);

vector <GDB_> competeRGFcontect (string projectname, string inputxyfilename, INPSET_ inSET);

double right_hand_rule_to_german (double corrDIPDIR);
double german_to_right_hand_rule (double corrDIPDIR);

string cGc_datagroup (string DATATYPE);

vector <GDB_> cGc_NDS (vector <GDB_> inGDB);

GDB_ cGc_NCDCSC_LINEATION_SC (GDB_ inGDB);
GDB_ cGc_NCDCSC_PITCH (GDB_ inGDB);
vector <GDB_> manipulate_N (vector <GDB_> inGDB);

vector <GDB_> cGc_NDS_DCNCSC (vector <GDB_> inGDB);
vector <GDB_> cGc_PITCHANGLE (vector <GDB_> inGDB);

vector <GDB_> cGc_MISFIT (vector <GDB_> inGDB);

vector <GDB_> cGc_UP (vector <GDB_> inGDB);
vector <GDB_> cGc_tilted_UP (vector <GDB_> inGDB);

vector <GDB_> cGc_OFFSET (vector <GDB_> inGDB);
vector <GDB_> cGc_LAMBDA_STRESSVECTOR_ESTIMATORS (vector <GDB_> inGDB);

CORRECTSTRIAE cGc_correct_striae_DIPcor (GDB_ inGDB);
CORRECTSTRIAE cGc_correct_striae_DIPDIRcor (GDB_ inGDB);
vector <GDB_> cGc_striae_correction (vector <GDB_> inGDB);

bool byLocType(const GDB_& x, const GDB_& y);
bool byLocTypeGc(const GDB_& x, const GDB_& y);
bool byiID(const GDB_& x, const GDB_& y);
bool byeigenvalue(const sort_jacobi& x, const sort_jacobi& y);
vector <GDB_> sort_by_iID (vector <GDB_> inGDB);

bool stopcriteria (string prevDATATYPE, string DATATYPE, string prevLOC, string LOC, string prevGC, string GC, INPSET_ inset);
bool stopcriteria (string prevDATATYPE, string DATATYPE, string prevLOC, string LOC);

void fold_from_planes (vector <GDB_> inGDB, ofstream& o, INPSET_ inset, CENTER center, PAPER P);

vector <GDB_> cGc_average (vector <GDB_> inGDB);
vector <GDB_> cGc_s0_average (vector <GDB_> inGDB);

GDB_ plane_tilt (GDB_ inGDB, bool paleonorth);
GDB_ lineation_tilt (GDB_ inGDB, bool paleonorth);
GDB_ SC_tilt (GDB_ inGDB, bool paleonorth);
GDB_ striae_tilt (GDB_ inGDB, bool paleonorth);
GDB_ S0_TILT (GDB_ inGDB, bool paleonorth);
vector <GDB_> cGc_RETILT (vector <GDB_> inGDB, INPSET_ inSET);

vector <GDB_> ptn (vector <GDB_> inGDB, INPSET_ inset);

vector <GDB_> clustering_GBD (INPSET_ inset, vector <GDB_> inGDB);

void process_rgf (string inputfilename, string XY_filename, INPSET_ inset);

#endif
