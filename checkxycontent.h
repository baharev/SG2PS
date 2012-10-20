#ifndef CHECKXYCONTENT_HPP_
#define CHECKXYCONTENT_HPP_

#include <string>
#include <vector>

#include "common.h"

using namespace std;

bool needxyfile ();
string inputxyfilename ();

bool xyEXISTENCEcheck (string xyname);
bool xyTABcheck (string xyname);
bool xyIDcheck (string xyname);
bool xyCOORDcheck (string xyname);
vector <LOC_X_Y> competeXYdatabase (string xyname);

GDB_ insertxy (GDB_ inGDB, string xyfilename);

bool xyfile_correct (string projectname);
string check_xy_inputs (string inputfilename, bool batch);

#endif
