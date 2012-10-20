#ifndef DATA_IO_HPP_
#define DATA_IO_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <map>
#include <cstdlib>

#include "common.h"

using namespace std;

PFN createprojectfoldernames (string projectname);
bool createprojectfolders (PFN output, vector <GDB_> inGDB);
bool copyoriginalfile (PFN output);

void outputrgfheader (ofstream& o, INPSET_ inset);
void outputaverageheader (ofstream& o);

void outputrecord (GDB_ i, ofstream& o, INPSET_ inpset);
void outputveragerecord (GDB_ i, ofstream& o);

void outputresultrgf (PFN output, vector <GDB_> outGDB, bool tilted, INPSET_ inset);

void outputaveragergf (PFN output, vector <GDB_> outGDB);

void outputselected_ps_rgf (PFN output, vector <GDB_> outGDB, vector <GDB_> tiltoutGDB, INPSET_ inset);

void output_to_rgf (PFN output, vector <GDB_> processGDB, INPSET_ inset, bool tilted);
void output_to_ps (PFN output, vector <GDB_> processGDB, vector <GDB_> tiltprocessGDB, INPSET_ inset, PAPER P, CENTER center);

void process_group_by_group (vector <GDB_> outGDB, vector <GDB_> tiltoutGDB, ofstream& o, INPSET_ inset, CENTER center, PAPER P);
void process_one_by_one (GDB_ processGDB, GDB_ tiltprocessGDB, ofstream& o, INPSET_ inset, CENTER center, PAPER P);

void output_elapsed_time (double elapsed_time);

#endif
