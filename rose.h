#ifndef ROSE_HPP_
#define ROSE_HPP_

#include <vector>

#include "structs.h"
#include "common.h"

ROSENUMBER compute_data_number_DIPDIR (vector <GDB_> inGDB, double strike_begin, double strike_end);
ROSENUMBER compute_data_number_DIP (vector <GDB_> inGDB, double strike_begin, double strike_end);

void PS_draw_rose_PLANE (GDB_ inGDB, ofstream& o, INPSET_ inset, CENTER center, ROSENUMBER percent, double begin_angle, bool vertical);
void PS_draw_rose_LINEATION (GDB_ inGDB, ofstream& o, INPSET_ inset, CENTER center, ROSENUMBER percent, double begin_angle, bool vertical);
void PS_draw_rose_STRIAE (GDB_ inGDB, ofstream& o, INPSET_ inset, CENTER center, ROSENUMBER percent, double begin_angle, bool vertical);

void PS_draw_rose_DIP_DIR (vector <GDB_> inGDB, ofstream& o, INPSET_ inset, CENTER center);
void PS_draw_rose_DIP (vector <GDB_> inGDB, ofstream& o, INPSET_ inset, CENTER center);

void PS_draw_rose (vector <GDB_> roseGDB, vector <GDB_> tiltroseGDB, ofstream& o, INPSET_ inset, CENTER center, PAPER P);

#endif
