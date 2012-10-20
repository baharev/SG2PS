#ifndef PS_HPP_
#define PS_HPP_

#include <string>
#include <vector>

#include "structs.h"

using namespace std;

PAPER PS_dimensions (INPSET_ inset);

void PS_header (string DATATYPE, string LOC, string GROUP, ofstream& o, INPSET_ inset, PAPER P);
void PS_border (GDB_ inGDB, ofstream& o, INPSET_ inset, PAPER P);
void PS_stress_scale (ofstream& o, PAPER P);
void ps_clusteringresult (ofstream& o, INPSET_ inset, PAPER P, int clusternumber, double error);

void PS_net (string DATATYPE, string LOC, ofstream& o, INPSET_ inset, PAPER P);
void PS_stressdata (vector <GDB_> inGDB, INPSET_ i, ofstream& o, CENTER center, PAPER P, STRESSFIELD sf, string method);
void PS_stressarrows (ofstream& o, CENTER center, PAPER P, STRESSFIELD sf);

void PS_mohr_circle (vector <GDB_> inGDB, ofstream& o, CENTER center, PAPER P, STRESSFIELD sf, STRESSTENSOR st, bool compression_positive);
void PS_RUP_distribution (vector <GDB_>  inGDB, ofstream& o, CENTER center, PAPER P);
void PS_ANG_distribution (vector <GDB_>  inGDB, ofstream& o, CENTER center, PAPER P);

void PS_stress_state (ofstream& o, PAPER P, CENTER center, STRESSFIELD sf);

void PS_folddata (GDB_ in, ofstream& o, CENTER center, PAPER P);

void PS_lineation (GDB_ i, ofstream& o, INPSET_ inset, CENTER center, STRESSFIELD sf, bool label, string type);
void PS_plane     (GDB_ i, ofstream& o, INPSET_ inset, CENTER center, bool label, string type);
void PS_polepoint (GDB_ i, ofstream& o, INPSET_ inset, CENTER center, bool label, string type);

void PS_sc_arrow (GDB_ i, ofstream& o, INPSET_ inset, CENTER center, VCTR d);
void PS_striaearrow (GDB_ i, ofstream& o, INPSET_ inset, CENTER center, bool label, string offset);

void PS_rosesegment (GDB_ i, ofstream& o, INPSET_ inset, CENTER center, double percentage, double degree, bool c_plane);
void PS_draw_rose_circle_horizontal (ofstream& o, INPSET_ inset, CENTER center, ROSENUMBER percent);
void PS_draw_rose_circle_vertical (ofstream& o, INPSET_ inset, CENTER center, ROSENUMBER percent);

void PS_datanumber_averagebedding (GDB_ i, ofstream& o, INPSET_ inset, PAPER P, CENTER center, size_t j);

void PS_DRAW_PTN (GDB_ i, ofstream& o, INPSET_ inset, CENTER center);
void PS_DRAW_plane (GDB_ i, ofstream& o, INPSET_ inset, CENTER center);
void PS_DRAW_lineation (GDB_ i, ofstream& o, INPSET_ inset, CENTER center);
void PS_DRAW_striae (GDB_ i, ofstream& o, INPSET_ inset, CENTER center);
void PS_DRAW_sc (GDB_ i, ofstream& o, INPSET_ inset, CENTER center);

void PS_idealmovement (vector <GDB_> inGDB, ofstream& o, INPSET_ inset, CENTER center, STRESSTENSOR st);

void PS_DRAW_record (GDB_ i, ofstream& o, INPSET_ inset, CENTER center);

void PS_SYMBOLS_border (ofstream& o, INPSET_ inset, PAPER P);
void PS_SYMBOL_draw_plane (double X, double Y, ofstream& o, PAPER P, string type, string group);

void PS_SYMBOLS_STRIAE (vector <GDB_> inGDB, ofstream& o, INPSET_ inset, PAPER P);
void PS_SYMBOLS_SC (vector <GDB_> inGDB, ofstream& o, INPSET_ inset, PAPER P);
void PS_SYMBOLS_PLANE (vector <GDB_> inGDB, ofstream& o, INPSET_ inset, PAPER P);
void PS_SYMBOLS_LINEATION (vector <GDB_> inGDB, ofstream& o, INPSET_ inset, PAPER P);
void PS_SYMBOLS_HOEPPNER (vector <GDB_> inGDB, ofstream& o, INPSET_ inset, PAPER P);

void PS_SYMBOLS_INVERSION (vector <GDB_> inGDB, ofstream& o, INPSET_ inset, PAPER P);
void PS_SYMBOLS_GROUPS (ofstream& o, INPSET_ inset, PAPER P);
void PS_SYMBOLS_BINGHAM (vector <GDB_> inGDB, ofstream& o, INPSET_ inset, PAPER P);
void PS_SYMBOLS_ROSE (ofstream& o, PAPER P);
void PS_SYMBOLS_LABEL (ofstream& o, INPSET_ inset, PAPER P);

void PS_SYMBOLS (vector <GDB_> inGDB, ofstream& o, INPSET_ inset, PAPER P);

#endif
