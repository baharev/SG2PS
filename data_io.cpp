#include <iomanip>
#include <iostream>
#include <algorithm>

#include "data_io.h"
#include "ps.h"
#include "rose.h"
#include "rgf.h"
#include "angelier.h"
#include "cluster.h"

using namespace std;

PFN createprojectfoldernames (string projectname) {

	time_t current_time;
	struct tm * TM;
	const string bs = char_to_string ('\\');

	time ( &current_time );
	TM = localtime ( &current_time );

	PFN output;

	output.datetime = int_to_string (TM->tm_year + 1900);

	if (TM->tm_mon < 9) output.datetime = output.datetime + "0";
	output.datetime = output.datetime + int_to_string (TM->tm_mon + 1);

	if (TM->tm_mday < 10) output.datetime = output.datetime + "0";
	output.datetime = output.datetime + int_to_string (TM->tm_mday) + "-";

	if (TM->tm_hour < 10) output.datetime = output.datetime + "0";
	output.datetime = output.datetime + int_to_string (TM->tm_hour);

	if (TM->tm_min < 10) output.datetime = output.datetime + "0";
	output.datetime = output.datetime + int_to_string (TM->tm_min);

	if (TM->tm_sec < 10) output.datetime = output.datetime + "0";
	output.datetime = output.datetime + int_to_string (TM->tm_sec);

	output.projectfolder 	= output.datetime + "_-_" + capslock(projectname);
	output.projectname 		= capslock(projectname);
	output.original			= output.projectfolder +  bs + "1_original";
	output.completed		= output.projectfolder +  bs + "2_completed";
	output.average			= output.projectfolder +  bs + "3_average";
	output.rgfsep			= output.projectfolder +  bs + "4_rgf_separated";
	output.pssep			= output.projectfolder +  bs + "5_ps_separated";

	return output;
}

bool createprojectfolders (PFN output, vector <GDB_> inGDB) {

	const string bs = char_to_string ('\\');
	int returncode = 0;

	returncode = system (("md " + output.projectfolder).c_str());
	returncode = system (("md " + output.original).c_str());
	returncode = system (("md " + output.completed).c_str());
	returncode = system (("md " + output.average).c_str());
	returncode = system (("md " + output.rgfsep).c_str());
	returncode = system (("md " + output.pssep).c_str());

	if (existence ("BOUDAIN", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "boudain").c_str());
		returncode = system (("md " + output.pssep + bs + "boudain").c_str());
	}

	if (existence ("CONTACT", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "contact").c_str());
		returncode = system (("md " + output.pssep + bs + "contact").c_str());
	}

	if (existence ("FOLDAXIS", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "foldaxis").c_str());
		returncode = system (("md " + output.pssep + bs + "foldaxis").c_str());
	}

	if (existence ("FOLDPLANE", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "foldplane").c_str());
		returncode = system (("md " + output.pssep + bs + "foldplane").c_str());
	}

	if (existence ("KINK", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "kink").c_str());
		returncode = system (("md " + output.pssep + bs + "kink").c_str());
	}

	if (existence ("LINEATION", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "lineation").c_str());
		returncode = system (("md " + output.pssep + bs + "lineation").c_str());
	}

	if (existence ("LITHOCLASE", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "lithoclase").c_str());
		returncode = system (("md " + output.pssep + bs + "lithoclase").c_str());
	}
	if (existence ("SC", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "sc").c_str());
		returncode = system (("md " + output.pssep + bs + "sc").c_str());
	}
	if (existence ("BEDDING", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "bedding").c_str());
		returncode = system (("md " + output.pssep + bs + "bedding").c_str());
	}
	if (existence ("S1", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "s1").c_str());
		returncode = system (("md " + output.pssep + bs + "s1").c_str());
	}

	if (existence ("S2", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "s2").c_str());
		returncode = system (("md " + output.pssep + bs + "s2").c_str());
	}

	if (existence ("S3", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "s3").c_str());
		returncode = system (("md " + output.pssep + bs + "s3").c_str());
	}

	if (existence ("S4", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "s4").c_str());
		returncode = system (("md " + output.pssep + bs + "s4").c_str());
	}

	if (existence ("S5", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "s5").c_str());
		returncode = system (("md " + output.pssep + bs + "s5").c_str());
	}

	if (existence ("FRACTURE", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "fracture").c_str());
		returncode = system (("md " + output.pssep + bs + "fracture").c_str());
	}

	if (existence ("STRIAE", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "striae").c_str());
		returncode = system (("md " + output.pssep + bs + "striae").c_str());
	}

	if (existence ("CROSSBEDDING", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "crossbedding").c_str());
		returncode = system (("md " + output.pssep + bs + "crossbedding").c_str());
	}

	if (existence ("VEIN", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "vein").c_str());
		returncode = system (("md " + output.pssep + bs + "vein").c_str());
	}

	if (existence ("FOLDSURFACE", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "foldsurface").c_str());
		returncode = system (("md " + output.pssep + bs + "foldsurface").c_str());
	}

	if (existence ("USERPLANE4", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "userplane4").c_str());
		returncode = system (("md " + output.pssep + bs + "userplane4").c_str());
	}

	if (existence ("USERPLANE5", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "userplane5").c_str());
		returncode = system (("md " + output.pssep + bs + "userplane5").c_str());
	}

	if (existence ("USERLINEATION1", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "userlineation1").c_str());
		returncode = system (("md " + output.pssep + bs + "userlineation1").c_str());
	}

	if (existence ("USERLINEATION2", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "userlineation2").c_str());
		returncode = system (("md " + output.pssep + bs + "userlineation2").c_str());
	}

	if (existence ("USERLINEATION3", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "userlineation3").c_str());
		returncode = system (("md " + output.pssep + bs + "userlineation3").c_str());
	}

	if (existence ("USERLINEATION4", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "userlineation4").c_str());
		returncode = system (("md " + output.pssep + bs + "userlineation4").c_str());
	}

	if (existence ("USERLINEATION5", inGDB)) {

		returncode = system (("md " + output.rgfsep + bs + "userlineation5").c_str());
		returncode = system (("md " + output.pssep + bs + "userlineation5").c_str());
	}

	if (returncode != 0) {

		cout << "Cannot create project folder." << endl;

		return false;
	}

	return true;
}

bool copyoriginalfile (PFN output) {

	ifstream inrgffile;
	ofstream outrgffile;

	const string bs = char_to_string ('\\');
	string buffer;

	string infilename = (output.projectname + ".rgf").c_str();
	string outfilename = (output.original + bs + output.projectname + ".rgf").c_str();

	inrgffile.open (infilename.c_str());
	outrgffile.open (outfilename.c_str());

	if (!(inrgffile.is_open())) {

		cout << "  - ERROR: cannot open input file.";
		return false;
	}

	if (!(outrgffile.is_open())) {

		cout << "  - ERROR: cannot create output file in project destination folder " + output.original + "." << endl;
		return false;
	}

	do {

		getline (inrgffile, buffer);
		if (inrgffile.eof ()) {

			outrgffile << buffer;
			break;
		}

		else {

			outrgffile << buffer << endl;
		}
	}

	while (!(inrgffile.eof()));

	inrgffile.close();
	outrgffile.close();

	return true;
}

void outputrgfheader (ofstream& o, INPSET_ inset) {

	o
	<< "ID" 		<< '\t'
	<< "GC" 		<< '\t'
	<< "COLOR" 		<< '\t'
	<< "LOC" 		<< '\t'
	<< "LOCX" 		<< '\t'
	<< "LOCY" 		<< '\t'
	<< "FORMATION" 	<< '\t'
	<< "DATATYPE"	<< '\t' << flush;

	if (inset.datarule == "G") 	o << "corrSTRIKEDIR" 	<< '\t' << flush;
	else 						o << "corrDIPDIR" 		<< '\t' << flush;

	o
	<< "corrDIP" 	<< '\t'
	<< "corrLDIR" 	<< '\t'
	<< "corrLDIP" 	<< '\t'
	<< "SENSE" 		<< '\t'
	<< "PALEONORTH"	<< '\t'
	<< "COMMENT"

	<< endl;
}

void outputaverageheader (ofstream& o) {

	o
	<< "ID" << '\t'
	<< "GC" << '\t'
	<< "COLOR" << '\t'
	<< "LOC" << '\t'
	<< "LOCX" << '\t'
	<< "LOCY" << '\t'
	<< "FORMATION" << '\t'
	<< "DATATYPE" << '\t'
	<< "averageDIPDIR" << '\t'
	<< "averageDIP" << '\t'
	<< "NODATA" << '\t'
	<< "NODATA" << '\t'
	<< "NODATA" << '\t'
	<< "NODATA" << '\t'
	<< "COMMENT" << endl;
}

void outputrecord (GDB_ i, ofstream& o, INPSET_ inpset) {

	o
	<< i.ID << '\t'
	<< i.GC << '\t'
	<< i.COLOR << '\t'
	<< i.LOC << '\t'
	<< fixed << setprecision (6) << i.LOCX << '\t'
	<< fixed << setprecision (6) << i.LOCY << '\t'
	<< i.FORMATION << '\t'
	<< i.DATATYPE << '\t' << flush;

	o << fixed << setprecision (1) << flush;

	if ((i.corr.DIPDIR > 361.0) || (i.DATATYPE == "LITHOLOGY")) o << "" << '\t' << flush;
	else {

		if (inpset.datarule == "R" ) 	o << german_to_right_hand_rule (i.corr.DIPDIR) << '\t' << flush;
		else							o << i.corr.DIPDIR << '\t' << flush;
	}

	if ((i.corr.DIP > 360.0) || (i.DATATYPE == "LITHOLOGY")) o << "" << '\t' << flush;
	else o << fixed << setprecision (1) << i.corr.DIP << '\t' << flush;

	if (i.corrL.DIPDIR > 360.0) o << "" << '\t' << flush;
	else o << i.corrL.DIPDIR << '\t' << flush;

	if (i.corrL.DIP > 360.0) o << "" << '\t' << flush;
	else o << i.corrL.DIP << '\t' << flush;

	if ((i.DATATYPE == "STRIAE") || (i.DATATYPE == "BEDDING")) {

		if (!((i.DATATYPE == "BEDDING")) && (i.OFFSET == "NONE")) o << i.OFFSET << '\t' << flush;
	}

	else o << "" << '\t' << flush;

	o
	<< i.PALEON	<< '\t'
	<< i.COMMENT<<

	flush;
}

void outputveragerecord (GDB_ i, ofstream& o) {

	o
	<< i.ID << '\t'
	<< i.GC << '\t'
	<< i.COLOR << '\t'
	<< i.LOC << '\t'
	<< setprecision (6) << i.LOCX << '\t'
	<< setprecision (6) << i.LOCY << '\t'
	<< i.FORMATION << '\t'
	<< i.DATATYPE << '\t' << flush;

	if (i.avd.DIPDIR > 361.0) o << "" << '\t' << flush;
	else o << i.avd.DIPDIR << '\t' << flush;

	if (i.avd.DIP > 361.0) o << "" << '\t' << flush;
	else o << i.avd.DIP << '\t' << flush;

	o << "" << '\t' << flush;
	o << "" << '\t' << flush;
	o << "" << '\t' << flush;
	o << "" << '\t' << flush;
	o << i.COMMENT<< '\t' << flush;
}

void outputresultrgf (PFN output, vector <GDB_> outGDB, bool tilted, INPSET_ inset) {

	ofstream outputfile;
	string outputfilename;
	string bs = char_to_string ('\\');

	size_t i = 0;

	outputfilename = output.completed + bs + output.projectname + "_completed";
	if (tilted) outputfilename = outputfilename + "_tilted";
	outputfilename = outputfilename + ".rgf";

	outputfile.open (outputfilename.c_str());

	outputrgfheader (outputfile, inset);

	while (i < outGDB.size()) {

		outputrecord (outGDB.at(i), outputfile, inset);

		if (i < outGDB.size() - 1) outputfile << endl;

		i++;
	}

	outputfile.close();

	if (tilted) 	cout << "  - Tilted RGF file exported." <<  endl;
	else 			cout << "  - Completed RGF file exported." <<  endl;
}

void outputaveragergf (PFN output, vector <GDB_> outGDB) {

	ofstream outputfile;
	string bs = char_to_string ('\\');
	string outputfilename = output.average + bs + output.projectname + "_average.rgf";

	size_t i = 0;
	size_t independentrecordcounter = 0;

	outputfile.open (outputfilename.c_str());
	outputaverageheader (outputfile);

	if ((outGDB.size() == 1) && (!((outGDB[0].DATATYPE == "STRIAE") || (outGDB[0].DATATYPE == "SC")))) {

		outputveragerecord (outGDB.at(0), outputfile);
		outputfile.close();
		cout << "  - Average RGF output completed." <<  endl;
	}

	do {

		do {

			i++;
			if (i == outGDB.size()) break;
		}

		while ((outGDB.at(i-1).DATATYPE == outGDB.at(i).DATATYPE) && (outGDB.at(i-1).LOC == outGDB.at(i).LOC));

		independentrecordcounter++;

		if (!((outGDB[i-1].DATATYPE == "STRIAE") || (outGDB[i-1].DATATYPE == "SC"))) {

			outputveragerecord (outGDB[i - 1], outputfile);
		}

		if (i < outGDB.size()) 	outputfile << endl;
	}

	while (i < outGDB.size());

	outputfile.close();

	cout << "  - Average RGF output completed." <<  endl;
}

void outputselected_ps_rgf (PFN output, vector <GDB_> outGDB, vector <GDB_> tiltoutGDB, INPSET_ inset) {

	vector <GDB_> processGDB, tiltprocessGDB;
	size_t i = 0;
	size_t independentrecordcounter = 0;
	CENTER center;
	PAPER P = PS_dimensions (inset);

	center.X = P.O1X;
	center.Y = P.O1Y;

	center.radius = P.R;

	do {

		processGDB.clear();
		tiltprocessGDB.clear();

		do {

			processGDB.push_back(outGDB[i]);
			tiltprocessGDB.push_back(tiltoutGDB[i]);

			i++;

			if (i == outGDB.size()) break;

		} while (stopcriteria (outGDB.at(i-1).DATATYPE, outGDB.at(i).DATATYPE, outGDB.at(i-1).LOC, outGDB.at(i).LOC, outGDB.at(i-1).GC, outGDB.at(i).GC, inset));

		independentrecordcounter++;

		if (processGDB.at(0).DATATYPE != "LITHOLOGY") {

			if (inset.group == "N") {

				if (existence_of_groupcodes (processGDB)) {

					processGDB = 		colorcode_grom_groupcode (processGDB);
					tiltprocessGDB = 	colorcode_grom_groupcode (tiltprocessGDB);
				}
			}

			else {

				processGDB = 		black_colorcode (processGDB);
				tiltprocessGDB = 	black_colorcode (tiltprocessGDB);
			}

			output_to_rgf (output, processGDB, inset, false);

			output_to_rgf (output, tiltprocessGDB, inset, true);

			output_to_ps (output, processGDB, tiltprocessGDB, inset, P, center);
		}

	} while (i < outGDB.size());

	cout << "DATA EXPORT FROM '" << capslock(output.projectname) << ".RGF' DATABASE FILE" << endl;
	cout << "  - Postscript output completed for " << independentrecordcounter << " file containing " << i << " records." <<  endl;
	cout << "  - RGF output completed for " << independentrecordcounter << " file containing " << i << " records." <<  endl;
	cout << "  - Tilted RGF output completed for " << independentrecordcounter << " file containing " << i << " records." <<  endl;
}

void output_to_rgf (PFN output, vector <GDB_> processGDB, INPSET_ inset, bool tilted) {

	ofstream output_rgf_file;
	string output_rgf_filename;
	string bs = char_to_string ('\\');
	size_t j = 0;

	if (inset.group == "Y") output_rgf_filename = output.rgfsep + bs + processGDB.at(0).DATATYPE + bs + processGDB.at(0).LOC + "_" + processGDB.at(0).DATATYPE + "_" + processGDB.at(0).GC;
	else 					output_rgf_filename = output.rgfsep + bs + processGDB.at(0).DATATYPE + bs + processGDB.at(0).LOC + "_" + processGDB.at(0).DATATYPE;

	if (tilted) output_rgf_filename = output_rgf_filename + "_tilted";
	else output_rgf_filename = output_rgf_filename + "_original";

	output_rgf_filename = output_rgf_filename + ".rgf";

	output_rgf_file.open (output_rgf_filename.c_str());

	sort(processGDB.begin(), processGDB.end(), byiID);

	outputrgfheader (output_rgf_file, inset);

	do {

		outputrecord (processGDB.at(j), output_rgf_file, inset);
		if (j < processGDB.size()-1) output_rgf_file << endl;
		j++;

	} while (j < processGDB.size());

	output_rgf_file.close();
}

void output_to_ps (PFN output, vector <GDB_> processGDB, vector <GDB_> tiltprocessGDB, INPSET_ inset, PAPER P, CENTER center) {

	ofstream output_ps_file, output_rgf_file, output_tiltedrgf_file;
	string output_rgf_filename,  output_tiltedrgf_filename, output_ps_filename;
	string bs = char_to_string ('\\');
	size_t j = 0;

	if (inset.group == "Y") output_ps_filename = output.pssep +  bs + processGDB.at(0).DATATYPE + bs + processGDB.at(0).LOC + "_" + processGDB.at(0).DATATYPE + "_" + processGDB.at(0).GC + ".ps";
	else 					output_ps_filename = output.pssep +  bs + processGDB.at(0).DATATYPE + bs + processGDB.at(0).LOC + "_" + processGDB.at(0).DATATYPE + ".ps";

	output_ps_file.open (output_ps_filename.c_str());
	PS_header (processGDB.at(0).DATATYPE, processGDB.at(0).LOC, processGDB.at(0).GC, output_ps_file, inset, P);
	PS_SYMBOLS(processGDB, output_ps_file, inset, P);

	if (processGDB.at(0).DATATYPE == "STRIAE") PS_stress_scale (output_ps_file, P);

	PS_border (processGDB.at(0), output_ps_file, inset, P);

	do {

		process_one_by_one (processGDB.at(j), tiltprocessGDB.at(j), output_ps_file, inset, center, P);
		if (j < processGDB.size()-1) output_ps_file << endl;
		j++;

	} while (j < processGDB.size());

	process_group_by_group (processGDB, tiltprocessGDB, output_ps_file, inset, center, P);

	PS_datanumber_averagebedding (processGDB.at(0), output_ps_file, inset, P, center, processGDB.size());

	PS_net (processGDB.at(0).DATATYPE, processGDB.at(0).LOC, output_ps_file, inset, P);
	output_ps_file.close();
}

void process_group_by_group (vector <GDB_> outGDB, vector <GDB_> tiltoutGDB, ofstream& o, INPSET_ inset, CENTER center, PAPER P) {

	CENTER mohr_center;

	PS_draw_rose (outGDB, tiltoutGDB, o, inset, center, P);

	if ((inset.fracture == "B") && (outGDB.at(0).DATATYPE == "FRACTURE") && (outGDB.size() < 2)) return;
	if ((inset.fracture == "B") && (outGDB.at(0).DATATYPE == "FRACTURE") && (tiltoutGDB.size() < 2)) return;

	if ((inset.fracture == "B") && (outGDB.at(0).DATATYPE == "FRACTURE")) {

		cout << "  - For '" << outGDB.at(0).LOC << "' location" << flush;
		if (inset.group == "Y")	cout << ", '"<< outGDB.at(0).GC << "', " << flush;
		else 					cout << "," << flush;
		cout << "fracture statistics after Bingham (1964): " << endl;

		if (outGDB.size() > 1) {

			cout << "    - Original : " << flush;
			center.X = P.O1X;
			center.Y = P.O1Y;
			outGDB = inversion ("BINGHAM", outGDB, o, inset, center, mohr_center, P);

			cout << "    - Corrected: " << flush;
			center.X = P.O2X;
			center.Y = P.O2Y;
			tiltoutGDB = inversion ("BINGHAM", tiltoutGDB, o, inset, center, mohr_center, P);
		}

		else cout << "less data than required to the statistics." << endl;
	}

	else if ((inset.inversion == "F") && ((outGDB.at(0).DATATYPE == "STRIAE"))) {

		cout << "  - For '" << outGDB.at(0).LOC << "' location" << flush;
		if (inset.group == "Y")	cout << ", '"<< outGDB.at(0).GC << "', " << flush;
		else 					cout << "," << flush;
		cout << "regression after Fry (1999): " << endl;

		if (outGDB.size() > 5) {

			cout << "    - Original : " << flush;
			center.X = P.O1X;
			center.Y = P.O1Y;
			mohr_center.X = P.O7X;
			mohr_center.Y = P.O7Y;
			outGDB = inversion ("FRY", outGDB, o, inset, center, mohr_center, P);

			cout << "    - Corrected: " << flush;
			center.X = P.O2X;
			center.Y = P.O2Y;
			mohr_center.X = P.O8X;
			mohr_center.Y = P.O8Y;
			tiltoutGDB = inversion ("FRY", tiltoutGDB, o, inset, center, mohr_center, P);
		}

		else cout << "less data than required to the statistics." << endl;
	}

	else if ((inset.inversion == "M") && ((outGDB.at(0).DATATYPE == "STRIAE"))) {

		cout << "  - For '" << outGDB.at(0).LOC << "' location" << flush;
		if (inset.group == "Y")	cout << ", '"<< outGDB.at(0).GC << "', " << flush;
		else 					cout << "," << flush;
		cout << "regression after Michael (1984): " << endl;

		if (outGDB.size() > 4) {

			cout << "    - Original : " << flush;
			center.X = P.O1X;
			center.Y = P.O1Y;
			mohr_center.X = P.O7X;
			mohr_center.Y = P.O7Y;
			outGDB = inversion ("MICHAEL", outGDB, o, inset, center, mohr_center, P);

			cout << "    - Corrected: " << flush;
			center.X = P.O2X;
			center.Y = P.O2Y;
			mohr_center.X = P.O8X;
			mohr_center.Y = P.O8Y;
			tiltoutGDB = inversion ("MICHAEL", tiltoutGDB, o, inset, center, mohr_center, P);
		}

		else cout << "less data than required to the statistics." << endl;
	}

	else if ((inset.inversion == "S") && ((outGDB.at(0).DATATYPE == "STRIAE"))) {

		cout << "  - For '" << outGDB.at(0).LOC << "' location" << flush;
		if (inset.group == "Y")	cout << ", '"<< outGDB.at(0).GC << "', " << flush;
		else 					cout << "," << flush;
		cout << "regression after Shan et al. (2003): " << endl;

		if (outGDB.size() > 4) {

			cout << "    - Original : " << flush;
			center.X = P.O1X;
			center.Y = P.O1Y;
			mohr_center.X = P.O7X;
			mohr_center.Y = P.O7Y;
			outGDB = inversion ("SHAN", outGDB, o, inset, center, mohr_center, P);

			cout << "    - Corrected: " << flush;
			center.X = P.O2X;
			center.Y = P.O2Y;
			mohr_center.X = P.O8X;
			mohr_center.Y = P.O8Y;
			tiltoutGDB = inversion ("SHAN", tiltoutGDB, o, inset, center, mohr_center, P);
		}

		else cout << "less data than required to the statistics." << endl;
	}

	else if ((inset.inversion == "A") && ((outGDB.at(0).DATATYPE == "STRIAE"))) {

		cout << "  - For '" << outGDB.at(0).LOC << "' location" << flush;
		if (inset.group == "Y")	cout << ", '"<< outGDB.at(0).GC << "', " << flush;
		else 					cout << "," << flush;
		cout << "inversion after Angelier (1990): " << endl;

		if (outGDB.size() > 3) {

			cout << "    - Original : " << flush;
			center.X = P.O1X;
			center.Y = P.O1Y;
			mohr_center.X = P.O7X;
			mohr_center.Y = P.O7Y;
			outGDB = inversion ("ANGELIER", outGDB, o, inset, center, mohr_center, P);

			cout << "    - Corrected: " << flush;
			center.X = P.O2X;
			center.Y = P.O2Y;
			mohr_center.X = P.O8X;
			mohr_center.Y = P.O8Y;
			tiltoutGDB = inversion ("ANGELIER", tiltoutGDB, o, inset, center, mohr_center, P);
		}

		else cout << "less data than required to the statistics." << endl;
	}

	else if ((inset.inversion == "O") && ((outGDB.at(0).DATATYPE == "STRIAE"))) {

		cout << "  - For '" << outGDB.at(0).LOC << "' location" << flush;
		if (inset.group == "Y")	cout << ", '"<< outGDB.at(0).GC << "', " << flush;
		else 					cout << "," << flush;
		cout << "inversion after Mostafa (2005): " << endl;

		if (outGDB.size() > 3) {

			cout << "    - Original : " << flush;
			center.X = P.O1X;
			center.Y = P.O1Y;
			mohr_center.X = P.O7X;
			mohr_center.Y = P.O7Y;
			outGDB = inversion ("MOSTAFA", outGDB, o, inset, center, mohr_center, P);

			cout << "    - Corrected: " << flush;
			center.X = P.O2X;
			center.Y = P.O2Y;
			mohr_center.X = P.O8X;
			mohr_center.Y = P.O8Y;
			tiltoutGDB = inversion ("MOSTAFA", tiltoutGDB, o, inset, center, mohr_center, P);
		}

		else cout << "less data than required to the statistics." << endl;
	}

	else if ((inset.inversion == "D") && ((outGDB.at(0).DATATYPE == "STRIAE"))) {

		cout << "  - For '" << outGDB.at(0).LOC << "' location" << flush;
		if (inset.group == "Y")	cout << ", '"<< outGDB.at(0).GC << "', " << flush;
		else 					cout << "," << flush;
		cout << "regression after Sprang (1972): " << endl;

		cout << "    - Original : " << flush;
		center.X = P.O1X;
		center.Y = P.O1Y;
		mohr_center.X = P.O7X;
		mohr_center.Y = P.O7Y;
		outGDB = inversion ("NDA", outGDB, o, inset, center, mohr_center, P);

		cout << "    - Corrected: " << flush;
		center.X = P.O2X;
		center.Y = P.O2Y;
		mohr_center.X = P.O8X;
		mohr_center.Y = P.O8Y;
		tiltoutGDB = inversion ("NDA", tiltoutGDB, o, inset, center, mohr_center, P);
	}

	else if ((inset.inversion == "P") && ((outGDB.at(0).DATATYPE == "STRIAE"))) {

		cout << "  - For '" << outGDB.at(0).LOC << "' location" << flush;
		if (inset.group == "Y")	cout << ", '"<< outGDB.at(0).GC << "', " << flush;
		else 					cout << "," << flush;
		cout << " regression after Turner (1953): " << endl;

		cout << "    - Original : " << flush;
		center.X = P.O1X;
		center.Y = P.O1Y;
		mohr_center.X = P.O7X;
		mohr_center.Y = P.O7Y;
		outGDB = inversion ("PTN", outGDB, o, inset, center, mohr_center, P);

		cout << "    - Corrected: " << flush;
		center.X = P.O2X;
		center.Y = P.O2Y;
		mohr_center.X = P.O8X;
		mohr_center.Y = P.O8Y;
		tiltoutGDB = inversion ("PTN", tiltoutGDB, o, inset, center, mohr_center, P);
	}

	else {}

	if  (outGDB[0].DATATYPE == "FOLDSURFACE") {

		cout << "  - For '" << outGDB.at(0).LOC << "' location" << flush;
		if (inset.group == "Y")	cout << ", '"<< outGDB.at(0).GC << "', " << flush;
		else 					cout << "," << flush;
		cout << " fold axis calculation: " << endl;

		if (outGDB.size() > 1) {

			cout << "    - Original : " << flush;
			center.X = P.O1X;
			center.Y = P.O1Y;
			fold_from_planes (outGDB, o, inset, center, P);

			cout << "    - Corrected: " << flush;
			center.X = P.O2X;
			center.Y = P.O2Y;
			fold_from_planes (tiltoutGDB, o, inset, center, P);
		}

		else cout << "less data than required to the statistics." << endl;
	}

}

void process_one_by_one (GDB_ processGDB, GDB_ tiltprocessGDB, ofstream& o, INPSET_ inset, CENTER center, PAPER P) {

	center.X = P.O1X;
	center.Y = P.O1Y;

	PS_DRAW_record (processGDB, o, inset, center);

	center.X = P.O2X;
	center.Y = P.O2Y;
	PS_DRAW_record (tiltprocessGDB, o, inset, center);
}


void output_elapsed_time (double elapsed_time) {

	if (elapsed_time < 1 * 1000.0) cout << "  - Elapsed time: " << fixed << setprecision (2) << elapsed_time << " milliseconds." << endl;
	else {
		elapsed_time = elapsed_time / 1000.0;
		if (elapsed_time < 1 * 60.0) cout << "  - Elapsed time: " << fixed << setprecision (2) << elapsed_time << " seconds." << endl;
		else {
			elapsed_time = elapsed_time / 60.0;
			if (elapsed_time < 1 * 60.0) cout << "  - Elapsed time: " << fixed << setprecision (2) << elapsed_time << " minutes." << endl;
			else {
				elapsed_time = elapsed_time / 60.0;
				if (elapsed_time < 1 * 60.0) cout << "  - Elapsed time: " << fixed << setprecision (2) << elapsed_time << " hours." << endl;
				else cout << "  - Elapsed time: " << fixed << setprecision (1) << elapsed_time / 60.0 << " days." << endl;
			}
		}
	}
}