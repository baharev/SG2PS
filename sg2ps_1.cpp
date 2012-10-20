#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <algorithm>


#include "common.h"
#include "rgf.h"
#include "data_io.h"
#include "checksettingfilecontent.h"
#include "checkrgffilecontent.h"
#include "checkxycontent.h"
#include "cluster.h"



using namespace std;

int main (int argument_number, char *argv[]) {

	string inputrgfname, xy_filename, inputrgfname_only, temp;
	vector <GDB_> geodatabase, tiltgeodatabase;
	INPSET_ inset;
	size_t finishtime, starttime;
	size_t j = 1;
	bool batch = false;
	bool using_xy_files = false;
	vector <string> inputfilename_vector;
	vector <bool> xy_file_ok;
	double elapsed_time;


	header ();

	inputfilename_vector = create_inputfilename_vector (argument_number, argv);

	if (argument_number > 1) batch = true;

	else {

		inputrgfname = inputfilename();
		inputfilename_vector.push_back(inputrgfname);
	}

	inputfilename_vector = check_rgf_inputs (inputfilename_vector, batch);

	if (batch) using_xy_files = true;
	else using_xy_files = needxyfile ();

	starttime = clock ();

	do {

		inset = manage_settings (batch, inputfilename_vector[j]);

		if (using_xy_files) xy_filename = check_xy_inputs (inputfilename_vector[j], batch);

		process_rgf (inputfilename_vector[j], xy_filename, inset);

		cout << endl;

		cout << "EVALUATION OF " << capslock(inputfilename_vector[j]) << ".RGF FILE COMPLETED." << endl;

		j++;

	} while (j < inputfilename_vector.size());

	finishtime = clock();

	elapsed_time = finishtime - starttime;

	output_elapsed_time (elapsed_time);

	return 0;
}
