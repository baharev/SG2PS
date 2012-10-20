#ifndef CHECKSETTINGFILECONTENT_HPP_
#define CHECKSETTINGFILECONTENT_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include "structs.h"
#include "common.h"

#include "checksettingfilecontent.h"

using namespace std;

void header ();
bool settingfilecorrect (string settingfilename);
INPSET_ loadsettingsfromsettingfile (string settingfilename);
void printsettingsonscreen (INPSET_ settings);
INPSET_ inputsettings_manually (string projectname);
bool outputsettingfile (INPSET_ _outputsettingfile, string projectname);
INPSET_ input_hardcoded ();
INPSET_ manage_settings (bool batch, string projectname);
string input_setting_decision ();

#endif
