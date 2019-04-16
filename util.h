//////////////////////////////////////////////////////////////////////
// util.h  Copyright (c) 2018 Dario Ghersi and Sean West            //
// Version: 20181223                                                //
// Goal: Survival analysis with the minimum p-value method and      //
//       empirically estimated null distribution                    //
//                                                                  //
// This file is part of the NEEP suite.                             //
// NEEP is free software: you can redistribute it and/or            //
// modify it under the terms of the GNU General Public License as   //
// published by the Free Software Foundation, either version 3 of   //
// the License, or (at your option) any later version.              //
//                                                                  //
// NEEP is distributed in the hope that it will be useful,          //
// but WITHOUT ANY WARRANTY; without even the implied warranty of   //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    //
// GNU General Public License for more details.                     //
//                                                                  //
// You should have received a copy of the GNU General Public        //
// License along with NEEP.                                         //
// If not, see <http://www.gnu.org/licenses/>.                      //
//////////////////////////////////////////////////////////////////////

#ifndef _util_
#define _util_

#include "neep.h"

#define MANTEL 1 // Mantel-Cox test
//#define EXPR_THRESH 0.85 // at least x% of transcripts have to be expressed

#define USAGE "\nUsage: neep -c CLINICAL -e EXPRESSION -O OUTPUT -n NUM_ITERATIONS -t EXP_THRESHOLD (-u)\n\n"

struct LrResult
{
	double stat;
	string direction;
	double hr;	 // Hazard ratio
	double mr1y; // Mortality ratio at year 1 (365 days)
	double mr2y; // Mortality ratio at year 2 (730 days)
	double mr5y; // Mortality ratio at year 5 (1825 days)
};

typedef pair<unsigned int, double> intDouble;

//////////////////////////////////////////////////////////////////////
// PROTOTYPES                                                       //
//////////////////////////////////////////////////////////////////////

void checkCommandLineArgs(char **, int);
bool cmdOptionExists(char **, char **, const string&);
vector<unsigned int> createFailureTimes(vector<unsigned int>,
					vector<bool>,
					vector<unsigned int>,
					vector<bool>);
bool comparator(const intDouble &, const intDouble &);
char *getCmdOption(char **, char **, const string &);
LrResult logrank(vector<unsigned int> &, vector<bool> &,
	       vector<unsigned int> &, vector<bool> &);
void printProgBar(unsigned int);
void storeClinicalData(vector<ClinicalSample> &, string);
void storeExpression(vector<ClinicalSample> &, vector<ExpressionData> &,
                     vector<unsigned int> &, string, double);

#endif
