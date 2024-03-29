//////////////////////////////////////////////////////////////////////
// util.C  Copyright (c) 2018 Dario Ghersi and Sean West            //
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

#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std;

#include "neep.h"
#include "util.h"


void checkCommandLineArgs(char **argv, int argc)
{
  // check all the parameters have been provided

  bool err = false;

  if (!cmdOptionExists(argv, argv+argc, "-c")) {
    cerr << "Clinical file missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-e")) {
    cerr << "Expression file missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-o")) {
    cerr << "Output file missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-n")) {
    cerr << "Number of randomizations missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-t")) {
	cerr << "Expression threshold missing\n";
	err = true;
  }
 
  if (err) {
    cout << USAGE;
    exit(1);
  }
}

//////////////////////////////////////////////////////////////////////

bool cmdOptionExists(char **begin, char **end, const string& option)
{
  return find(begin, end, option) != end;
}

//////////////////////////////////////////////////////////////////////

bool comparator(const intDouble &pair1, const intDouble &pair2)
{
  // comparison function

  return pair1.second < pair2.second;
}


//////////////////////////////////////////////////////////////////////

vector<unsigned int> createFailureTimes(vector<unsigned int> timesA,
					vector<bool> eventA,
					vector<unsigned int> timesB,
					vector<bool> eventB)
{
 // create a vector of failure times
  vector<unsigned int> failures;
  for (unsigned int i = 0; i < timesA.size(); i++) {
    if (eventA[i]) {
      failures.push_back(timesA[i]);
    }
  }

  for (unsigned int i = 0; i < timesB.size(); i++) {
    if (eventB[i]) {
      failures.push_back(timesB[i]);
    }
  }
  
  // remove duplicates and sort
  sort(failures.begin(), failures.end());
  failures.erase(unique(failures.begin(), failures.end()),
		 failures.end());

  return failures;
}

//////////////////////////////////////////////////////////////////////

char *getCmdOption(char **begin, char **end, const string & option)
{
  char **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) {
    return *itr;
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////

LrResult logrank(vector<unsigned int> &timesA, vector<bool> &eventA,
	       vector<unsigned int> &timesB, vector<bool> &eventB)
{
  // calculate the logrank test statistics

  double lrStat = 0.0;
  unsigned int sizeA = timesA.size(), sizeB = timesB.size();

  // sort the times and their indices for group A
  vector<intDouble> pPairsA;
  intDouble foo;
  for (unsigned int i = 0; i < sizeA; i++) {
    foo = make_pair(i, timesA[i]);
    pPairsA.push_back(foo);
  }
  sort(pPairsA.begin(), pPairsA.end(), comparator);

  // sort the times and their indices for group B
  vector<intDouble> pPairsB;
  for (unsigned int i = 0; i < sizeB; i++) {
    foo = make_pair(i, timesB[i]);
    pPairsB.push_back(foo);
  }
  sort(pPairsB.begin(), pPairsB.end(), comparator);

  // create a vector of failure times
  vector<unsigned int> failures = createFailureTimes(timesA, eventA,
						     timesB, eventB);
  unsigned int fSize = failures.size();

  // compute the observed and expected events
  double obsA = 0.0, obsB = 0.0, atRiskA = sizeA, atRiskB = sizeB;
  unsigned int indexA = 0.0, indexB = 0.0;
  double currAtRiskA = atRiskA, currAtRiskB = atRiskB, totFailures = 0.0;
  unsigned int f;
  double totObsA = 0.0, totObsB = 0.0, totAtRisk = 0.0;
  double expA = 0.0, expB = 0.0;
  double totExpA = 0.0, totExpB = 0.0;
  double variance = 0.0, num = 0.0; // Mantel-Cox
  for (unsigned int i = 0; i < fSize; i++) {

    f = failures[i];
    currAtRiskA = atRiskA; currAtRiskB = atRiskB;
    totFailures = 0.0, totAtRisk = 0.0;

    // handle group A
    while (pPairsA[indexA].second <= f && indexA < sizeA) {
      if (eventA[pPairsA[indexA].first]) {
	      obsA++;
      }
      currAtRiskA--;
      indexA++;
    }

    if (i < (fSize - 1)) {
      while (pPairsA[indexA].second < failures[i + 1] && indexA < sizeA) {
        indexA++;
        currAtRiskA--;
      }
    }

    // handle group B
    while (pPairsB[indexB].second <= f && indexB < sizeB) {
      if (eventB[pPairsB[indexB].first]) {
        obsB++;
      }
      currAtRiskB--;
      indexB++;
    }

    // handle the situation where the patient with the latest followup did not survive
    // This bit of code may be producing errors
    /*
    if (i < (fSize - 1)) {
      while (pPairsB[indexB].second < failures[i + 1] && indexB < sizeB) {
        indexB++;
        currAtRiskB--;
      }
    }
    */
    // It needs more testing and failure-proofing

    // calculate the total failures
    totFailures = obsA + obsB;

    // calculate the expected cases
    expA = totFailures * atRiskA / (atRiskA + atRiskB);
    expB = totFailures - expA;

    // update the totals
    totObsA += obsA; totObsB += obsB;
    totExpA += expA;
    totExpB += expB;

    if (MANTEL) { // calculate the variance
      totAtRisk = atRiskA + atRiskB;
      if (totAtRisk > 1){
	    variance += totFailures * (atRiskA / totAtRisk) * (1.0 - atRiskA / totAtRisk) *
	  			    (totAtRisk - totFailures) / (totAtRisk - 1.0);
	    num += obsA - expA;
      }

    }

    // reset the observed cases
    obsA = obsB = 0; totAtRisk = 0;
    atRiskA = currAtRiskA;
    atRiskB = currAtRiskB;
  }

  // calculate the hazard ratio
  double l = (totObsA - totExpA) / variance;

  // calculate the mortality ratios
  double eventsA1y = 0;
  double eventsA2y = 0;
  double eventsA5y = 0;
  for (unsigned int i = 0; i < sizeA; i++){
    if (pPairsA[i].second < 1825){
    	eventsA5y += 1.0;
    }
    if (pPairsA[i].second < 730){
    	eventsA2y += 1.0;
    }
    if (pPairsA[i].second < 365){
    	eventsA1y += 1.0;
    }
  }
  double eventsB1y = 0;
  double eventsB2y = 0;
  double eventsB5y = 0;
  for (unsigned int i = 0; i < sizeB; i++){
    if (pPairsB[i].second < 1825){
    	eventsB5y += 1.0;
    }
    if (pPairsB[i].second < 730){
    	eventsB2y += 1.0;
    }
    if (pPairsB[i].second < 365){
    	eventsB1y += 1.0;
    }
  }
  // MR 1 year (365 days)
  double mr1y = (eventsB1y * sizeA) / (sizeB * eventsA1y);
  // MR 2 year (730 days)
  double mr2y = (eventsB2y * sizeA) / (sizeB * eventsA2y);
  // MR 3 year (1825 days)
  double mr5y = (eventsB5y * sizeA) / (sizeB * eventsA5y);

  // calculate the chi squared statistics
  if (MANTEL) {
    lrStat = pow(num / sqrt(variance), 2);
  }
  else {
    lrStat = pow(totObsA - totExpA, 2) / totExpA +
             pow(totObsB - totExpB, 2) / totExpB;
  }

  // set direction of survival
  string direction = "high expression survived longer";
  if (totObsA - totExpA < 0){
	  direction = "low expression survived longer";
  }

  LrResult result;
  result.stat = lrStat;
  result.direction = direction;
  result.hr = exp(l);
  result.mr1y = mr1y;
  result.mr2y = mr2y;
  result.mr5y = mr5y;

  return result;
}

//////////////////////////////////////////////////////////////////////

void storeClinicalData(vector<ClinicalSample> &clinical, string fileName)
{
  // store the sample ID, days to event, and event status into
  // a vector of structures ('clinical')

  string line, temp;
  string barcode, s_days, s_event;

  // open the input file
  fstream infile;
  infile.open(fileName, fstream::in);

  // complain if the file doesn't exist
  if (! infile.good()) {
    cerr << "Can't open " << fileName << endl;
    exit(1);
  }

  // process each sample
  while (getline(infile, line)) {

    ClinicalSample cs;

    // split the string into three fields
    stringstream linestream(line);
    getline(linestream, barcode, ',');
    getline(linestream, s_days, ',');
    getline(linestream, s_event, ',');

    // add the sample to the vector
    cs.barcode = barcode;
    cs.days = stoi(s_days);
    cs.event = stoi(s_event);
    clinical.push_back(cs);
  }

  infile.close();
}

//////////////////////////////////////////////////////////////////////

void printProgBar(unsigned int percent) {
  // print a progress bar

  string bar;

  for (unsigned int i = 0; i < 50; i++) {
    if (i < (percent / 2)) {
      bar.replace(i, 1, "=");
    }
    else if (i == (percent / 2)) {
      bar.replace(i, 1, ">");
    }
    else {
      bar.replace(i, 1, " ");
    }
  }

  cout << "\r" "[" << bar << "] ";
  cout.width(3);
  cout << percent << "%     " << flush;
}

//////////////////////////////////////////////////////////////////////

void storeExpression(vector<ClinicalSample> &clinical,
                     vector<ExpressionData> &expression,
                     vector<unsigned int> &index,
					 string fileName,
					 double expressionThreshold)
{
  // store the expression values of transcripts (or genes) that
  // have at least 1.0-expressionThreshold expression in all clinical samples

  string line, temp, id;

    // open the input file
  fstream infile;
  infile.open(fileName, fstream::in);

  // complain if the file doesn't exist
  if (! infile.good()) {
    cerr << "Can't open " << fileName << endl;
    exit(1);
  }

  // process the header and extract the indices of the matching
  // clinical samples
  unsigned int numSamples = clinical.size();
  getline(infile, line);
  stringstream linestream(line);
  getline(linestream, temp, ',');
  for (unsigned int i = 0; i < numSamples; i++) {
    getline(linestream, temp, ',');
    for (unsigned int j  = 0; j < numSamples; j++) {
      if (temp == clinical[j].barcode) {
        index.push_back(j);
      }
    }
  }

  // process each gene/isoform
  unsigned int isExpr = 0; // no. of transcripts/genes with non-zero expr.
  while (getline(infile, line)) {
    vector<double> expr;
    isExpr = 0;

    // get the sequence id
    stringstream linestream(line);
    getline(linestream, id, ',');

    // get all the expression values
    for (unsigned int i = 0; i < numSamples; i++) {
      double value;
      string s_value;
      getline(linestream, s_value, ',');
      value = stod(s_value);
      if (value > 0) {
      	isExpr += 1;
      }
      expr.push_back(value);
    }

    // if the fraction of expressed samples is >= threshold
    if (double(isExpr) / numSamples >= (1.0 - expressionThreshold)) {
      struct ExpressionData ed;
      ed.id = id;
      ed.exprVect = expr;
      expression.push_back(ed);
    }
  }

  infile.close();
}
