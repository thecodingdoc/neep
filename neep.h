//////////////////////////////////////////////////////////////////////
// neep.h  Copyright (c) 2018 Dario Ghersi and Sean West            //
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

#ifndef _neep_
#define _neep_

///////////////////////////////////////////////////////////////////////
// STRUCTURES                                                        //
///////////////////////////////////////////////////////////////////////

struct ClinicalSample {
  string barcode;
  unsigned int days;
  bool event;
};

//////////////////////////////////////////////////////////////////////

struct ExpressionData {
  string id;
  vector<double> exprVect;
};

///////////////////////////////////////////////////////////////////////

struct BestLogRank {
  double stat;
  unsigned int bestPos;
  string direction;
};

///////////////////////////////////////////////////////////////////////
// CLASSES                                                          //
//////////////////////////////////////////////////////////////////////

class Parameters {

 public:
  string clinicalFileName;
  string expressionFileName;
  string outFileName;
  unsigned int numIter;

  Parameters(char **, int);
};

//////////////////////////////////////////////////////////////////////

class Results {

 public:
  vector<double> rawP;
  vector<double> adjustedP;
  vector<unsigned int> sortedOrder;

  Results(vector<double> &);
};

#endif

//////////////////////////////////////////////////////////////////////
// PROTOTYPES                                                       //
//////////////////////////////////////////////////////////////////////

void calculateBestLogRank(vector<ExpressionData> &,
                          vector<ClinicalSample> &,
                          vector<unsigned int> &,
                          vector<BestLogRank> &);
void calculateNull(vector<ClinicalSample> &, vector<double> &,
                   unsigned int);
void calculatePValues(vector<double> &, vector<BestLogRank> &,
                      vector<double> &);
void printResults(string, Results &, vector<BestLogRank> &);
