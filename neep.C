//////////////////////////////////////////////////////////////////////
// neep.C  Copyright (c) 2018 Dario Ghersi and Sean West            //
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
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>

using namespace std;

#include "util.h"
#include "neep.h"

double epsilon = numeric_limits<double>::epsilon();

///////////////////////////////////////////////////////////////////////
// CONSTRUCTORS                                                      //
///////////////////////////////////////////////////////////////////////

Parameters::Parameters(char **argv, int argc)
{
  // parse the command-line arguments

  clinicalFileName = getCmdOption(argv, argv + argc, "-c");
  expressionFileName = getCmdOption(argv, argv + argc, "-e");
  outFileName = getCmdOption(argv, argv + argc, "-o");
  numIter = stoi(getCmdOption(argv, argv + argc, "-n"));
  expressionThreshold = stod(getCmdOption(argv, argv + argc, "-t"));
  isUniform = cmdOptionExists(argv, argv + argc, "-u");
}

//////////////////////////////////////////////////////////////////////

Results::Results(vector<double> &pvalues)
{
  // perform FDR correction using the Benjamini-Hochberg approach

  // combine the p-values wih their index for sorting
  vector<intDouble> pPairs;
  intDouble foo;
  for (unsigned int i = 0; i < pvalues.size(); i++) {
    foo = make_pair(i, pvalues[i]);
    pPairs.push_back(foo);
  }

  // sort the p-values
  sort(pPairs.begin(), pPairs.end(), comparator);
  
  // apply the Benjamini-Hochberg procedure to correct the p-values
  // for multiple hypothesis testing
  double currMin, value;
  vector<double> minValues;
  unsigned int m = pPairs.size();
  for (unsigned int i = 0; i < m; i++) {
    currMin = pPairs[i].second * m / (i + 1);
    for (unsigned int j = i + 1; j < m; j++) {
      value = pPairs[j].second * m / (j + 1);
      if (value < currMin) {
        currMin = value;
      }
    }
    minValues.push_back(min(currMin, 1.0));
  }

  // put the adjusted values in the original p-value order
  adjustedP.resize(pPairs.size());
  sortedOrder.resize(pPairs.size());
  rawP.resize(pPairs.size());

  for (unsigned int i = 0; i < pPairs.size(); i++) {
    adjustedP[pPairs[i].first] = minValues[i];
    rawP[pPairs[i].first] = pPairs[i].second;
    sortedOrder[i] = pPairs[i].first;
  }
}

//////////////////////////////////////////////////////////////////////
// FUNCTIONS                                                        //
//////////////////////////////////////////////////////////////////////

void calculateBestLogRank(vector<ExpressionData> &expression,
                          vector<ClinicalSample> &clinical,
                          vector<unsigned int> &index,
                          vector<BestLogRank> &bestLogRank,
						  double expressionThreshold)
{
  // process each transcript/gene and return the best logrank
  // statistics and the corresponding best split

  unsigned int numSamples = clinical.size();
  unsigned int numExpr = expression.size();

  // calculate the minimum and maximum indices
  unsigned int minInd, maxInd;
  minInd = floor(expressionThreshold * numSamples);
  maxInd = floor((1.0 - expressionThreshold) * numSamples);

  // process each transcript (or gene)
  double oldPercentage = 0.0;
  printProgBar(0.0);
  for (unsigned int i = 0; i < numExpr; i++) {
	//cout << "\n\n" << i << "\n";
	//cout << index[i] << "\n";
	//cout << expression[i].id << "\n";



    // call the progress bar every 100 transcripts
    if ((i % 100) == 0 || numExpr < 100) {
      double percentage = 100.0 * i / numExpr;
      if (percentage > oldPercentage) {
        oldPercentage = percentage;
        printProgBar(percentage);
      }
    }
    
    // sort the expression vector
    vector<intDouble> pPairs;
    intDouble foo;
    for (unsigned int j = 0; j < numSamples; j++) {
      foo = make_pair(j, expression[i].exprVect[j]);
      pPairs.push_back(foo);
    }
    sort(pPairs.begin(), pPairs.end(), comparator);

   // identify the best threshold and corresponding max statistics
   double maxStat = 0.0;
   unsigned bestPos = 0;
   string maxDir = "";
   double hr = 0.0;
   double mr1y = 0.0;
   double mr2y = 0.0;
   double mr5y = 0.0;
   
   //for (unsigned int currPos = minInd; currPos < maxInd; currPos++) {
   unsigned int currPos = minInd;
   while (currPos <= maxInd) {
     vector<unsigned int> timesA, timesB;
     vector<bool> eventA, eventB;

     // take care of ties
     while (fabs(expression[i].exprVect[index[pPairs[currPos].first]] -
                expression[i].exprVect[index[pPairs[currPos + 1].first]]) <
            epsilon && currPos <= maxInd) {
       currPos++;
     }

     // populate group A
     for (unsigned int j = 0; j < currPos; j++) {
       timesA.push_back(clinical[index[pPairs[j].first]].days);
       eventA.push_back(clinical[index[pPairs[j].first]].event);
     }

     // populate group B
     for (unsigned int j = currPos; j < numSamples; j++) {
       timesB.push_back(clinical[index[pPairs[j].first]].days);
       eventB.push_back(clinical[index[pPairs[j].first]].event);
     }


     // calculate the logrank statistics
     LrResult result = logrank(timesA, eventA, timesB, eventB);

     if (result.stat > maxStat) {
       maxStat = result.stat;
       bestPos = currPos;
       maxDir = result.direction;
       hr = result.hr;
       mr1y = result.mr1y;
       mr2y = result.mr2y;
       mr5y = result.mr5y;
     }

     currPos++;
   }
   
   struct BestLogRank blogrank;
   blogrank.stat = maxStat;
   blogrank.bestPos = bestPos;
   blogrank.direction = maxDir;
   blogrank.hr = hr;
   blogrank.mr1y = mr1y;
   blogrank.mr2y = mr2y;
   blogrank.mr5y = mr5y;
   bestLogRank[i] = blogrank;

  }


  printProgBar(100.0);

}

//////////////////////////////////////////////////////////////////////

void calculateNull(vector<ClinicalSample> &clinical,
                   vector<double> &nullDist,
				   unsigned int numIter,
				   double expressionThreshold,
				   bool isUniform)
{
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);

  // calculate the null distribution for the minimum p-value stat
  unsigned int numSamples = clinical.size();
  unsigned int minInd, maxInd;
  minInd = floor(expressionThreshold * numSamples);
  maxInd = floor((1.0 - expressionThreshold) * numSamples);

  double oldPercentage = 0.0;
  unsigned itCompleted = 0;
  printProgBar(0.0);
  #pragma omp parallel for
  for (unsigned int it = 0; it < numIter; it++) {

    // call the progress bar every 100 iterations
    if ((itCompleted % 500) == 0 || numIter < 500) {
      double percentage = 100.0 * itCompleted / numIter;
      if (percentage > oldPercentage) {
        oldPercentage = percentage;
        printProgBar(percentage);
      }
    }
    
    vector<unsigned int> myV(numSamples);

    // construct the randomized int vector myV
    if (isUniform){
      // expression calculation
      vector<double> expression(numSamples);
      for (unsigned int m = 0; m < numSamples; m++){
    	  expression[m] = uniform(generator);
      }

      // sort spots
      vector<intDouble> pPairs;
      intDouble foo;
      for (unsigned int j = 0; j < numSamples; j++) {
        foo = make_pair(j, expression[j]);
        pPairs.push_back(foo);
      }
      sort(pPairs.begin(), pPairs.end(), comparator);


      for (unsigned int m = 0; m < numSamples; m++){
    	  myV[m] = pPairs[m].first;
      }
    } else {
	  // initialize the vector of indices to use in the permutations
	  for (unsigned int i = 0; i < numSamples; i++) {
	    myV[i] = i;
  	  }
      // shuffle the vector of indices
      random_shuffle(myV.begin(), myV.end());
    }


    // process each threshold
    double maxStat = 0.0;

    //    printProgBar(0.0);
    for (unsigned int currPos = minInd; currPos <= maxInd; currPos++) {
      vector<unsigned int> timesA, timesB;
      vector<bool> eventA, eventB;

      // populate group A
      for (unsigned int j = 0; j < currPos; j++) {
        timesA.push_back(clinical[myV[j]].days);
        eventA.push_back(clinical[myV[j]].event);
      }

      // populate group B
      for (unsigned int j = currPos; j < numSamples; j++) {
        timesB.push_back(clinical[myV[j]].days);
        eventB.push_back(clinical[myV[j]].event);
      }

      // calculate the logrank statistics
      LrResult result = logrank(timesA, eventA, timesB, eventB);
      if (result.stat > maxStat) {
        maxStat = result.stat;
      }
    }
    nullDist[it] = maxStat;

    #pragma omp atomic
    itCompleted++;
  }

  printProgBar(100.0);
}

//////////////////////////////////////////////////////////////////////

void calculatePValues(vector<double> &nullDist,
                      vector<BestLogRank> &bestLogRank,
                      vector<double> &empiricalP)
{
  // calculate the empirical p-values for each transcript (or gene)

  unsigned int numPerm = nullDist.size();
  
  // sort the null distribution
  sort(nullDist.begin(), nullDist.end());

  // calculate the empirical p-value
  for (unsigned int i = 0; i < bestLogRank.size(); i++) {
    for (unsigned int j = 0; j < numPerm; j++) {
      if (bestLogRank[i].stat <= nullDist[j]) {
        empiricalP[i] = 1.0 - double(j) / numPerm;
        break; // no need to check further
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void printResults(string outFileName, Results &results,
		  vector<ExpressionData> &expression,
		  vector<BestLogRank> &bestLogRank)
{
  // print the results

  fstream outFile;
  outFile.open(outFileName, fstream::out);
  
  // print the header
  outFile << "ID\tBEST_STATISTICS\tBEST_SPLIT\tP-VALUE\tFDR\tDIRECTION\tHAZARD_RATIO\tMORTALITY_RATIO_1YR\tMORTALITY_RATIO_2YR\tMORTALITY_RATIO_5YR\n";

  // print all calculated values
  for (unsigned int i = 0; i < results.sortedOrder.size(); i++) {
    outFile << expression[results.sortedOrder[i]].id << "\t"
    		<< fixed << bestLogRank[results.sortedOrder[i]].stat << "\t"
			<< bestLogRank[results.sortedOrder[i]].bestPos << "\t"
			<< scientific << results.rawP[results.sortedOrder[i]] << "\t"
			<< results.adjustedP[results.sortedOrder[i]] << "\t"
			<< bestLogRank[results.sortedOrder[i]].direction << "\t"
			<< bestLogRank[results.sortedOrder[i]].hr << "\t"
			<< bestLogRank[results.sortedOrder[i]].mr1y << "\t"
			<< bestLogRank[results.sortedOrder[i]].mr2y << "\t"
			<< bestLogRank[results.sortedOrder[i]].mr5y << endl;
  }

  outFile.close();
}

//////////////////////////////////////////////////////////////////////
// MAIN PROGRAM                                                     //
//////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // check the command-line arguments
  checkCommandLineArgs(argv, argc);

  // get the parameters
  Parameters p(argv, argc);

  // store the clinical variables
  cout << "Storing clinical variables..." << flush;
  vector<ClinicalSample> clinical;
  storeClinicalData(clinical, p.clinicalFileName);
  cout << "done\n" << flush;

  // store the expression data
  cout << "Storing expression data..." << flush;
  vector<ExpressionData> expression;
  vector<unsigned int> index; // the index of the clinical samples
                              // matching the expression samples
  storeExpression(clinical, expression, index, p.expressionFileName, p.expressionThreshold);
  cout << "done\n" << flush;

  // calculate the best logrank statistics and corresponding split
  // for each isoform
  vector<BestLogRank> bestLogRank(expression.size());
  cout << "Computing the minimum p-value for each threshold...\n" << flush;
  calculateBestLogRank(expression, clinical, index, bestLogRank, p.expressionThreshold);
  cout << endl << flush;

  // calculate the null distribution
  vector<double> nullDist(p.numIter);
  cout << "Computing the null distribution...\n" << flush;
  calculateNull(clinical, nullDist, p.numIter, p.expressionThreshold, p.isUniform);
  cout << endl << flush;

  // calculate the empirical p-values
  cout << "Calculating the empirical p-values..." << flush;
  vector<double> empiricalP(bestLogRank.size());
  calculatePValues(nullDist, bestLogRank, empiricalP);
  cout << "done\n" << flush;

  // perform p-value adjustment
  cout << "Performing FDR correction..." << flush;
  Results results(empiricalP);
  cout << "done\n" << flush;
  
  // print the results
  printResults(p.outFileName, results, expression, bestLogRank);

  return 0;
}
