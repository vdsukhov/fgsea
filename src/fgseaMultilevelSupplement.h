#ifndef FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H
#define FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H

#include <vector>
#include <random>
#include <set>
#include <cmath>
#include <algorithm>
#include <boost/math/special_functions/digamma.hpp>

#include "esCalculation.h"
#include "util.h"
#include "Rcpp.h"

using namespace std;

class EsRuler {
private:

    const vector<double> &ranks;
    const unsigned int sampleSize;
    const unsigned int pathwaySize;
    const double movesScale;
    const bool logStatus;


    vector<double> enrichmentScores;
    vector<vector<int> > currentSamples;

    vector<unsigned int> probCorrector;

    void duplicateSamples(mt19937 &rng);

    vector<int> chunkLastElement;
    int chunksNumber;

    struct SampleChunks {
        vector<double> chunkSum;
        vector<vector<int>> chunks;
        SampleChunks(int);
    };

    int perturbate(const vector<double> &ranks, int k, SampleChunks &cusSampleChunks,
               double bound, mt19937 &rng);

    int chunkLen(int ind);

public:

    EsRuler(const vector<double> &inpRanks,
            unsigned int inpSampleSize,
            unsigned int inpPathwaySize,
            double inpMovesScale,
            bool inpLog);

    ~EsRuler();

    void extend(double ES, int seed, double eps);

    pair<double, bool> getPvalue(double ES, double eps, bool sign);
};

double betaMeanLog(unsigned long a, unsigned long b);

pair<double, bool> calcLogCorrection(const vector<unsigned int> &probCorrector,
                                     long probCorrIndx, unsigned int sampleSize);


#endif //FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H
