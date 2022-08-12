#ifndef FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H
#define FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H

#include <vector>
#include <random>
#include <set>
#include <cmath>
#include <algorithm>
#include <boost/math/special_functions/digamma.hpp>

using namespace std;

class EsRuler {
private:

    const vector<double> &ranks;
    const unsigned int sampleSize;
    const unsigned int pathwaySize;


    vector<double> enrichmentScores;
    vector<vector<int> > currentSamples;

    vector<unsigned int> probCorrector;

    void duplicateSamples();
    void updateSample();

    // vector<int> chunkLastElement;
    // int chunksNumber;

    // struct SampleChunks {
    //     vector<double> chunkSum;
    //     vector<vector<int>> chunks;
    //     SampleChunks(int);
    // };

    // int perturbate(const vector<double> &ranks, int k, SampleChunks &cusSampleChunks,
    //            double bound, mt19937 &rng);

    // int chunkLen(int ind);

public:

    EsRuler(const vector<double> &inpRanks, unsigned int inpSampleSize, unsigned int inpPathwaySize);

    ~EsRuler();

    void extend(double ES, int seed, double eps);

    pair<double, bool> getPvalue(double ES, double eps, bool sign);
};

double betaMeanLog(unsigned long a, unsigned long b);

pair<double, bool> calcLogCorrection(const vector<unsigned int> &probCorrector,
                                     long probCorrIndx, unsigned int sampleSize);

int perturbate(const vector<double> &ranks, vector<int> &sample, double bound, mt19937 &rng);


#endif //FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H
