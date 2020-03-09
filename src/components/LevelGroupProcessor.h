#ifndef _LEVEL_GROUP_PROCESSOR_H
#define _LEVEL_GROUP_PROCESSOR_H

#include <string>
#include <vector>
#include <vector>
#include "rnaseq.h"
#include "Chain.h"
#include "Config.h"
#include "Weights.h"
using namespace std;

// A processor to find the optimal structure in a LevelGroup. A Level Group is a set of rnas 
// stacked in some order.
class LevelGroupProcessor
{
private:
    // k (1/eps = max num of levels per group)
    int k;
    Chain * chain;
    vector<vector<int> > allGroupings;
    Weights* weights;
    bool loopAround;

public:
    LevelGroupProcessor(int k, Chain*, Weights*, bool loopAround);

    vector<pair<int,int>> prepareWithSubsets(int groupingNum, Config&);
    vector<vector<int>> GetAllGroupings();
    Config findMinEnergyStructure();
};

#endif
