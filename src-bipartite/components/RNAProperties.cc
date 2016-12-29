#include <iostream>
#include <sstream>
#include <fstream>
#include "Weights.h"
#include "RNAProperties.h"
#include "../utils/Arrays.h"
using std::string;

RNAProperties::RNAProperties(Weights * weights, vector<OrigRNASeq> * v) :
	weights(weights), origRNASequences(v)
{
	int numOfRNA = origRNASequences->size();
    // int xi;

    int i;
    short evenSeq[10], e = 0;
    short oddSeq[10], o = 0;

    
    // find even or odd
    for(i=0;i<numOfRNA;i++)
    {
        if((origRNASequences->at(i)).type == 0)
            evenSeq[e++] = i;
        else if((origRNASequences->at(i)).type == 1)
            oddSeq[o++] = i;    
    }

    //do RNAup for each even/odd pair
    int pairs = o * e;

    numEven = e;
    numOdd = o;

    // fill up rnaLengths
    // first fill up evens, then odds
    for(int eId=0; eId<e; eId++)
    {
    	int id = evenSeq[eId];
    	cout << "in even, id = " << id << endl;
    	rnaLengths.push_back(origRNASequences->at(id).originalLength); 
    }

    for(int oId=0; oId<o; oId++)
    {
    	int id = oddSeq[oId];
    	cout << "in odd, id = " << id << endl;
    	rnaLengths.push_back(origRNASequences->at(id).originalLength); 
    }

}


int RNAProperties::getNumEven()
{
    return numEven;
}

int RNAProperties::getNumOdd()
{
    return numOdd;
}

Weights * RNAProperties::getWeights()
{
	return weights;
}