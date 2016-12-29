// Originally written in C, with a major rewrite in C++, 
// with a lot of the old C code still thrown around.
// Sorry ^_^  !

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cstring>
#include <string>
#include <float.h> 
#include <cstdarg>
#include <limits>
#include <cstddef>
#include <fstream>
#include "PRBDPCore.h"
#include "RNAProperties.h"
#include "rnaseq.h"
#include "Window.h"
#include "Weights.h"
#include "Config.h"
#include "../utils/MultiIndexArray.h"

using namespace std;

#define USE_CACHE 0
#define TEST 0
#define TEST_TO_FILE 0


PRBDPCore::PRBDPCore(RNAProperties * rnaProps)
{
    weights = rnaProps->getWeights();

    numEven = rnaProps->getNumEven();
    numOdd = rnaProps->getNumOdd();

    cout << "New constructor" << endl;
    totalLevels = numEven + numOdd;

    cout << "total levels = " << totalLevels << endl;
    cout << "numEven = " << numEven << endl;
    cout << "numOdd = " << numOdd << endl;

    realLengths.resize(totalLevels);
    lengths.resize(totalLevels);
    rnaLengths.resize(totalLevels);

// cout << "S" << endl;
    // Recall that order in rnaProps->rnaLengths is all evens then all odds
    for(int i=0;i<totalLevels;i++)
    {
        rnaLengths[i] = rnaProps->rnaLengths[i];
        lengths[i] = rnaProps->rnaLengths[i] + 1;
    }
        

    initLinear(totalLevels);
    allForLoops();

}
   

PRBDPCore::PRBDPCore(SingleRunConfig * data, int gapSize, int wType, bool wrapAround) 
    : data(data), wType(wType), wrapAround(wrapAround)
{

#if TEST_TO_FILE == 1
    testfp = fopen("testfp.txt", "w");
#endif

    double temp_double;
    double * returnVal = &temp_double;
    gap = gapSize;

    int num = data->numOfRNA;
    
    for(int i=0; i<10*10; i++)
        interLocs[0][i] = data->interLocs[0][i];
    
    totalLevels = data->totalLevels;
    
    // dimension are given by number of strings numOfRNA
    realLengths.resize(totalLevels);
    lengths.resize(totalLevels);
    cutoffs.resize(num+1);
    treeNodeLens.resize(num+1); 
    // elemToRNAmap.resize(totalLevels);
    rnaLengths.resize(num+1);   

    for(int i=0;i<num;i++)
        rnaLengths[i] = data->rnaSequences[i]->originalLength;

    treeNodeLens[0] = 0;

    initMetaTrees(num);
    initLinear(num);


    allForLoops();
        
#if TEST_TO_FILE == 1
            fclose(testfp);
#endif
}





void PRBDPCore::initMetaTrees(int num)
{

#if TEST == 1
    printf("data->setup[i][j] looks like: \n");
    for(int xa=0;xa<10;xa++) {
        for(int xb=0;xb<10;xb++) {
            printf("%d ", data->setup[xa][xb]);
        }
        printf("\n");
    }
#endif  

    // setup the cuttoffs array. Cutoffs specify the cumulative length of all RNAs prior to this one
    int i; 
    for(i=0;i<totalLevels;i++)
    {
        int j;
        //if there is no node here, move on
        if(data->setup[i][0] == -9) continue;
        
        //the first element of each row is zero. This is because the cutoff the first rna on a level will be 0 (no previous RNA!)
        cutoffs[data->setup[i][0]] = 0;
        
        //look at all RNAs on this level
        for(j=1;j<5 && data->setup[i][j] != -9;j++)
        {   
            int rna = data->setup[i][j];    //get RNA Id
            int prevRNA = data->setup[i][j-1]; // Id of RNA behind this one

            // set cutoff this RNA to length of last RNA + cutoff of that RNA
            // cutoffs[rna] = strlen(str1[prevRNA - 1]) + cutoffs[prevRNA]; 
            cutoffs[rna] = data->rnaSequences[prevRNA]->compressedLength + cutoffs[prevRNA];    
        }
    }

    // to be safe, set cutoff of Null RNA to 0
    cutoffs[0] = 0;
    
    // does 3 things: (i = level #)
    // 1. calculates the lengths of entire levels
    // 2. saves lengths of each node/RNA in treeNodeLens
    // 3. fills up elemToRNAmap
    for(i=0;i<totalLevels;i++)
        realLengths[i] = 0;
    
    for(i=0;i<totalLevels;i++)
    {
        int j;      //column iterator
        int numInThisLevel = 0;     //counter of how many RNAs we have seen
        
        for(j=0;j<5 && data->setup[i][j] != -9;j++)
        {
            numInThisLevel++;   //update counter
            int rna = data->setup[i][j];    //get the Id of the RNA we are looking at
            int len = data->rnaSequences[rna]->compressedLength;    //get length of this RNA, from string table
        
            // update length of the whole level 
            realLengths[i] += len;
        
            // save the length of this RNA/node
            treeNodeLens[rna] = len;
        }
        
        // Why are we doing this?
        lengths[i] = realLengths[i] + 1;
        
        // initialize the array for this level
        elemToRNAmap[i].resize(lengths[i]); // = (int*) malloc(sizeof(int) * lengths[i]);
        
        // set the 0th elem to 0, because first elem of level is 1
        elemToRNAmap[i][0] = 0;
        
        int count = 1;  //This iterator iterates over the columns in this row
        
        // for all RNAs in this level
        for(j=0;j<numInThisLevel;j++)
        {
            int rnaNum = data->setup[i][j]; //get Id of RNA given by level, j.
            
            // for all positions occupied by this RNA
            int k; for(k=0;k<treeNodeLens[rnaNum];k++)
            {
                // map this elem to rnaNum, inc elem #
                elemToRNAmap[i][count++] = rnaNum;
            }
        }
        
    }
}

// All stuff related to linearization initialization
void PRBDPCore::initLinear(int num)
{
    // find out the total number of elements the grand array is going to hold
    // this is = len(level1) * ... * len(level_k)
    totalDimSize = 1;
    for(int i=0;i<totalLevels;i++)
        totalDimSize *= lengths[i];
    
    // make hash multipliers, i.e., the bases of the hash function. bases[0] = len(Level1), base[1] = base[0]*len(level2), and so on...
    // (ADD more details in hash function comments...)
    bases.resize(totalLevels);
    
    bases[0] = lengths[0];
    for(int i=1;i<totalLevels;i++)
        bases[i] = lengths[i] * bases[i-1];
    
    H = MultiIndexArray<double>::CreateMultiIndexArray(lengths);
    chosenByArr = MultiIndexArray<short>::CreateMultiIndexArray(lengths);
    chosenW = MultiIndexArray<int>::CreateMultiIndexArray(lengths);
    chosenW2 = MultiIndexArray<int>::CreateMultiIndexArray(lengths);
}


/*
void PRBDPCore::allForLoops()
{

    auto begin = chrono::high_resolution_clock::now();  

    H[0] = 0;   // base case
    // start i from where all indices >= 1
    int block = totalDimSize/10, percent = 0;
    int * indices = new int[totalLevels];
    
    printf("Completed: ");
    for(int i=1;i<totalDimSize;i++)
    {
        
        H.revHash(i, indices);
        atStepH(i, indices);
        
        if(i % block == 0) 
        { 
            percent += 10;  
            printf("%d %% .. ", percent); 
            fflush(stdout);
        };

    }
    printf("\n");

#if TEST == 1
    printf("iterCount = %ld\n", iterCount++);
    printf("H_arrayHits = %d\n", H_arrayHits++);
    printf("W_arrayHits = %d\n", W_arrayHits++);
    
#endif
    int * finalIndices = new int[totalLevels];
    
    for(int i=0;i<totalLevels;i++)
        finalIndices[i] = realLengths[i];

    double h = H[finalIndices];
    printf("Sub Total: %f  \n", h);
        
    backtrackNR(finalIndices);      

    delete[] finalIndices;
    delete[] indices;

    auto end = chrono::high_resolution_clock::now();    
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    
    // cout << "Time after all H:" << (double)ms << "ms" << endl;

}
*/


void PRBDPCore::allForLoops()
{

    auto begin = chrono::high_resolution_clock::now();  

    H[0] = 0;   // base case
    // start i from where all indices >= 1
    int block = totalDimSize/10, percent = 0;
    int * indices = new int[totalLevels];
    
    cout << "totalDimSize = " << totalDimSize << endl;
    printf("Completed: ");
    for(int i=1;i<totalDimSize;i++)
    {
        
        H.revHash(i, indices);
        
        atStepH(i, indices);
        
        // if(i % block == 0) 
        // { 
        //     percent += 10;  
        //     printf("%d %% .. ", percent); 
        //     fflush(stdout);
        // };

    }
    // printf("\n");

#if TEST == 1
    printf("iterCount = %ld\n", iterCount++);
    printf("H_arrayHits = %d\n", H_arrayHits++);
    printf("W_arrayHits = %d\n", W_arrayHits++);
    
#endif
    int * finalIndices = new int[totalLevels];
    
    for(int i=0;i<totalLevels;i++)
        // finalIndices[i] = realLengths[i];
        finalIndices[i] = lengths[i];

    cout << "totalDimSize = " << totalDimSize << endl; 
    double h = H[totalDimSize-1];
    printf("Sub Total: %f  \n", h);
        
    // backtrackNR(finalIndices);      

    delete[] finalIndices;
    delete[] indices;

    auto end = chrono::high_resolution_clock::now();    
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    
    cout << "Time after all H:" << (double)ms << "ms" << endl;

}





void PRBDPCore::atStepH(int theIndex, int * indices)
{
    int doPrint = 0;
    
    int ci;
    for(ci=0;ci<totalLevels;ci++)
        if(indices[ci] < 0)
            return;
    
#if TEST == 1     
    for(ci=0;ci<totalLevels;ci++)
        cout << indices[ci] << " ";
    cout << endl;
#endif
    // return;


    int indd[4];
    indd[0] = indices[0];
    indd[1] = indices[1];
    indd[2] = indices[2];
    indd[3] = indices[3];
    
    
    double min = std::numeric_limits<double>::max(), val;
    int argMinI, argMinW, argMinW2;
    int i , chosenCount = 0, chosenBy;

#if TEST == 1 
    printf("Computing for: ");  
    printf("[%d %d %d %d]  \n", indd[0], indd[1], indd[2], indd[3]);
#endif
#if TEST_TO_FILE == 1 
    fprintf(testfp, "Computing for: "); 
    fprintf(testfp, "[%d %d %d %d]  \n", indd[0], indd[1], indd[2], indd[3]);
#endif


    // First part: move one peg on each level: 
    // H(i-1,j,k,..), H(i,j-1,k,...), H(i,j,k-1,...), etc
    for(i=0;i<totalLevels;i++)
    {
        iterCount++;
        chosenCount++;
    
        //for r_i, do r_i = r_i - 1 just for the sake of passing down
        indices[i]--;
        if(indices[i] < 0)
        {
            indices[i]++;
            continue;
        }
        val = H[indices];
        

#if TEST == 1
        H_arrayHits++;
#endif

        if(val < min)
        {
            min = val;
            argMinI = i;
            chosenBy = chosenCount;

#if TEST_TO_FILE == 1           
            fprintf(testfp, "Found min in part 1 \n");
            fprintf(testfp, "[%d %d %d %d]  \n", indd[0], indd[1], indd[2], indd[3]);
            fprintf(testfp, "With min = %f \n", min);
            // fprintf(testfp, "top: %d, bot: %d, H(indices): %f, weight: %f, val: %f, ind[%d]: %d, ind[%d]: %d \n", topSeq, botSeq, H[indices], weight, val);
            fprintf(testfp, "chosen by = %d \n", chosenBy);
            fprintf(testfp, ".\n");
            fflush(testfp);
#endif

#if TEST == 1
        printf("Found min in part 1 \n");
        printf("[%d %d %d %d]  \n", indd[0], indd[1], indd[2], indd[3]);
        printf("With min = %f \n", min);
        printf("chosen by = %d \n", chosenBy);
        printf(".\n");
#endif

            
        }
        // restore original value so we can reuse this.
        indices[i]++;
    }

    double wt = 0, lastI = 0;

    // Second part: 
    // here we look at all possible (even, odd) pairs

    int itersForLoop = numEven * numOdd;
    for(int pairID=0; pairID<itersForLoop; pairID++)  
    {

        int odd = pairID % numOdd;
        int even = pairID / numOdd;

        int eIndex = even;
        int oIndex = numEven + odd;

        int currEvenPos = indices[eIndex];
        int currOddPos =  indices[oIndex];

        // cout << "even: " << e << "   ,   odd: " << o << endl;

        for(int w1=1; w1<25; w1++)
        {
            for(int w2=1; w2<25; w2++)
            {
                if(w1 != w2)
                    continue;

                int substrEvenPos = currEvenPos - w1 - 1;
                int substrOddPos  = currOddPos  - w2 - 1;

                if(substrEvenPos >= 0 && substrOddPos >= 0)
                {
                    double windowWeight = weights->
                            getWeight(even, odd,  currEvenPos, currOddPos, w1, w2);

                    indices[eIndex] = substrEvenPos;
                    indices[oIndex] = substrOddPos;

                    val = windowWeight + H[indices];

                    // cout << "windowWeight = " << windowWeight << endl;
                    if(val < min)
                    {
#if TEST == 1 
                        cout << "val < min = " << min << endl;
#endif
                        min = val;
                    }
                }
            }
        }

        indices[eIndex] = currEvenPos;
        indices[oIndex] = currOddPos;

    }

    // Wrap around

#if TEST == 1 
    cout << "Setting H[" << theIndex << "] = " << min << endl;
#endif

    H[theIndex] = min;
    chosenByArr[theIndex] = chosenBy;
    chosenW[theIndex] = argMinW;
    chosenW2[theIndex] = argMinW2;


#if TEST == 1
    printf("min5 = %f\n", min);
    printf("Set final H[%d %d %d] = %f or %f\n", indices[0], indices[1], indices[2], H[indices], min);
    // printf("Set final H[%d %d %d] = %f or %f\n", indices[0], indices[1], indices[2], HARRAY[theIndex], min);
    H_arrayHits++;
#endif

    for(int ci=0; ci<totalLevels; ci++)
        matrixStream << indices[ci] << "\t";
    matrixStream << min << "\n";
}




void PRBDPCore::backtrackNR(int * indices)
{   //cout << "baktrack" << endl;
    for(int ci=0;ci<data->totalLevels;ci++)
        if(indices[ci] < 0)
            return;
    
    short cb = chosenByArr[indices];

    // if chosenBy is between 1 and  n (that means it's in the first half of the algo, the one with single RNAs only - no interactions)
    if(cb >= 1 && cb <= data->totalLevels)
    {
        indices[cb-1]--;
        backtrackNR(indices);
    }
    // else if we are dealing with an interaction (if cb == 0 return)
    else if(cb > data->totalLevels) //  ( k < cb <= 2*numOfRNA - 1)
    {
        // subopt b/w 2 RNA
        short w = chosenW[indices];
        short w2 = chosenW2[indices];

        // recalculate indices. cb is >= n whereas indieces \in [0,...,n-1], there for do cb - n and cb - n -1.
        // We will put values on these indices later.
        int next1 = cb - data->totalLevels - 1;
        int next2 = cb - data->totalLevels ;
        if(next2 >= data->totalLevels)
            next2 = 0;

        int topSeq = elemToRNAmap[next1][indices[next1]];
        int botSeq = elemToRNAmap[next2][indices[next2]];

        int index1MinusOffset = indices[next1] - cutoffs[topSeq];
        int index2MinusOffset = indices[next2] - cutoffs[botSeq];

        int indexOnTopSeq = data->rnaSequences[topSeq]->compressedRNA[index1MinusOffset];
        int indexOnBotSeq = data->rnaSequences[botSeq]->compressedRNA[index2MinusOffset];

        int startOnTopSeq = indexOnTopSeq - w;
        int startOnBotSeq = indexOnBotSeq - w2;
        
        //if(gap > startOnTopSeq) 
        //startOnTopSeq -= gap;
        //if(gap > startOnBotSeq) 
        //startOnBotSeq -= gap;
        
        int uGap1 = 0, uGap2 = 0;
        if(startOnTopSeq - gap > 0) 
        {
            startOnTopSeq -= gap;
            uGap1 = 1;
        }
        
        if(startOnBotSeq - gap > 0) 
        {
            startOnBotSeq -= gap;
            uGap2 = 1;
        }
    
        int compTopSeqStart = data->rnaSequences[topSeq]->expandedRNAmap[startOnTopSeq] - 1;
        int compBotSeqStart = data->rnaSequences[botSeq]->expandedRNAmap[startOnBotSeq] - 1;

        indices[next1] = compTopSeqStart + cutoffs[topSeq];
        indices[next2] = compBotSeqStart + cutoffs[botSeq];
            
        if(uGap1)
            startOnTopSeq += gap;
            
        if(uGap2)
            startOnBotSeq += gap;
            
        printf("Interaction b/w RNA %d & %d : [%d,%d] & [%d,%d] \n", topSeq, botSeq, startOnTopSeq, indexOnTopSeq, startOnBotSeq, indexOnBotSeq);
        
        int rel = data->interLocs[topSeq][botSeq];
        double weight = (data->rnaSequences[topSeq]->type == 0) ? 
                data->rnaCollections[rel][indexOnTopSeq][indexOnBotSeq][w][w2] :
                data->rnaCollections[rel][indexOnBotSeq][indexOnTopSeq][w][w2] ;

        Window win(data->newId2OldId[topSeq], data->newId2OldId[botSeq],
                   indexOnTopSeq, indexOnBotSeq, w, w2, weight);
        
        config.add(win);

        backtrackNR(indices);
    }
}



PRBDPCore::~PRBDPCore()
{
    matrixStream.close();
    // cout << "Config: " << endl;
    // for(auto win : config)
    // {
    //  cout << win << "  " << win.weight << endl;
    // }
}
