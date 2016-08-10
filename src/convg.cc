#include <iostream>
#include "components/rnaseq.h"
#include "components/Weights.h"
#include "components/Cache.h"
#include "components/Chain.h"
#include "components/Config.h"
#include "components/LevelGroupProcessor.h"
#include "utils/parser.h"

#define PROG_TYPE "findone_"

void printMilliSecs(){}

void printTime() {}

int main(int argc, char *argv[])
{
	
	CmdLineArgs * args = Parser::ParseArgs(argc, argv);
	vector<OrigRNASeq> origRNASequences = Parser::GetAndParseInput();
	int numOfRNA = origRNASequences.size();

	Cache matchingCache;

	// string dir = "../rnaup_weights/output/";
	// string dir = "../output/";
	string dir = "output/";
	Weights weights(dir, &origRNASequences);
	weights.Read();

	Chain chain(&origRNASequences, numOfRNA);
	chain.printFlatStruct();
	chain.determineStruct();

	if(args->k > chain.getBBHeight())
		args->k = chain.getBBHeight();
	LevelGroupProcessor lgProc(args->k, &chain, &weights);
	
	vector<vector<int> > x = lgProc.GetAllGroupings();
	for(int groupNum=0; groupNum<x.size(); groupNum++)
	{
		Config c;
		lgProc.prepareWithSubsets(groupNum, c);
	}
	
	return 0;
}



