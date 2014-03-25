/**
*   usage :
*       (single center is not considered here)
*
*   Problems:
*       1. add seed in log file to make it possible to reproduce the calculation.
*       2.
*
====================================================================
*/
#ifndef CPMP_H
#define CPMP_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include "CPPutilibs/Graph.h"
#include "CPPutilibs/Random.h"
#include "CPPutilibs/RangeRand.h"
#include "CPPutilibs/RandSelect.h"
#include "CPPutilibs/Timer.h"

template <typename T_DIST>
class CPMP
{
public:
    const unsigned medianNum;

    typedef std::vector<int> Assignment;
    typedef std::vector<int> Demand;

    struct Solution
    {
        typename TopologicalGraph<T_DIST>::Distance totalDist;

        typename TopologicalGraph<T_DIST>::VertexSet median;
        Assignment assign;

        int iterCount;
        double duration;
    };

    CPMP( UndirectedGraph<T_DIST> &ug, Demand demand, unsigned medianNum, unsigned medianCap, int maxIterCount );
    ~CPMP();

    void solve();
    bool check() const;

    void printResult( std::ostream &os ) const;

    static void initResultSheet( std::ofstream &csvFile );
    void appendResultToSheet( const std::string &instanceFileName, std::ofstream &csvFile ) const;

private:
    // known conditions
    UndirectedGraph<T_DIST> graph;
    Demand demand;
    unsigned medianCap;

    // representation of the solution
    // typename TopologicalGraph<T_DIST>::VertexSet center;
    // Assignment assign;


    int maxIterCount;
    Solution bestSolution;
    Timer timer;
    std::string solvingAlgorithm;
};












using namespace std;

template <typename T_DIST>
CPMP<T_DIST>::CPMP( UndirectedGraph<T_DIST> &ug, Demand d, unsigned mn, unsigned mc, int mic )
: graph( ug ), demand( d ), medianNum( mn ), medianCap( mc ), maxIterCount( mic )
{
    // graph.getDistSeqTable();
}


template <typename T_DIST>
CPMP<T_DIST>::~CPMP()
{
}

template <typename T_DIST>
void CPMP<T_DIST>::solve()
{
    ostringstream ss;
    ss << '[' << typeid(T_DIST).name() << ']' << "VNS";
    solvingAlgorithm = ss.str();


}

template <typename T_DIST>
bool CPMP<T_DIST>::check() const
{
    // check if the customers are assigned to the median in the median set

    // fulfill the capacity constrain

    // recalculate the objective function from scratch
    for (Assignment::const_iterator iter = bestSolution.assign.begin( ); iter != bestSolution.assign.end( ); iter++) {
        
    }

    return true;
}

template <typename T_DIST>
void CPMP<T_DIST>::printResult( ostream &os ) const
{
    os << bestSolution.totalDist << std::endl;
}

template <typename T_DIST>
void CPMP<T_DIST>::initResultSheet( std::ofstream &csvFile )
{
    csvFile << "Date, Instance, Algorithm, TotalIter, RandSeed, Duration, IterCount, TotalDist, Centers, Assignment" << std::endl;
}

template <typename T_DIST>
void CPMP<T_DIST>::appendResultToSheet( const string &instanceFileName, ofstream &csvFile ) const
{
    csvFile << Timer::getLocalTime() << ", "
        << instanceFileName << ", "
        << solvingAlgorithm << ", "
        << maxIterCount << ", "
        << Random::seed << ','
        << bestSolution.duration << ", "
        << bestSolution.iterCount << ", "
        << bestSolution.totalDist << ", ";

    for (TopologicalGraph<T_DIST>::VertexSet::const_iterator iter = bestSolution.median.begin(); iter != bestSolution.median.end(); iter++) {
        csvFile << *iter << '|';
    }

    csvFile << ", ";

    for (Assignment::const_iterator iter = bestSolution.assign.begin(); iter != bestSolution.assign.end(); iter++) {
        csvFile << *iter << '|';
    }
    csvFile << std::endl;
}


#endif