/**
*   usage :
*       (single median is not considered here)
*
*   Problems:
*       1. add seed in log file to make it possible to reproduce the calculation.
*       2.
*
*
******** algorithm ********
//solve()
//{
//    gen_init_solution();    // it can be gen_init_solution_randomly or gen_init_solution_greedily
//
//    while (stop condition not reached) {
//        loop until many no improvement search occured {
//            local_search();     // it can be local_search_with_variable_neighborhood() etc.
//            relocate_median();
//        }
//        perturbation();     // it can be restart, crossover, path_relinking, scatter_search or other diversification procedures
//    }
//}
//
//gen_init_solution_greedily()
//{
//    select P medians randomly and set themselves as the medians of themselves.
//    for (each customer) {
//        assign it to the closest median with capacity left.
//            update the capacity.
//            add the distance to the objective function value.
//    }
//
//    while (there are nodes not assigned to median) {
//        assign it to the closest median.
//    }
//    fix_to_fit_constraint();
//}
//
//local_search_with_variable_neighborhood()
//{
//    while (true) {
//        if (search_neighborhood_of_reassign_customer_to_new_median() find an improvement) {
//            update_current_and_best_solution();
//        } else if (search_neighborhood_of_swap_customers_from_two_medians() find an improvement) {
//            update_current_and_best_solution();
//        } else {
//            local search stop;
//        }
//    }
//}
//
//search_neighborhood_of_reassign_customer_to_new_median()
//{
//    for (each customer) {
//        for (each median other than its original median) {
//            if (the capacity constraint is not violated) {
//                calculate the delta value of the objective function by assuming that the customer is reassigned to the new median.
//                if (the delta is less than the minimal delta found before) {
//                    record the index of the customer and the new median.
//                } else if (the delta is equal to the minimal delta found before) {
//                    decide whether to record the index of the customer and the new median by the probability of 1 / TotalOccurrenceOfThisMinima
//                }
//            }
//        }
//    }
//
//    return the minimal delta;
//}
******** algorithm ********
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

    typedef int Unit;    // unit of the demand and capacity
    typedef std::vector<Unit> DemandList;
    typedef std::vector<Unit> CapacityList;

    static const int INVALID_INDEX = -1;
    typedef std::vector<int> MedianList;
    typedef std::vector<int> AssignList;
    typedef std::vector< std::vector<int> > AssignMat;
    struct CustomerIndex
    {
        int med;    // the median serving this customer(the line number of the assignMat)
        int seq;    // the row number of the assignMat
    };

    struct Solution
    {
        Solution() {}
        Solution( const CPMP &cpmp ) : medianList( cpmp.medianNum ), medianIndex( cpmp.graph.vertexAllocNum, INVALID_INDEX ),
            customerIndex( cpmp.graph.vertexAllocNum ), assignMat( cpmp.graph.vertexAllocNum, std::vector<int>( cpmp.graph.vertexNum, INVALID_INDEX ) ),
            restCap( cpmp.graph.vertexAllocNum, 0 ), customerCount( cpmp.graph.vertexAllocNum, 0 ), totalDist( 0 )
        {
        }

        bool isMedian( int newMedian )
        {
            return (medianIndex[newMedian] != INVALID_INDEX);
        }

        typename TopologicalGraph<T_DIST>::Distance totalDist;

        MedianList medianList;          // has a fixed length of P
        std::vector<int> medianIndex;   // i == medianList[medianIndex[i]]
        AssignMat assignMat;            // the median itself will not appear in this matrix as a costomer
        std::vector<int> customerCount; // customerCount[i] indicate the number of customers ( other than itself ) served by vertex i if it is an median
        std::vector<CustomerIndex> customerIndex;   // i == assignMat[customerIndex[i].med][customerIndex[i].seq]
        CapacityList restCap;
    };

    struct Output
    {
        Output() {}
        Output( CPMP &cpmp, const Solution& s, int iterationCount ) : totalDist( s.totalDist ),
            iterCount( iterationCount ), median( s.medianList ), assign( cpmp.graph.vertexNum )
        {
            cpmp.timer.record();
            duration = cpmp.timer.getTotalDuration();

            int j = 0;
            for (int i = cpmp.graph.minVertexIndex; i <= cpmp.graph.maxVertexIndex; i++) {
                assign[j++] = s.customerIndex[i].med;
            }
        }

        typename TopologicalGraph<T_DIST>::Distance totalDist;

        MedianList median;
        AssignList assign;

        int iterCount;
        double duration;
    };

    // known conditions
    const UndirectedGraph<T_DIST> graph;
    const DemandList demandList;
    const CapacityList capList;

    // the map from indices to vertex of demand should be the same as the ug.
    CPMP( UndirectedGraph<T_DIST> &ug, const DemandList &demand, unsigned medianNum, unsigned medianCap, int maxIterCount );
    ~CPMP();

    void solve();
    bool check() const;

    void printResult( std::ostream &os ) const;

    static void initResultSheet( std::ofstream &csvFile );
    void appendResultToSheet( const std::string &instanceFileName, std::ofstream &csvFile ) const;

private:
    void genInitSolution();
    void genFeasibleInitSolutionByRestarting();
    void genFeasibleInitSolutionByRepairing();  // not finished

    void localSearchOnReassignCustomerNeighborhood();
    bool exploreShiftCustomerNeighborhood();
    bool exploreSwapCustomerNeighborhood();

    Solution curSln;

    // solution output and statistic data
    int maxIterCount;
    Output optima;
    Timer timer;
    std::string solvingAlgorithm;
};













template <typename T_DIST>
CPMP<T_DIST>::CPMP( UndirectedGraph<T_DIST> &ug, const DemandList &dl, unsigned mn, unsigned mc, int mic )
: graph( ug ), demandList( dl ), medianNum( mn ), capList( ug.vertexAllocNum, mc ),
maxIterCount( mic ), curSln( *this )
{
}

template <typename T_DIST>
CPMP<T_DIST>::~CPMP()
{
}

template <typename T_DIST>
void CPMP<T_DIST>::solve()
{
    std::ostringstream ss;
    ss << '[' << typeid(T_DIST).name() << ']' << "ShiftCustomerNeighborhood";
    solvingAlgorithm = ss.str();

    genInitSolution();
    optima = Output( *this, curSln, 0 );

    for (int iterCount = 0; iterCount < maxIterCount; iterCount++) {
        localSearchOnReassignCustomerNeighborhood();
        // localSearchOnRelocateMedianNeighborhood();
    }
}






template <typename T_DIST>
void CPMP<T_DIST>::localSearchOnReassignCustomerNeighborhood()
{
    bool isImproved = false;
    while (true) {
        isImproved = exploreShiftCustomerNeighborhood();
        if (!isImproved) {
            isImproved = exploreSwapCustomerNeighborhood();
            if (!isImproved) {  // local optima of the two neighborhoods found
                return;
            }
        }
    }
}

template <typename T_DIST>
bool CPMP<T_DIST>::exploreShiftCustomerNeighborhood()
{
    
}

template <typename T_DIST>
bool CPMP<T_DIST>::exploreSwapCustomerNeighborhood()
{
    return false;
}





template <typename T_DIST>
void CPMP<T_DIST>::genInitSolution()
{
    genFeasibleInitSolutionByRestarting();
    // genFeasibleInitSolutionByRepairing();    // not finished
}

template <typename T_DIST>
void CPMP<T_DIST>::genFeasibleInitSolutionByRestarting()
{
    RangeRand rr( graph.minVertexIndex, graph.maxVertexIndex );

    bool isFeasible;
    Solution sln;

    do {    // loop until a feasible init solution is generated
        isFeasible = true;
        sln = curSln;

        // select P median randomly
        for (unsigned curMedianNum = 0; curMedianNum < medianNum; curMedianNum++) {
            int newMedian;
            do {
                newMedian = rr();
            } while (sln.isMedian( newMedian ));
            // add this median to current solution
            sln.medianList[curMedianNum] = newMedian;
            sln.medianIndex[newMedian] = curMedianNum;
            sln.customerIndex[newMedian].med = newMedian;   // for convenience to generate output
            // set itself as its median and update the capacity
            sln.restCap[newMedian] = capList[newMedian] - demandList[newMedian];
        }

        // assign customers to the closest median with enough capacity
        for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
            if (!sln.isMedian( i )) {
                int med;
                int j = graph.minVertexIndex;
                for (; j <= graph.maxVertexIndex; j++) {
                    med = graph.nthClosestVertex( i, j );
                    if (sln.isMedian( med ) && (sln.restCap[med] >= demandList[i])) {
                        sln.assignMat[med][sln.customerCount[med]] = med;
                        sln.customerIndex[i].med = med;
                        sln.customerIndex[i].seq = sln.customerCount[med];
                        sln.customerCount[med]++;
                        sln.restCap[med] -= demandList[i];
                        sln.totalDist += graph.distance( med, i );
                        break;
                    }
                }
                if (j == graph.maxVertexIndex) {    // did not find a median
                    isFeasible = false;
                    break;
                }
            }
        }
    } while (!isFeasible);

    curSln = sln;
}

template <typename T_DIST>
void CPMP<T_DIST>::genFeasibleInitSolutionByRepairing()
{
    RangeRand rr( graph.minVertexIndex, graph.maxVertexIndex );

    // select P median randomly
    for (unsigned curMedianNum = 0; curMedianNum < medianNum; curMedianNum++) {
        int newMedian;
        do {
            newMedian = rr();
        } while (curSln.isMedian( newMedian ));
        // add this median to current solution
        curSln.medianList[curMedianNum] = newMedian;
        curSln.medianIndex[newMedian] = curMedianNum;
        sln.customerIndex[newMedian].med = newMedian;   // for convenience to generate output
        // set itself as its median and update the capacity
        curSln.restCap[newMedian] = capList[newMedian] - demandList[newMedian];
    }

    // assign customers to the closest median with enough capacity
    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        if (!curSln.isMedian( i )) {
            int med;
            int j = graph.minVertexIndex;
            for (; j <= graph.maxVertexIndex; j++) {
                med = graph.nthClosestVertex( i, j );
                if (curSln.isMedian( med ) && (curSln.restCap[med] >= demandList[i])) {
                    curSln.assignMat[med][curSln.customerCount[med]] = med;
                    curSln.customerIndex[i].med = med;
                    curSln.customerIndex[i].seq = curSln.customerCount[med];
                    curSln.customerCount[med]++;
                    curSln.restCap[med] -= demandList[i];
                    curSln.totalDist += graph.distance( med, i );
                    break;
                }
            }
            if (j == graph.maxVertexIndex) {    // did not find a median
                for (j = graph.minVertexIndex; j <= graph.maxVertexIndex; j++) {
                    med = graph.nthClosestVertex( i, j );
                    if (curSln.isMedian( med )) {
                        curSln.assignMat[med][curSln.customerCount[med]] = med;
                        curSln.customerIndex[i].med = med;
                        curSln.customerIndex[i].seq = curSln.customerCount[med];
                        curSln.customerCount[med]++;
                        curSln.restCap[med] -= demandList[i];
                        curSln.totalDist += graph.distance( med, i );
                        break;
                    }
                }
            }
        }
    }

    // repair current solution to a feasible solution
    //repairCurSolution();
}





template <typename T_DIST>
bool CPMP<T_DIST>::check() const
{
    const double ERROR = 0.01;
    typename TopologicalGraph<T_DIST>::Distance totalDist = 0;
    CapacityList restCap( capList );

    // recalculate the objective function from scratch
    int i = graph.minVertexIndex;
    for (AssignList::const_iterator iter = optima.assign.begin(); iter != optima.assign.end(); iter++, i++) {
        // check if the median really exist
        bool isMedian = false;
        for (MedianList::const_iterator iterm = optima.median.begin(); iterm != optima.median.end(); iterm++) {
            if (*iterm == *iter) {
                isMedian = true;
                break;
            }
        }
        // check the capacity
        if (isMedian && ((restCap[*iter] -= demandList[i]) >= 0)) {
            totalDist += graph.distance( *iter, i );
        } else {
            return false;
        }
    }

    // check the objective function value
    if ((totalDist - optima.totalDist)*(totalDist - optima.totalDist) > ERROR) {
        return false;
    }

    return true;
}






template <typename T_DIST>
void CPMP<T_DIST>::printResult( std::ostream &os ) const
{
    os << optima.totalDist << std::endl;
}

template <typename T_DIST>
void CPMP<T_DIST>::initResultSheet( std::ofstream &csvFile )
{
    csvFile << "Date, Instance, Algorithm, TotalIter, RandSeed, Duration, IterCount, TotalDist, Centers, Assignment" << std::endl;
}

template <typename T_DIST>
void CPMP<T_DIST>::appendResultToSheet( const std::string &instanceFileName, std::ofstream &csvFile ) const
{
    csvFile << Timer::getLocalTime() << ", "
        << instanceFileName << ", "
        << solvingAlgorithm << ", "
        << maxIterCount << ", "
        << Random::seed << ", "
        << optima.duration << ", "
        << optima.iterCount << ", "
        << optima.totalDist << ", ";

    for (MedianList::const_iterator iter = optima.median.begin(); iter != optima.median.end(); iter++) {
        csvFile << *iter << '|';
    }

    csvFile << ", ";

    for (AssignList::const_iterator iter = optima.assign.begin(); iter != optima.assign.end(); iter++) {
        csvFile << *iter << '|';
    }
    csvFile << std::endl;
}


#endif