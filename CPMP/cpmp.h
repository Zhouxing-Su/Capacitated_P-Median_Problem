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
//        loop until many no improvement search occured{
//        local_search();     // it can be local_search_with_variable_neighborhood() etc.
//        relocate_median();
//    }
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
        Solution( unsigned medianNum, unsigned vertexAllocNum, unsigned vertexNum )
        : medianList( medianNum ), medianIndex( vertexAllocNum, INVALID_INDEX ), customerIndex( vertexAllocNum ),
        assignMat( vertexAllocNum, std::vector<int>( vertexNum, INVALID_INDEX ) ),
        restCap( vertexAllocNum, 0 ), customerCount( vertexAllocNum, 0 )
        {
        }

        Solution() {}

        bool isMedian( int newMedian )
        {
            return (medianIndex[newMedian] != INVALID_INDEX);
        }

        MedianList medianList;          // has a fixed length of P
        std::vector<int> medianIndex;   // i == medianList[medianIndex[i]]
        AssignMat assignMat;            // the median itself will not appear in this matrix as a costomer
        std::vector<int> customerCount; // customerCount[i] indicate the number of customers ( other than itself ) served by vertex i if it is an median
        std::vector<CustomerIndex> customerIndex;   // i == assignMat[customerIndex[i].med][customerIndex[i].seq]
        CapacityList restCap;
    };

    struct Output
    {
        typename TopologicalGraph<T_DIST>::Distance totalDist;

        MedianList median;
        AssignList assign;

        int iterCount;
        double duration;
    };

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

    Solution curSln;

    // known conditions
    UndirectedGraph<T_DIST> graph;
    const DemandList demandList;
    const CapacityList capList;

    // solution output and statistic data
    int maxIterCount;
    Output bestSolution;
    Timer timer;
    std::string solvingAlgorithm;
};













template <typename T_DIST>
CPMP<T_DIST>::CPMP( UndirectedGraph<T_DIST> &ug, const DemandList &dl, unsigned mn, unsigned mc, int mic )
: graph( ug ), demandList( dl ), medianNum( mn ), capList( ug.vertexAllocNum, mc ),
maxIterCount( mic ), curSln( mn, ug.vertexAllocNum, ug.vertexNum )
{
    graph.getDistSeqTable();
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


}

template <typename T_DIST>
bool CPMP<T_DIST>::check() const
{
    // check if the customers are assigned to the median in the median set

    // fulfill the capacity constraint

    // recalculate the objective function from scratch
    for (AssignList::const_iterator iter = bestSolution.assign.begin(); iter != bestSolution.assign.end(); iter++) {

    }

    return true;
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
void CPMP<T_DIST>::printResult( std::ostream &os ) const
{
    os << bestSolution.totalDist << std::endl;
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
        << bestSolution.duration << ", "
        << bestSolution.iterCount << ", "
        << bestSolution.totalDist << ", ";

    for (MedianList::const_iterator iter = bestSolution.median.begin(); iter != bestSolution.median.end(); iter++) {
        csvFile << *iter << '|';
    }

    csvFile << ", ";

    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        csvFile << bestSolution.assign[i] << '|';
    }
    csvFile << std::endl;
}


#endif