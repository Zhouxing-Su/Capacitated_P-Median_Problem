/**
*   Usage:
*       solve capacitated p-median problem.
*       (single median is not considered here)
*
*   Parameters:
*       DistType @ main.cpp
*       MAX_ITER_COUNT @ solve()
*       MAX_NO_IMPROVE_COUNT @ solve()
*       TABU_TENURE_ASSIGN @ solve()
*       TABU_TENURE_RELOCATE @ solve()
*       DEMAND_DISTRIBUTION_DAMPING @ selectOpenedMedian()
*
*   Problems:
*       1. add seed in log file to make it possible to reproduce the calculation.
*       2. relocateSingleMedian() can be only applied to instances with all medians having same capacity.
*       3. invoke some math function in standard library which may not compatible with T_DIST type.
*
**/

//******** algorithm ********
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
//******** algorithm ********


#ifndef CPMP_A_H
#define CPMP_A_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include "../CPPutilibs/Graph.h"
#include "../CPPutilibs/Timer.h"
#include "../CPPutilibs/Log.h"



template <typename T_DIST = int, int DIST_MULTIPLICATION = 1000>
class CPMP
{
public:
    static const int INVALID_INDEX = -1;
    static const T_DIST MAX_DIST_DELTA = 1000 * DIST_MULTIPLICATION;
    static const T_DIST MAX_TOTAL_DIST = 10000 * DIST_MULTIPLICATION;

    const unsigned medianNum;

    typedef int Unit;    // unit of the demand and capacity
    typedef std::vector<Unit> DemandList;
    typedef std::vector<Unit> CapacityList;

    typedef std::vector<int> MedianList;
    typedef std::vector<int> AssignList;

    struct Solution;
    struct Output
    {
        Output() : totalDist( MAX_TOTAL_DIST ) {}
        Output( CPMP &cpmp, const Solution& s );

        T_DIST totalDist;

        MedianList median;
        AssignList assign;

        double duration;
    };

    // known conditions
    const UndirectedGraph<T_DIST, DIST_MULTIPLICATION> graph;
    const DemandList demandList;
    const CapacityList capList;

    // the map from indices to vertex of demand should be the same as the ug.
    CPMP( UndirectedGraph<T_DIST, DIST_MULTIPLICATION> &ug, const DemandList &demand, unsigned medianNum, unsigned medianCap );
    ~CPMP();

    void solve();

    void solve_BranchAndCut();
    void locateMedian( int restMedNum, int firstAvaliableMed );
    void assignCustomer( int customer );

    bool check() const;

    void printResult( std::ostream &os ) const;

    static void initResultSheet( std::ofstream &csvFile );
    void appendResultToSheet( const std::string &instanceFileName, std::ofstream &csvFile ) const;

private:
    struct Solution
    {
        Solution() {}
        Solution( const CPMP &cpmp ) : median( cpmp.medianNum, INVALID_INDEX ), assign( cpmp.graph.vertexAllocNum, INVALID_INDEX ),
            restCap( cpmp.graph.vertexAllocNum, 0 ), totalDist( 0 )
        {
        }

        // judging by the vector assign
        bool isMedian( int vertex )
        {
            return (assign[vertex] == vertex);
        }

        T_DIST totalDist;

        MedianList median;
        AssignList assign;
        CapacityList restCap;
    };


    Solution curSln;

    Output optima;
    Timer timer;
    std::string solvingAlgorithm;
};













template <typename T_DIST, int DIST_MULTIPLICATION>
CPMP<T_DIST, DIST_MULTIPLICATION>::CPMP( UndirectedGraph<T_DIST, DIST_MULTIPLICATION> &ug, const DemandList &dl, unsigned mn, unsigned mc )
: graph( ug ), demandList( dl ), medianNum( mn ), capList( ug.vertexAllocNum, mc ), curSln( *this )
{
}

template <typename T_DIST, int DIST_MULTIPLICATION>
CPMP<T_DIST, DIST_MULTIPLICATION>::~CPMP()
{
}

template <typename T_DIST, int DIST_MULTIPLICATION>
CPMP<T_DIST, DIST_MULTIPLICATION>::Output::Output( CPMP &cpmp, const Solution& s ) :
totalDist( s.totalDist ), median( s.median ), assign( s.assign )
{
    cpmp.timer.record();
    duration = cpmp.timer.getTotalDuration();
}

template <typename T_DIST, int DIST_MULTIPLICATION>
void CPMP<T_DIST, DIST_MULTIPLICATION>::solve()
{
    solve_BranchAndCut();
}

template <typename T_DIST, int DIST_MULTIPLICATION>
void CPMP<T_DIST, DIST_MULTIPLICATION>::solve_BranchAndCut()
{
    std::ostringstream ss;
    ss << '[' << typeid(T_DIST).name() << "]BranchAndCut";
    solvingAlgorithm = ss.str();

    locateMedian( medianNum, graph.maxVertexIndex );
}

template <typename T_DIST, int DIST_MULTIPLICATION>
void CPMP<T_DIST, DIST_MULTIPLICATION>::locateMedian( int restMedNum, int lastAvaliableMed )
{
    if (restMedNum == 0) {
        assignCustomer( graph.maxVertexIndex );
    } else {
        int medCount = medianNum - restMedNum;
        restMedNum--;
        for (int med = lastAvaliableMed; med > graph.minVertexIndex; med--) {
            curSln.median[medCount] = med;
            curSln.assign[med] = med;
            curSln.restCap[med] = capList[med] - demandList[med];
            locateMedian( restMedNum, med - 1 );
        }
    }
}

template <typename T_DIST, int DIST_MULTIPLICATION>
void CPMP<T_DIST, DIST_MULTIPLICATION>::assignCustomer( int customer )
{
    if (customer < graph.minVertexIndex) { // finish assigning all customers to medians
        if (curSln.totalDist < optima.totalDist) {
            optima = Output( *this, curSln );
            Log<T_DIST>::writeln( optima.totalDist );
        }
    } else if (curSln.isMedian( customer )) {
        assignCustomer( customer - 1 );
    } else {
        for (unsigned i = 0; i < medianNum; i++) {
            int med = curSln.median[i];
            curSln.restCap[med] -= demandList[customer];
            if (curSln.restCap[med] >= 0) {
                curSln.assign[customer] = med;
                curSln.totalDist += graph.distance( med, customer );
                assignCustomer( customer - 1 );
                curSln.totalDist -= graph.distance( med, customer );
            }
            curSln.restCap[med] += demandList[customer];
        }
    }
}






template <typename T_DIST, int DIST_MULTIPLICATION>
bool CPMP<T_DIST, DIST_MULTIPLICATION>::check() const
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






template <typename T_DIST, int DIST_MULTIPLICATION>
void CPMP<T_DIST, DIST_MULTIPLICATION>::printResult( std::ostream &os ) const
{
    os << (graph.isMultiplied() ? (optima.totalDist / static_cast<double>(graph.DistMultiplication)) : optima.totalDist) << std::endl;
}

template <typename T_DIST, int DIST_MULTIPLICATION>
void CPMP<T_DIST, DIST_MULTIPLICATION>::initResultSheet( std::ofstream &csvFile )
{
    csvFile << "Date, Instance, Algorithm, MAX_ITER_COUNT, MAX_NO_IMPROVE_COUNT, RandSeed, Duration, IterCount, MoveCount, TotalDist, Medians, Assignment" << std::endl;
}

template <typename T_DIST, int DIST_MULTIPLICATION>
void CPMP<T_DIST, DIST_MULTIPLICATION>::appendResultToSheet( const std::string &instanceFileName, std::ofstream &csvFile ) const
{
    csvFile << Timer::getLocalTime() << ", "
        << instanceFileName << ", "
        << solvingAlgorithm << ", "
        << 0 << ", "
        << 0 << ", "
        << 0 << ", "
        << optima.duration << ", "
        << 0 << ", "
        << 0 << ", "
        << (graph.isMultiplied() ? (optima.totalDist / static_cast<double>(graph.DistMultiplication)) : optima.totalDist) << ", ";

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