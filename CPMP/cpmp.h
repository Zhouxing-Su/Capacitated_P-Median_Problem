/**
*   usage :
*       (single center is not considered here)
*
*   Problems:
*       1. add seed in log file to make it possible to reproduce the calculation.
*       2.
*
*
    ******** algorithm ********
    solve()
    {
        gen_init_solution();    // it can be gen_init_solution_randomly or gen_init_solution_greedily

        while (stop condition not reached) {
            loop until many no improvement search occured {
                local_search();     // it can be local_search_with_variable_neighborhood() etc.
                relocate_median();
            }
            perturbation();     // it can be restart, crossover, path_relinking, scatter_search or other diversification procedures
        }
    }

    gen_init_solution_greedily()
    {
        loop P times {
            select a node as a median randomly.
            while (the node has capacity left) {
                assign the closest node which has no median to this node.
                add the distance to the objective function value.
            }
        }

        while (there are nodes not assigned to median) {
            assign it to the closest median.
        }
        fix_to_fit_constraint();
    }

    local_search_with_variable_neighborhood()
    {
        while (true) {
            if (search_neighborhood_of_reassign_customer_to_new_median() find an improvement) {
                update_current_and_best_solution();
            } else if (search_neighborhood_of_swap_customers_from_two_medians() find an improvement) {
                update_current_and_best_solution();
            } else {
                local search stop;
            }
        }
    }

    search_neighborhood_of_reassign_customer_to_new_median()
    {
        for (each customer) {
            for (each median other than its original median) {
                if (the capacity constraint is not violated) {
                    calculate the delta value of the objective function by assuming that the customer is reassigned to the new median.
                    if (the delta is less than the minimal delta found before) {
                        record the index of the customer and the new median.
                    } else if (the delta is equal to the minimal delta found before) {
                        decide whether to record the index of the customer and the new median by the probability of 1/TotalOccurrenceOfThisMinima
                    }
                }
            }
        }

        return the minimal delta;
    }
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

    // representation of the current solution
    typename TopologicalGraph<T_DIST>::VertexSet center;
    Assignment assign;


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

    // fulfill the capacity constraint

    // recalculate the objective function from scratch
    for (Assignment::const_iterator iter = bestSolution.assign.begin(); iter != bestSolution.assign.end(); iter++) {

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