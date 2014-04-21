/**
*   Usage:
*       solve capacitated p-median problem.
*       (single median is not considered here)
*
*   Parameters:
*       DistType @ main.cpp
*       MAX_ITER_COUNT @ solve()
*       MAX_NO_IMPROVE_COUNT @ solve()
*       TABU_TENURE_BASE @ solve()
*
*   Problems:
*       1. add seed in log file to make it possible to reproduce the calculation.
*       2. relocateSingleMedian() can be only applied to instances with all medians having same capacity.
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


#ifndef CPMP_H
#define CPMP_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include "../CPPutilibs/Graph.h"
#include "../CPPutilibs/Random.h"
#include "../CPPutilibs/RangeRand.h"
#include "../CPPutilibs/RandSelect.h"
#include "../CPPutilibs/Timer.h"
#include "../CPPutilibs/Log.h"



template <typename T_DIST = int>
class CPMP
{
public:
    static const int INVALID_INDEX = -1;
    static const T_DIST MAX_DIST_DELTA = 100000;

    const unsigned medianNum;

    typedef int Unit;    // unit of the demand and capacity
    typedef std::vector<Unit> DemandList;
    typedef std::vector<Unit> CapacityList;

    typedef std::vector<int> MedianList;
    typedef std::vector<int> AssignList;
    typedef std::vector< std::vector<int> > AssignMat;
    struct CustomerIndex    // the median itself will only record a med but the seq field is undefined
    {
        int med;    // the median serving this customer (the line number of the assignMat)
        int seq;    // the row number of the assignMat
    };

    struct Solution;
    struct Output
    {
        Output() {}
        Output( CPMP &cpmp, const Solution& s, int iterationCount, int movementCount );

        T_DIST totalDist;

        MedianList median;
        AssignList assign;

        int iterCount;
        int moveCount;
        double duration;
    };

    // known conditions
    const UndirectedGraph<T_DIST> graph;
    const DemandList demandList;
    const CapacityList capList;

    // the map from indices to vertex of demand should be the same as the ug.
    CPMP( UndirectedGraph<T_DIST> &ug, const DemandList &demand, unsigned medianNum, unsigned medianCap );
    ~CPMP();

    void solve( int maxIterCount, int maxNoImproveCount, int tabuTenureBase );
    void solve_ShiftSwapTabu( int maxIterCount, int maxNoImproveCount, int tabuTenureBase );
    void solve_ShiftSwap();
    void solve_RandomInit();
    bool check() const;

    void printResult( std::ostream &os ) const;

    static void initResultSheet( std::ofstream &csvFile );
    void appendResultToSheet( const std::string &instanceFileName, std::ofstream &csvFile ) const;

private:
    struct Move
    {
        Move() : param1( INVALID_INDEX ), param2( INVALID_INDEX ), distDelta( MAX_DIST_DELTA ) {}
        Move( int p1, int p2, T_DIST delta ) : param1( p1 ), param2( p2 ), distDelta( delta ) {}

        bool isImproved()
        {
            return (distDelta < 0);
        }

        int param1;
        int param2;
        T_DIST distDelta;
    };

    struct Solution
    {
        Solution() {}
        Solution( const CPMP &cpmp ) : medianList( cpmp.medianNum ), medianIndex( cpmp.graph.vertexAllocNum, INVALID_INDEX ),
            customerIndex( cpmp.graph.vertexAllocNum ), assignMat( cpmp.graph.vertexAllocNum, std::vector<int>( cpmp.graph.vertexNum, INVALID_INDEX ) ),
            restCap( cpmp.graph.vertexAllocNum, 0 ), customerCount( cpmp.graph.vertexAllocNum, 0 ), totalDist( 0 )
        {
        }

        // judging by the vector medianIndex
        bool isMedian( int vertex )
        {
            return (medianIndex[vertex] != INVALID_INDEX);
        }

        // move customer from oldMedian to newMedian
        // (may make the capacity a negative value)
        void shiftCustomer( const CPMP &cpmp, int customer, int newMedian, T_DIST distDelta )
        {
            int oldMedian = customerIndex[customer].med;
            // reset the index of the last customer of the old median to the index of the median to be moved 
            // and update customerCount of the old median
            customerIndex[assignMat[oldMedian][--customerCount[oldMedian]]].seq = customerIndex[customer].seq;
            // replace the customer to be moved with the last customer of the old median
            assignMat[oldMedian][customerIndex[customer].seq] = assignMat[oldMedian][customerCount[oldMedian]];
            // assign the customer to the new median and update customerCount of the new median
            customerIndex[customer].med = newMedian;
            customerIndex[customer].seq = customerCount[newMedian];
            assignMat[newMedian][customerCount[newMedian]++] = customer;
            // update the capacity of the medians
            restCap[oldMedian] += cpmp.demandList[customer];
            restCap[newMedian] -= cpmp.demandList[customer];
            // update the total distance
            totalDist += distDelta;
        }

        // exchange the medians of cutsomer1 and customer2
        // (may make the capacity a negative value)
        void swapCustomer( const CPMP &cpmp, int customer1, int customer2, T_DIST distDelta )
        {
            int median1 = customerIndex[customer1].med;
            int median2 = customerIndex[customer2].med;
            int seq1 = customerIndex[customer1].seq;
            int seq2 = customerIndex[customer2].seq;
            // update assign matrix and its indices
            customerIndex[customer1].med = median2;
            customerIndex[customer2].med = median1;
            customerIndex[customer1].seq = seq2;
            customerIndex[customer2].seq = seq1;
            assignMat[median1][seq1] = customer2;
            assignMat[median2][seq2] = customer1;
            // update the capacity of the medians
            Unit capDelta = cpmp.demandList[customer1] - cpmp.demandList[customer2];
            restCap[median1] += capDelta;
            restCap[median2] -= capDelta;
            // update the total distance
            totalDist += distDelta;
        }

        T_DIST totalDist;

        MedianList medianList;          // has a fixed length of P
        std::vector<int> medianIndex;   // i == medianList[medianIndex[i]]
        // line i record customers assigned to median i (if it is a median).
        AssignMat assignMat;    // the median itself will not appear in this matrix as a costomer
        std::vector<int> customerCount; // customerCount[i] indicate the number of customers ( other than itself ) served by vertex i if it is an median
        std::vector<CustomerIndex> customerIndex;   // i == assignMat[customerIndex[i].med][customerIndex[i].seq]
        CapacityList restCap;
    };

    // the lines indicate customers and columns indicate medians
    typedef std::vector< std::vector<int> > ReassignTabu;


    void genInitSolution();
    void genFeasibleInitSolutionByRestarting();
    void genFeasibleInitSolutionByRepairing();  // not finished

    // end a local search when there is no valid move or local optima found
    void localSearchOnReassignCustomerNeighborhood();
    Move exploreShiftCustomerNeighborhood();
    Move exploreSwapCustomerNeighborhood();

    // end a tabu search when there is no valid move or a certain number of no improvement move occur
    void tabuSearchOnReassignCustomerNeighborhood( int noImproveCountDown );
    Move exploreTabuShiftCustomerNeighborhood();
    Move exploreTabuSwapCustomerNeighborhood();

    void relocateSingleMedian();
    void relocateSingleMedianWithTotallyReassign();
    void relocateSingleMedianWithMinimalReassign();
    int selectClosedMedian();   // return the median to be closed
    int selectOpenedMedian();   // return the median to be opened

    void perturbCustomerAssignment();

    void recoverByReassign();

    Solution curSln;
    Solution optimaOnCurrentMedianDistribution;
    bool medianDistributionChanged;
    // a N*N matrix where N is the vertex number.
    // reassignTabu[c][m] record if the customer c can be assigned to median m
    ReassignTabu reassignTabu;

    // solution output and statistic data
    int validMoveCount;
    int invalidMoveCount;
    int optimaReachCount;

    int MAX_ITER_COUNT;
    int MAX_NO_IMPROVE_COUNT;
    int TABU_TENURE_BASE;
    int iterCount;  // iteration count of local search or tabu search
    int moveCount;  // total move count of local search or tabu search

    Output optima;
    Timer timer;
    std::string solvingAlgorithm;
};













template <typename T_DIST>
CPMP<T_DIST>::CPMP( UndirectedGraph<T_DIST> &ug, const DemandList &dl, unsigned mn, unsigned mc )
: graph( ug ), demandList( dl ), medianNum( mn ), capList( ug.vertexAllocNum, mc ), curSln( *this ),
reassignTabu( ug.vertexAllocNum, std::vector<int>( ug.vertexAllocNum, 0 ) ), iterCount( 0 ), moveCount( 0 ),
MAX_ITER_COUNT( 1 ), MAX_NO_IMPROVE_COUNT( 1 ), validMoveCount( 0 ), invalidMoveCount( 0 ), medianDistributionChanged( false )
{
    Random::setSeed();
}

template <typename T_DIST>
CPMP<T_DIST>::~CPMP()
{
}

template <typename T_DIST>
CPMP<T_DIST>::Output::Output( CPMP &cpmp, const Solution& s, int iterationCount, int movementCount ) : totalDist( s.totalDist ),
iterCount( iterationCount ), median( s.medianList ), assign( cpmp.graph.vertexNum ), moveCount( movementCount )
{
    cpmp.timer.record();
    duration = cpmp.timer.getTotalDuration();

    int j = 0;
    for (int i = cpmp.graph.minVertexIndex; i <= cpmp.graph.maxVertexIndex; i++) {
        assign[j++] = s.customerIndex[i].med;
    }
}

template <typename T_DIST>
void CPMP<T_DIST>::solve( int maxIterCount, int maxNoImproveCount, int tabuTenureBase )
{
    MAX_ITER_COUNT = maxIterCount;
    MAX_NO_IMPROVE_COUNT = maxNoImproveCount;
    TABU_TENURE_BASE = tabuTenureBase;

    std::ostringstream ss;
    ss << '[' << typeid(T_DIST).name() << ']' << "ShiftSwapTabuRelocate(B=" << TABU_TENURE_BASE << ')';
    solvingAlgorithm = ss.str();

    genInitSolution();

    RangeRand rr( 1, 4 );

    while (iterCount++ < MAX_ITER_COUNT) {
        tabuSearchOnReassignCustomerNeighborhood( MAX_NO_IMPROVE_COUNT );
        if (rr() == 1) {    // relocate a median() with a probability of 1/4
            relocateSingleMedian();
        } else {    // perturb customer assignment with a probability of 3/4
            // perturbCustomerAssignment();
        }
    }
    //std::cout << "[invalidMoveCount] " << invalidMoveCount << " [validMoveCount] " << validMoveCount << std::endl;
    Log<>::write( "[optimaReachCount] " );
    Log<int>::writeln( optimaReachCount );
}

template <typename T_DIST>
void CPMP<T_DIST>::solve_ShiftSwapTabu( int maxIterCount, int maxNoImproveCount, int tabuTenureBase )
{
    MAX_ITER_COUNT = maxIterCount;
    MAX_NO_IMPROVE_COUNT = maxNoImproveCount;
    TABU_TENURE_BASE = tabuTenureBase;

    std::ostringstream ss;
    ss << '[' << typeid(T_DIST).name() << ']' << "ShiftSwapTabu(B=" << TABU_TENURE_BASE << ')';
    solvingAlgorithm = ss.str();

    genInitSolution();

    while (iterCount++ < MAX_ITER_COUNT) {
        tabuSearchOnReassignCustomerNeighborhood( MAX_NO_IMPROVE_COUNT );
    }
}

template <typename T_DIST>
void CPMP<T_DIST>::solve_ShiftSwap()
{
    std::ostringstream ss;
    ss << '[' << typeid(T_DIST).name() << ']' << "ShiftSwap";
    solvingAlgorithm = ss.str();

    genInitSolution();

    iterCount++;
    localSearchOnReassignCustomerNeighborhood();
}

template <typename T_DIST>
void CPMP<T_DIST>::solve_RandomInit()
{
    std::ostringstream ss;
    ss << '[' << typeid(T_DIST).name() << ']' << "RandomInit";
    solvingAlgorithm = ss.str();

    genInitSolution();
}




template <typename T_DIST>
void CPMP<T_DIST>::tabuSearchOnReassignCustomerNeighborhood( int noImproveCountDown )
{
    while (noImproveCountDown) {
        Move shift, swap;
        shift = exploreTabuShiftCustomerNeighborhood();
        if (!shift.isImproved()) {
            swap = exploreTabuSwapCustomerNeighborhood();
            if (!swap.isImproved()) {  // local optima of the two neighborhoods found
                noImproveCountDown--;
                // record the local optima if it is better
                if (curSln.totalDist <= optimaOnCurrentMedianDistribution.totalDist || medianDistributionChanged) {
                    optimaOnCurrentMedianDistribution = curSln;
                    medianDistributionChanged = false;
                    if (curSln.totalDist < optima.totalDist) {
                        optimaReachCount = 1;
                        optima = Output( *this, curSln, iterCount, moveCount );
                    } else if (curSln.totalDist == optima.totalDist) {
                        optimaReachCount++;
                    }
                }
            }
        }

        // execute the move by the better neighborhood
        if (shift.distDelta < swap.distDelta) {
            reassignTabu[shift.param1][shift.param2] = moveCount + TABU_TENURE_BASE;
            curSln.shiftCustomer( *this, shift.param1, shift.param2, shift.distDelta );
        } else if (swap.distDelta == MAX_DIST_DELTA) {  // shift >= swap == MAX_DIST_DELTA
            return;   // no valid move
        } else {
            reassignTabu[swap.param1][curSln.customerIndex[swap.param1].med] = moveCount + TABU_TENURE_BASE;
            reassignTabu[swap.param2][curSln.customerIndex[swap.param2].med] = moveCount + TABU_TENURE_BASE;
            curSln.swapCustomer( *this, swap.param1, swap.param2, swap.distDelta );
        }

        moveCount++;
    }
}

template <typename T_DIST>
typename CPMP<T_DIST>::Move CPMP<T_DIST>::exploreTabuShiftCustomerNeighborhood()
{
    T_DIST delta;
    T_DIST minDelta( MAX_DIST_DELTA );
    int customer = INVALID_INDEX;
    int newMedian = INVALID_INDEX;

    RandSelect rs;

    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        if (!curSln.isMedian( i )) {
            T_DIST originalDist = graph.distance( i, curSln.customerIndex[i].med );
            for (MedianList::const_iterator iter = curSln.medianList.begin();
                iter != curSln.medianList.end(); iter++) {
                if (curSln.customerIndex[i].med != *iter) {
                    if (curSln.restCap[*iter] >= demandList[i]) {
                        delta = graph.distance( i, *iter ) - originalDist;
                        if ((reassignTabu[i][*iter] < moveCount)
                            || ((curSln.totalDist + delta) < optima.totalDist)) {
                            if ((delta == minDelta) && rs.isSelected()) {
                                minDelta = delta;
                                customer = i;
                                newMedian = *iter;
                            } else if ((delta < minDelta)) {
                                rs.reset( 2 );
                                minDelta = delta;
                                customer = i;
                                newMedian = *iter;
                            }
                        }
                    }
                }
            }
        }
    }

    return Move( customer, newMedian, minDelta );
}

template <typename T_DIST>
typename CPMP<T_DIST>::Move CPMP<T_DIST>::exploreTabuSwapCustomerNeighborhood()
{
    T_DIST minDelta( MAX_DIST_DELTA );
    int c1 = INVALID_INDEX;
    int c2 = INVALID_INDEX;

    RandSelect rs;

    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        int medi = curSln.customerIndex[i].med;
        for (int j = graph.minVertexIndex; j <= graph.maxVertexIndex; j++) {
            int medj = curSln.customerIndex[j].med;
            if ((!curSln.isMedian( i )) && !curSln.isMedian( j ) && (medi != medj)) {
                Unit capDelta = demandList[i] - demandList[j];
                if (((curSln.restCap[medi] + capDelta) > 0)
                    && ((curSln.restCap[medj] - capDelta) > 0)) {
                    T_DIST distDelta = graph.distance( i, medj )
                        + graph.distance( j, medi )
                        - graph.distance( i, medi )
                        - graph.distance( j, medj );
                    if ((reassignTabu[i][medj] < moveCount)
                        || (reassignTabu[j][medi] < moveCount)
                        || ((curSln.totalDist + distDelta) < optima.totalDist)) {
                        if ((distDelta == minDelta) && rs.isSelected()) {
                            minDelta = distDelta;
                            c1 = i;
                            c2 = j;
                        } else if ((distDelta < minDelta)) {
                            rs.reset( 2 );
                            minDelta = distDelta;
                            c1 = i;
                            c2 = j;
                        }
                    }
                }
            }
        }
    }

    return Move( c1, c2, minDelta );
}

template <typename T_DIST>
void CPMP<T_DIST>::localSearchOnReassignCustomerNeighborhood()
{
    while (true) {
        Move shift, swap;
        shift = exploreShiftCustomerNeighborhood();
        if (!shift.isImproved()) {
            swap = exploreSwapCustomerNeighborhood();
            if (!swap.isImproved()) {  // local optima of the two neighborhoods found
                // record the local optima if it is better
                if (curSln.totalDist < optima.totalDist) {
                    optima = Output( *this, curSln, iterCount, moveCount );
                }
                return;
            }
        }

        // execute the move by the better neighborhood
        if (shift.distDelta < swap.distDelta) {
            curSln.shiftCustomer( *this, shift.param1, shift.param2, shift.distDelta );
        } else if (swap.distDelta == MAX_DIST_DELTA) {  // shift >= swap == MAX_DIST_DELTA
            return;   // no valid move
        } else {
            curSln.swapCustomer( *this, swap.param1, swap.param2, swap.distDelta );
        }
    }
}

template <typename T_DIST>
typename CPMP<T_DIST>::Move CPMP<T_DIST>::exploreShiftCustomerNeighborhood()
{
    T_DIST delta;
    T_DIST minDelta( MAX_DIST_DELTA );
    int customer = INVALID_INDEX;
    int newMedian = INVALID_INDEX;

    RandSelect rs;

    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        T_DIST originalDist = graph.distance( i, curSln.customerIndex[i].med );
        for (MedianList::const_iterator iter = curSln.medianList.begin();
            iter != curSln.medianList.end(); iter++) {
            if ((!curSln.isMedian( i )) && (curSln.customerIndex[i].med != *iter)) {
                if (curSln.restCap[*iter] >= demandList[i]) {
                    delta = graph.distance( i, *iter ) - originalDist;
                    if ((delta == minDelta) && rs.isSelected()) {
                        minDelta = delta;
                        customer = i;
                        newMedian = *iter;
                    } else if ((delta < minDelta)) {
                        rs.reset( 2 );
                        minDelta = delta;
                        customer = i;
                        newMedian = *iter;
                    }
                }
            }
        }
    }

    return Move( customer, newMedian, minDelta );
}

template <typename T_DIST>
typename CPMP<T_DIST>::Move CPMP<T_DIST>::exploreSwapCustomerNeighborhood()
{
    T_DIST minDelta( MAX_DIST_DELTA );
    int c1 = INVALID_INDEX;
    int c2 = INVALID_INDEX;

    RandSelect rs;

    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        int medi = curSln.customerIndex[i].med;
        for (int j = graph.minVertexIndex; j <= graph.maxVertexIndex; j++) {
            int medj = curSln.customerIndex[j].med;
            if ((!curSln.isMedian( i )) && !curSln.isMedian( j ) && (medi != medj)) {
                Unit capDelta = demandList[i] - demandList[j];
                if (((curSln.restCap[medi] + capDelta) > 0)
                    && ((curSln.restCap[medj] - capDelta) > 0)) {
                    T_DIST distDelta = graph.distance( i, medj )
                        + graph.distance( j, medi )
                        - graph.distance( i, medi )
                        - graph.distance( j, medj );
                    if ((distDelta == minDelta) && rs.isSelected()) {
                        minDelta = distDelta;
                        c1 = i;
                        c2 = j;
                    } else if ((distDelta < minDelta)) {
                        rs.reset( 2 );
                        minDelta = distDelta;
                        c1 = i;
                        c2 = j;
                    }
                }
            }
        }
    }

    return Move( c1, c2, minDelta );
}


template <typename T_DIST>
void CPMP<T_DIST>::relocateSingleMedian()
{
    relocateSingleMedianWithTotallyReassign();
    //relocateSingleMedianWithMinimalReassign();

    medianDistributionChanged = true;
}

template <typename T_DIST>
void CPMP<T_DIST>::relocateSingleMedianWithTotallyReassign()
{
    bool isFeasible = true;
    int loopCount = 1;

    int closedMedian = selectClosedMedian();
    int closedMedianIndex = optimaOnCurrentMedianDistribution.medianIndex[closedMedian];

    RangeRand orr( graph.minVertexIndex, graph.maxVertexIndex );

    while (true) {    // loop until a feasible solution is generated
        // select a median to be opened
        int openedMedian;
        do {
            openedMedian = orr();
        } while (curSln.isMedian( openedMedian ));

        // replace the closed median with the opened median in median list
        curSln.medianList[closedMedianIndex] = openedMedian;
        curSln.medianIndex[closedMedian] = INVALID_INDEX;
        curSln.medianIndex[openedMedian] = closedMedianIndex;
        curSln.customerIndex[openedMedian].med = openedMedian;

        // reproduce the assignment ( refer to genFeasibleInitSolutionByRestarting() )
        curSln.totalDist = 0;

        RangeRand rr( graph.minVertexIndex, graph.maxVertexIndex );

        for (unsigned i = 0; i < medianNum; i++) {
            int med = curSln.medianList[i];
            curSln.customerCount[med] = 0;
            curSln.restCap[med] = capList[med] - demandList[med];
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
                    isFeasible = false;
                    break;
                }
            }
        }

        if (isFeasible) {
            break;
        } else {    // recover the median distribution
            curSln.medianList[closedMedianIndex] = closedMedian;
            curSln.medianIndex[closedMedian] = closedMedianIndex;
            curSln.medianIndex[openedMedian] = INVALID_INDEX;
            isFeasible = true;
        }
    }
}

template <typename T_DIST>
void CPMP<T_DIST>::relocateSingleMedianWithMinimalReassign()
{
    // select a median to be closed
    RangeRand crr( 0, medianNum - 1 );
    int closedMedianIndex = crr();
    int closedMedian = curSln.medianList[closedMedianIndex];

    // select a median to be opened
    RangeRand orr( graph.minVertexIndex, graph.maxVertexIndex );
    int openedMedian;
    do {
        openedMedian = orr();
    } while (curSln.isMedian( openedMedian ));

    // remove the openedMedian from its old median
    int oldMedian = curSln.customerIndex[openedMedian].med;
    curSln.assignMat[oldMedian][curSln.customerIndex[openedMedian].seq] =
        curSln.assignMat[oldMedian][--curSln.customerCount[oldMedian]];
    curSln.customerIndex[openedMedian].med = openedMedian;
    // replace the closed median with the opened median in median list
    curSln.medianList[closedMedianIndex] = openedMedian;
    curSln.medianIndex[closedMedian] = INVALID_INDEX;
    curSln.medianIndex[openedMedian] = closedMedianIndex;
    // reassign all customers of the closed median(include itself) to the opened median
    curSln.totalDist = 0;
    curSln.assignMat[openedMedian][0] = closedMedian;
    curSln.customerIndex[closedMedian].med = openedMedian;
    curSln.customerIndex[closedMedian].seq = 0;
    for (int i = curSln.customerCount[closedMedian]; i > 0; i--) {
        int customer = curSln.assignMat[closedMedian][i - 1];
        curSln.assignMat[openedMedian][i] = customer;
        curSln.customerIndex[customer].med = openedMedian;
        curSln.customerIndex[customer].seq = i;
        curSln.totalDist -= graph.distance( closedMedian, customer );
        curSln.totalDist += graph.distance( openedMedian, customer );
    }
    curSln.customerCount[openedMedian] = curSln.customerCount[closedMedian] + 1;
    curSln.customerCount[closedMedian] = 0;
    // update rest capacity and deal with capacity overflow
    curSln.restCap[openedMedian] = curSln.restCap[closedMedian] - demandList[openedMedian];
    curSln.restCap[closedMedian] = 0;
    recoverByReassign();
}

template <typename T_DIST>
int CPMP<T_DIST>::selectClosedMedian()
{
    int closedMedian = optimaOnCurrentMedianDistribution.medianList[0];

    int med = optimaOnCurrentMedianDistribution.medianList[0];
    T_DIST v = 1;
    int customerNum = optimaOnCurrentMedianDistribution.customerCount[med];
    for (int j = 0; j < customerNum; j++) {
        v += graph.distance( med, optimaOnCurrentMedianDistribution.assignMat[med][j] );
    }
    v *= (optimaOnCurrentMedianDistribution.restCap[med] + 1);
    v /= (customerNum + 1); // in case the median only serve for itself

    T_DIST max = v;
    RandSelect rs;

    for (unsigned i = 1; i < medianNum; i++) {
        med = optimaOnCurrentMedianDistribution.medianList[i];
        v = 1;
        customerNum = optimaOnCurrentMedianDistribution.customerCount[med];
        for (int j = 0; j < customerNum; j++) {
            v += graph.distance( med, optimaOnCurrentMedianDistribution.assignMat[med][j] );
        }
        // because the capacity is the same for all medians,
        // use sln.restCap[med] instead of the ratio to reduce calculation.
        // otherwise it should be "sln.restCap[med] / capList[med];"
        v *= (optimaOnCurrentMedianDistribution.restCap[med] + 1);
        v /= (customerNum + 1); // in case the median only serve for itself

        if ((v == max) && rs.isSelected()) {
            closedMedian = med;
            max = v;
        } else if (v > max) {
            closedMedian = med;
            max = v;
            rs.reset( 2 );
        }
    }

    return closedMedian;
}

template <typename T_DIST>
int CPMP<T_DIST>::selectOpenedMedian()
{

}

template <typename T_DIST>
void CPMP<T_DIST>::perturbCustomerAssignment()
{

}

template <typename T_DIST>
void CPMP<T_DIST>::recoverByReassign()
{

}



template <typename T_DIST>
void CPMP<T_DIST>::genInitSolution()
{
    genFeasibleInitSolutionByRestarting();
    // genFeasibleInitSolutionByRepairing();    // not finished

    optimaOnCurrentMedianDistribution = curSln;
    optima = Output( *this, curSln, 0, 0 );
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
    os << (graph.isMultiplied() ? (optima.totalDist / static_cast<double>(graph.DistMultiplication)) : optima.totalDist) << std::endl;
}

template <typename T_DIST>
void CPMP<T_DIST>::initResultSheet( std::ofstream &csvFile )
{
    csvFile << "Date, Instance, Algorithm, MAX_ITER_COUNT, MAX_NO_IMPROVE_COUNT, RandSeed, Duration, IterCount, MoveCount, TotalDist, Medians, Assignment" << std::endl;
}

template <typename T_DIST>
void CPMP<T_DIST>::appendResultToSheet( const std::string &instanceFileName, std::ofstream &csvFile ) const
{
    csvFile << Timer::getLocalTime() << ", "
        << instanceFileName << ", "
        << solvingAlgorithm << ", "
        << MAX_ITER_COUNT << ", "
        << MAX_NO_IMPROVE_COUNT << ", "
        << Random::getSeed() << ", "
        << optima.duration << ", "
        << optima.iterCount << ", "
        << optima.moveCount << ", "
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