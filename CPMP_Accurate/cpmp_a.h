/**
*   Usage:
*       solve capacitated p-median problem.
*       (single median is not considered here)
*
*   Parameters:
*       DistType @ main.cpp
*
*   Problems:
*       1./add seed in log file to make it possible to reproduce the calculation.
*       2. relocateSingleMedian() can be only applied to instances with all medians having same capacity.
*       3. invoke some math function in standard library which may not compatible with T_DIST type.
*
**/


#ifndef CPMP_A_H
#define CPMP_A_H

#include "boost/thread/thread.hpp"

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

    const int medianNum;

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

    bool check() const;

    void printOptima( std::ostream &os ) const;

    static void initResultSheet( std::ofstream &csvFile );
    void appendResultToSheet( const std::string &instanceFileName, std::ofstream &csvFile ) const;

private:
    struct Solution
    {
        Solution() {}
        Solution( const CPMP *p ) : cpmp( p ), median( p->medianNum, INVALID_INDEX ),
            assign( p->graph.vertexAllocNum, INVALID_INDEX ),
            restCap( p->graph.vertexAllocNum, 0 ), totalDist( 0 ),
            closestMedian( p->graph.vertexAllocNum, MedianList( p->medianNum ) )
        {
        }

        // judging by the vector assign
        bool isMedian( int vertex ) const
        {
            return (assign[vertex] == vertex);
        }

        T_DIST estimateLowerBound() const;
        void genClosestMedianList();    // sort medians by distance for each customer

        const CPMP *cpmp;
        std::vector<MedianList> closestMedian;  // sorted medians by distance for each customer

        T_DIST totalDist;

        MedianList median;
        AssignList assign;
        CapacityList restCap;
    };

    void init();
    void initCustomerIndexByDescendDemand();

    void estimateUpperBound();
    void locateMedian( int restMedNum, int firstAvaliableMed );
    void assignCustomer( int index );


    std::vector<int> customerIndexByDescendDemand;
    Solution curSln;

    Output optima;
    Timer timer;
    std::string solvingAlgorithm;

    // for statistic
    friend class Daemon;
    class Daemon
    {
    public:
        // bind observed object
        Daemon( CPMP &p ) :cpmp( p ) {}

        void operator()()
        {
            while (true) {
                std::cin.ignore();  // wait for ENTER key press
                displayAll();
            }
        }

        void displayAll()
        {
            displayCutState();
            displayMedianDistribution();
        }

        void displayCutState( std::ostream &os = std::cout )
        {
            cpmp.cut[1][CUT_FAIL] = cpmp.cut[2][CUT_SUCCESS] + cpmp.cut[2][CUT_FAIL];

            for (int i = 0; i < cpmp.CUT_POINT; i++) {
                os << ' ';
                for (int j = 0; j < cpmp.CUT_STATE; j++) {
                    os << cpmp.cut[i][j] << ' ';
                }
                os << std::endl;
            }
        }

        void displayMedianDistribution( std::ostream &os = std::cout )
        {
            os << ' ';
            for (MedianList::const_iterator iter = cpmp.curSln.median.begin(); iter != cpmp.curSln.median.end(); iter++) {
                os << *iter << '|';
            }
            os << std::endl;
        }

        void displayMedianDistribution( T_DIST estimateLowerBount, std::ostream &os = std::cout )
        {
            os << " [" << estimateLowerBount << "] ";
            for (MedianList::const_iterator iter = cpmp.curSln.median.begin(); iter != cpmp.curSln.median.end(); iter++) {
                os << *iter << '|';
            }
            os << std::endl;
        }

    private:
        CPMP &cpmp;
    };

    Daemon daemon;
    static const int CUT_POINT = 3;
    static const int CUT_STATE = 2;
    static const int CUT_FAIL = 0;
    static const int CUT_SUCCESS = 1;
    std::vector<std::vector<long long> > cut;
    boost::thread tDaemon;
};










template <typename T_DIST, int DIST_MULTIPLICATION>
CPMP<T_DIST, DIST_MULTIPLICATION>::CPMP( UndirectedGraph<T_DIST, DIST_MULTIPLICATION> &ug, const DemandList &dl, unsigned mn, unsigned mc )
    : graph( ug ), demandList( dl ), medianNum( mn ), capList( ug.vertexAllocNum, mc ), curSln( this ),
    daemon( *this ), cut( CUT_POINT, std::vector<long long>( CUT_STATE, 0 ) ),
    customerIndexByDescendDemand( ug.vertexAllocNum )
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
void CPMP<T_DIST, DIST_MULTIPLICATION>::init()
{
    initCustomerIndexByDescendDemand();
    estimateUpperBound();
    tDaemon = boost::thread( daemon );
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

    init();

    locateMedian( medianNum, graph.maxVertexIndex );
}

template <typename T_DIST, int DIST_MULTIPLICATION>
void CPMP<T_DIST, DIST_MULTIPLICATION>::estimateUpperBound()
{
    
}

template <typename T_DIST, int DIST_MULTIPLICATION>
void CPMP<T_DIST, DIST_MULTIPLICATION>::locateMedian( int restMedNum, int lastAvaliableMed )
{
    if (restMedNum == 0) {
        curSln.genClosestMedianList();
        // [cut]0 lower bound of median distribution greater than optima
        if (curSln.estimateLowerBound() <= optima.totalDist) {
            cut[0][CUT_FAIL]++;
            //daemon.displayAll();
            curSln.totalDist = 0;
            assignCustomer( graph.maxVertexIndex );
        } else {
            cut[0][CUT_SUCCESS]++;
        }
    } else {
        int medCount = medianNum - restMedNum;
        restMedNum--;
        for (int med = lastAvaliableMed; med >= graph.minVertexIndex; med--) {
            curSln.median[medCount] = med;
            curSln.assign[med] = med;
            curSln.restCap[med] = capList[med] - demandList[med];
            locateMedian( restMedNum, med - 1 );
            curSln.assign[med] = INVALID_INDEX;
        }
    }
}

template <typename T_DIST, int DIST_MULTIPLICATION>
void CPMP<T_DIST, DIST_MULTIPLICATION>::assignCustomer( int index )
{
    static const int TRY_MEDIAN_NUM = medianNum / 2;

    if (index < graph.minVertexIndex) { // finish assigning all customers to medians
        if (curSln.totalDist < optima.totalDist) {
            optima = Output( *this, curSln );
            Log<T_DIST>::writeln( optima.totalDist );
        }
    } else {
        int customer = customerIndexByDescendDemand[index];
        if (curSln.isMedian( customer )) {
            assignCustomer( index - 1 );
        } else {
            for (int i = 0; i < TRY_MEDIAN_NUM; i++) {
                int med = curSln.closestMedian[customer][i];
                curSln.restCap[med] -= demandList[customer];
                // [cut]1 capacity overflow
                if (curSln.restCap[med] >= 0) {
                    curSln.assign[customer] = med;
                    curSln.totalDist += graph.distance( med, customer );
                    // [cut]2 partial distance sum greater than optima
                    if (curSln.totalDist <= optima.totalDist) {
                        cut[2][CUT_FAIL]++;
                        assignCustomer( index - 1 );
                    } else {
                        cut[2][CUT_SUCCESS]++;
                    }
                    curSln.totalDist -= graph.distance( med, customer );
                } else {
                    cut[1][CUT_SUCCESS]++;
                }
                curSln.restCap[med] += demandList[customer];
            }
        }
    }
}

template <typename T_DIST, int DIST_MULTIPLICATION>
T_DIST CPMP<T_DIST, DIST_MULTIPLICATION>::Solution::estimateLowerBound() const
{
    T_DIST dist = 0;
    CapacityList rc = cpmp->capList;

    std::vector<T_DIST> minDistDelta( cpmp->graph.vertexAllocNum, cpmp->MAX_TOTAL_DIST );

    // get the distance sum with relaxation
    for (int i = cpmp->graph.minVertexIndex; i <= cpmp->graph.maxVertexIndex; i++) {
        if (!isMedian( i )) {
            int med = closestMedian[i][0];
            dist += cpmp->graph.distance( i, med );
            rc[med] -= cpmp->demandList[i];
        }
    }

    // record the minimal distance delta of every over capacitatied medians
    for (int i = cpmp->graph.minVertexIndex; i <= cpmp->graph.maxVertexIndex; i++) {
        int med = closestMedian[i][0];
        if (restCap[med] < 0) {
            T_DIST distDelta = (cpmp->graph.distance( i, closestMedian[i][1] )
                - cpmp->graph.distance( i, med ));
            if (minDistDelta[med] > distDelta) {
                minDistDelta[med] = distDelta;
            }
        }
    }
    // compensate(increase) the lower bound by making minimal distance delta shifting
    for (int i = cpmp->graph.minVertexIndex; i <= cpmp->graph.maxVertexIndex; i++) {
        if (minDistDelta[i] != cpmp->MAX_TOTAL_DIST) {
            dist += minDistDelta[i];
        }
    }

    return dist;
}

template <typename T_DIST, int DIST_MULTIPLICATION>
void CPMP<T_DIST, DIST_MULTIPLICATION>::Solution::genClosestMedianList()
{
    for (int c = cpmp->graph.minVertexIndex; c <= cpmp->graph.maxVertexIndex; c++) {
        closestMedian[c][0] = median[0];
        for (int i = 1; i < cpmp->medianNum; i++) {
            int j = i;
            T_DIST dist = cpmp->graph.distance( c, median[i] );
            for (; j > 0; j--) {
                int k = j - 1;
                if (cpmp->graph.distance( c, median[k] ) > dist) {
                    closestMedian[c][j] = closestMedian[c][k];
                } else {
                    break;
                }
            }
            closestMedian[c][j] = median[i];
        }
    }
}

template <typename T_DIST, int DIST_MULTIPLICATION>
void CPMP<T_DIST, DIST_MULTIPLICATION>::initCustomerIndexByDescendDemand()
{
    // insertion sort
    customerIndexByDescendDemand[graph.minVertexIndex] = graph.minVertexIndex;
    for (int i = graph.minVertexIndex + 1; i <= graph.maxVertexIndex; i++) {
        int j = i;
        for (; j > graph.minVertexIndex; j--) {
            int k = j - 1;
            if (demandList[i] > demandList[customerIndexByDescendDemand[k]]) {
                customerIndexByDescendDemand[j] = customerIndexByDescendDemand[k];
            } else {
                break;
            }
        }
        customerIndexByDescendDemand[j] = i;
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
void CPMP<T_DIST, DIST_MULTIPLICATION>::printOptima( std::ostream &os ) const
{
    os << "Opt. "
        << (graph.isMultiplied() ? (optima.totalDist / static_cast<double>(graph.DistMultiplication)) : optima.totalDist)
        << std::endl;
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