#include "cpmp.h"

using namespace std;

typedef int DistType;

int solve_pmedcap1( ofstream &csvFile, int instanceNum )
{

    ostringstream fname;
    fname << "../Instances/pmedcap1(" << instanceNum << ").txt";
    ifstream ifs( fname.str() );

    unsigned problemNum, optima;
    unsigned vertexNum, medianNum, medianCap;

    int nodeSeqNum;
    GeometricalGraph::Coord x, y;
    int demand;

    GeometricalGraph::PointList pl;
    CPMP<DistType>::DemandList dl;

    ifs >> problemNum >> optima;
    ifs >> vertexNum >> medianNum >> medianCap;

    while (vertexNum--) {
        ifs >> nodeSeqNum >> x >> y >> demand;
        pl.push_back( GeometricalGraph::Point( x, y ) );
        dl.push_back( demand );
    }

    ifs.close();


    GeometricalGraph gg( pl );
    //UndirectedGraph<double> dug( gg );
    //dug.getDistSeqTable();
    UndirectedGraph<DistType> uug( gg );
    uug.getDistSeqTable();

    // for each instance, run some times for judging average performance
    const int runTime = 10;
    const int maxIterCountBase = 800;
    const int tabuTenureBase = uug.vertexNum * medianNum / 4;
    const int maxNoImproveCount = tabuTenureBase * 32;
    for (int i = 1; i <= runTime; i++) {
        {
            CPMP<DistType> cpmp( uug, dl, medianNum, medianCap );
            cpmp.solve( maxIterCountBase, maxNoImproveCount, tabuTenureBase );
            cpmp.printResult( cout );
            if (!cpmp.check()) {
                csvFile << "[LogicError] ";
            }
            cpmp.appendResultToSheet( fname.str(), csvFile );
        }
    }

    return 0;
}

int main()
{
    ofstream ofs( "../Instances/log.csv", ios::app );
    //CPMP<DistType>::initResultSheet( ofs ); // call if log.csv is not exist

    solve_pmedcap1( ofs, 1 );

    ofs.close();
    return 0;
}