#include "cpmp.h"

using namespace std;


int solve_pmedcap1( ofstream &csvFile, int instanceNum )
{
    ostringstream fname;
    fname << "pmedcap1(" << instanceNum << ").txt";
    ifstream ifs( fname.str() );

    unsigned problemNum, optima;
    unsigned vertexNum, medianNum, medianCap;

    int nodeSeqNum;
    GeometricalGraph::Coord x, y;
    int demand;

    GeometricalGraph::PointList pl;
    CPMP<double>::DemandList dl;

    ifs >> problemNum >> optima;
    ifs >> vertexNum >> medianNum >> medianCap;

    while (vertexNum--) {
        ifs >> nodeSeqNum >> x >> y >> demand;
        pl.push_back( GeometricalGraph::Point( x, y ) );
        dl.push_back( demand );
    }

    ifs.close();


    GeometricalGraph gg( pl );
    UndirectedGraph<double> dug( gg );
    dug.getDistSeqTable();
    // UndirectedGraph<unsigned> uug( gg );

    // for each instance, run some times for judging average performance
    const int runTime = 2;
    const int maxIterCountBase = 12800;
    for (int i = 1; i <= runTime; i++) {
        {
            CPMP<double> cpmp( dug, dl, medianNum, medianCap, i * maxIterCountBase );
            cpmp.solve();
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
    ofstream ofs( "log.csv", ios::app );
    // CPMP<double>::initResultSheet( ofs ); // call if log.csv is not exist

    solve_pmedcap1( ofs, 1 );

    ofs.close();
    return 0;
}