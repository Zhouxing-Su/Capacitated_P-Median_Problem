#include "cpmp.h"


int solve_pmedcap1( ofstream &csvFile, int instanceNum )
{
    ostringstream fname;
    fname << "pmedcap1(" << instanceNum << ").txt";
    ifstream ifs( fname.str() );

    unsigned problemNum, bestSolutionValue;
    unsigned nodeNum, medianNum, medianCap;

    int nodeSeqNum;
    GeometricalGraph::Coord x, y;
    int demand;

    GeometricalGraph::PointList pl;
    CPMP<double>::Demand dl;

    ifs >> problemNum >> bestSolutionValue;
    ifs >> nodeNum >> medianNum >> medianCap;

    while (nodeNum--) {
        ifs >> nodeSeqNum >> x >> y >> demand;
        pl.push_back( GeometricalGraph::Point( x, y ) );
        dl.push_back( demand );
    }

    ifs.close();


    GeometricalGraph gg( pl );
    UndirectedGraph<double> dug( gg );
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
            //cpmp.appendResultToSheet( fname.str(), csvFile );
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