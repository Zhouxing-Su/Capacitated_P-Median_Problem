#include "cpmp_a.h"

using namespace std;

const int ACCURACY = 1;
typedef int DistType;
typedef CPMP<DistType, ACCURACY> Problem;

int solve_pmedcap1( ofstream &csvFile, int group, int instanceNum )
{
    ostringstream fname;
    fname << "../Instances/pmedcap" << group << "(" << instanceNum << ").txt";
    ifstream ifs( fname.str() );

    unsigned problemNum, optima;
    unsigned vertexNum, medianNum, medianCap;

    int nodeSeqNum;
    GeometricalGraph::Coord x, y;
    int demand;

    GeometricalGraph::PointList pl;
    Problem::DemandList dl;

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
    UndirectedGraph<DistType, ACCURACY> uug( gg );
    uug.getDistSeqTable();

    Problem cpmp( uug, dl, medianNum, medianCap );
    cpmp.solve();
    cpmp.printOptima( cout );
    if (!cpmp.check()) {
        csvFile << "[LogicError] ";
    }
    cpmp.appendResultToSheet( fname.str(), csvFile );

    return 0;
}

int main()
{
    ofstream ofs( "../Instances/log.csv", ios::app );
    //CPMP<DistType>::initResultSheet( ofs ); // call if log.csv is not exist

    solve_pmedcap1( ofs, 1, 1 );
    //for (int i = 1; i <= 4; i++) {
    //    for (int j = 1; j <= 20; j++) {
    //        solve_pmedcap1( ofs, i, j );
    //    }
    //}

    ofs.close();
    return 0;
}