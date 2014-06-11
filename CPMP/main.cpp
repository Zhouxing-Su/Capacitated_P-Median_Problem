#include "cpmp.h"

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

    // for each instance, run some times for judging average performance
    const int runTime = 4;
    const int maxIterCountBase = 3200;
    const int tabuTenureAssign = uug.vertexNum * medianNum / 5;
    const int tabuTenureOpenMedian = uug.vertexNum / 4;
    const int tabuTenureCloseMedian = medianNum / 3;
    const int maxNoImproveCount = uug.vertexNum * medianNum * 4;
    const int randomWalkStep = uug.vertexNum * 2;
    const DistType demandDistributionDamping = (uug.DistMultiplication * (gg.getMinCoverRect().right - gg.getMinCoverRect().left) / 4);
    for (int i = 1; i <= runTime; i++) {
        {
            Problem cpmp( uug, dl, medianNum, medianCap );
            cpmp.solve( maxIterCountBase, maxNoImproveCount, tabuTenureAssign,
                tabuTenureOpenMedian, tabuTenureCloseMedian, demandDistributionDamping, randomWalkStep );
            cpmp.printOptima( cout );
            if (!cpmp.check()) {
                csvFile << "[LogicError] ";
            }
            cpmp.appendResultToSheet( fname.str(), csvFile );
        }
        ;
        {
            Problem cpmp( uug, dl, medianNum, medianCap );
            cpmp.solve_ShiftSwapTabuRelocateTabu( maxIterCountBase, maxNoImproveCount, tabuTenureAssign,
                tabuTenureOpenMedian, tabuTenureCloseMedian, demandDistributionDamping );
            cpmp.printOptima( cout );
            if (!cpmp.check()) {
                csvFile << "[LogicError] ";
            }
            cpmp.appendResultToSheet( fname.str(), csvFile );
        }
        ;
        {
            Problem cpmp( uug, dl, medianNum, medianCap );
            cpmp.solve_ShiftSwapTabuRelocate( maxIterCountBase, maxNoImproveCount, tabuTenureAssign, demandDistributionDamping );
            cpmp.printOptima( cout );
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

    //solve_pmedcap1( ofs, 1, 3 );
    //1.加上服务节点的扰动
    //2.换服务节点时，对被删服务节点的所有用户，按距离由远到近分配给有剩余容量的服务节点，
    //    而不是直接给新增的服务节点
    for (int i = 1; i <= 4; i++) {
        for (int j = 1; j <= 20; j++) {
            solve_pmedcap1( ofs, i, j );
        }
    }

    ofs.close();
    return 0;
}