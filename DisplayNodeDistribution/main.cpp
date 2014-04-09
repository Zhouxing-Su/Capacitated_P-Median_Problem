#include "DisplayNodeDistribution.h"

using namespace std;
using namespace cv;

int main( void )
{
    const int MAX_STR_LEN = 256;

    ifstream logFile( "../Instances/log.csv" );

    char buf[MAX_STR_LEN];
    char date[MAX_STR_LEN], instance[MAX_STR_LEN], algorithm[MAX_STR_LEN],
        totalIter[MAX_STR_LEN], randSeed[MAX_STR_LEN], duration[MAX_STR_LEN],
        iterCount[MAX_STR_LEN], totalDist[MAX_STR_LEN];

    logFile.getline( buf, MAX_STR_LEN, '\n' );  // read sheet header

    for (int testCase = 1; true; testCase++) {
        // read calculation infomation
        logFile.getline( date, MAX_STR_LEN, ',' );
        logFile.getline( instance, MAX_STR_LEN, ',' );
        logFile.getline( algorithm, MAX_STR_LEN, ',' );
        logFile.getline( totalIter, MAX_STR_LEN, ',' );
        logFile.getline( randSeed, MAX_STR_LEN, ',' );
        logFile.getline( duration, MAX_STR_LEN, ',' );
        logFile.getline( iterCount, MAX_STR_LEN, ',' );
        logFile.getline( totalDist, MAX_STR_LEN, ',' );

        if (logFile.eof()) {
            break;
        }

        // read instance to get the median number, node number and coordinate of nodes
        ifstream instFile( instance + 1 );    // skip a space

        unsigned problemNum, optima;
        unsigned vertexNum, medianNum, medianCap;

        int nodeSeqNum;
        GeometricalGraph::Coord x, y;
        int demand;

        GeometricalGraph::PointList pl;

        instFile >> problemNum >> optima;
        instFile >> vertexNum >> medianNum >> medianCap;

        for (unsigned i = 0; i < vertexNum; i++) {
            instFile >> nodeSeqNum >> x >> y >> demand;
            pl.push_back( GeometricalGraph::Point( x, y ) );
        }

        instFile.close();

        // read the result
        char c;
        vector<int> medians, assignment;
        for (unsigned i = 0; i < medianNum; i++) {
            int median;
            logFile >> median >> c;
            medians.push_back( median );
        }
        logFile >> c;   // read the comma
        for (unsigned i = 0; i < vertexNum; i++) {
            int node;
            logFile >> node >> c;
            assignment.push_back( node );
        }

        // calculate the apropriate width and height, then regularize the coordinates
        GeometricalGraph gg( pl );
        const double GRAPH_AMP = 8;
        gg.stretch( GRAPH_AMP );
        const GeometricalGraph::Coord MARGIN_LEN = 50;
        gg.shift( GeometricalGraph::Point( MARGIN_LEN, MARGIN_LEN ) );
        GeometricalGraph::Rectangle minCoverRect( gg.getMinCoverRect() );

        // Create a white background image with margin
        int width = minCoverRect.right + MARGIN_LEN;
        int height = minCoverRect.top + MARGIN_LEN;
        Mat distributionImg( width, height, CV_8UC3, RGBcolor::WHITE );

        // draw assignments as lines and customers
        for (int i = 0; i < gg.vertexNum; i++) {
            Point customer( gg.point( i ).x, gg.point( i ).y );
            Point median( gg.point( assignment[i] ).x, gg.point( assignment[i] ).y );
            circle( distributionImg, customer, 2, RGBcolor::BLACK, 1, CV_AA );
            line( distributionImg, customer, median, RGBcolor::BLACK, 1, CV_AA );
        }

        // draw medians as circles
        for (int i = 0; i < medianNum; i++) {
            Point median( gg.point( medians[i] ).x, gg.point( medians[i] ).y );
            circle( distributionImg, median, 4, RGBcolor::RED, 1, CV_AA );
        }

        // imshow( date, distributionImg );
        sprintf( buf, "../Instances/results/%d.png", testCase );
        imwrite( buf, distributionImg );
        waitKey();
    }

    logFile.close();

    waitKey();
    return 0;
}