#include "DisplayNodeDistribution.h"

using namespace std;
using namespace cv;

int main( void )
{
    const int MAX_STR_LEN = 256;

    ifstream logFile( "../CPMP/log.csv" );

    char c;
    char buf[MAX_STR_LEN];

    char date[MAX_STR_LEN], instance[MAX_STR_LEN], algorithm[MAX_STR_LEN],
        totalIter[MAX_STR_LEN], randSeed[MAX_STR_LEN], duration[MAX_STR_LEN],
        iterCount[MAX_STR_LEN], totalDist[MAX_STR_LEN];

    logFile.getline( buf, MAX_STR_LEN, '\n' );  // read sheet head

    while (!logFile.eof()) {
        // read calculation infomation
        logFile.getline( date, MAX_STR_LEN, ',' );
        logFile.getline( instance, MAX_STR_LEN, ',' );
        logFile.getline( algorithm, MAX_STR_LEN, ',' );
        logFile.getline( totalIter, MAX_STR_LEN, ',' );
        logFile.getline( randSeed, MAX_STR_LEN, ',' );
        logFile.getline( duration, MAX_STR_LEN, ',' );
        logFile.getline( iterCount, MAX_STR_LEN, ',' );
        logFile.getline( totalDist, MAX_STR_LEN, ',' );

        // read instance to get the median number, node number and coordinate of nodes
        ifstream instFile( instance );
        int medianNum, nodeNum;

        // read the result
        vector<int> medians, assignment;
        for (int i = 0; i < medianNum; i++) {
            int median;
            logFile >> median >> c;



        }

        for (int i = 0; i < nodeNum; i++) {
            int node;
            logFile >> node >> c;



        }

        /// Create white empty images
        int row, col;
        Mat distributionImage( row, col, CV_8UC3, RGBcolor::WHITE );

        //distributionImage.at<Vec3b>( row, col ) = Vec3b( 0, 255, 0 );
        //line( distributionImage, Point( row, col ), Point( row, col ), RGBcolor::BLACK );
        //waitKey( 1 );
        //imshow( date, distributionImage );

        instFile.close();
    }

    logFile.close();

    waitKey();
    return 0;
}