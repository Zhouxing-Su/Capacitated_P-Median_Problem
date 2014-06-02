#include <iostream>
#include <fstream>

using namespace std;

const int GROUP_NUM = 4;
const int INST_NUM = 20;
const int MAX_STR_LEN = 1024;
const double MAX_DURATION = 36000;
const int MAX_DIST = 10000;
const int RUN_TIME = 4;

const int optima[4][40] = {
    { 713, 740, 751, 651, 664, 778, 787, 820, 715, 829, 1006, 966, 1026, 982, 1091, 954, 1034, 1043, 1031, 1005, 1288, 1256, 1279, 1220, 1193, 1264, 1323, 1233, 1219, 1201, 1378, 1424, 1367, 1385, 1437, 1382, 1458, 1382, 1374, 1416 },
    { 383, 412, 405, 384, 429, 482, 445, 403, 436, 461, 544, 504, 555, 544, 583, 534, 542, 508, 551, 565, 681, 660, 663, 594, 629, 653, 736, 644, 649, 630, 727, 878, 713, 827, 746, 701, 753, 746, 722, 764 },
    { 298, 336, 314, 303, 351, 390, 361, 353, 373, 390, 414, 391, 446, 447, 474, 447, 431, 456, 445, 460, 599, 561, 564, 505, 488, 540, 579, 503, 545, 502, 575, 777, 620, 710, 611, 580, 631, 608, 589, 630 },
    { 266, 298, 311, 277, 356, 370, 358, 312, 412, 458, 415, 377, 412, 421, 496, 428, 440, 450, 450, 486, 552, 601, 555, 487, 436, 512, 745, 471, 494, 444, 541, 816, 557, 851, 552, 551, 594, 592, 540, 588 }
};


ifstream logFile( "../Instances/log.csv" );

char c;
char buf[MAX_STR_LEN];
char date[MAX_STR_LEN], instance[MAX_STR_LEN], algorithm[MAX_STR_LEN];
char maxIterCount[MAX_STR_LEN], maxNoImproveCount[MAX_STR_LEN], randSeed[MAX_STR_LEN];
char iterCount[MAX_STR_LEN], moveCount[MAX_STR_LEN];
double duration;
int totalDist;

// read calculation infomation
bool readLine()
{
    logFile.getline( date, MAX_STR_LEN, ',' );
    logFile.getline( instance, MAX_STR_LEN, ',' );
    logFile.getline( algorithm, MAX_STR_LEN, ',' );
    logFile.getline( maxIterCount, MAX_STR_LEN, ',' );
    logFile.getline( maxNoImproveCount, MAX_STR_LEN, ',' );
    logFile.getline( randSeed, MAX_STR_LEN, ',' );
    logFile >> duration >> c;
    logFile.getline( iterCount, MAX_STR_LEN, ',' );
    logFile.getline( moveCount, MAX_STR_LEN, ',' );
    logFile >> totalDist >> c;
    logFile.getline( buf, MAX_STR_LEN );    // read the medians and assignment( no use )

    return logFile.eof();
}

int main()
{
    ofstream ofs( "out.csv" );
    ofs << "Instance, Optima, Method1_Dist, Method1_Duration, Method2_Dist, Method2_Duration, Method3_Dist, Method3_Duration" << endl;

    logFile.getline( buf, MAX_STR_LEN );  // read sheet header

    for (int k = 1; k <= GROUP_NUM; k++) {
        for (int j = 1; j <= INST_NUM; j++) {
            double minTime3 = MAX_DURATION;
            int minDist3 = MAX_DIST;
            double minTime2 = MAX_DURATION;
            int minDist2 = MAX_DIST;
            double minTime1 = MAX_DURATION;
            int minDist1 = MAX_DIST;

            for (int i = 0; i < RUN_TIME; i++) {
                if (readLine()) {
                    return 0;
                }
                if (minDist3 > totalDist) {
                    minDist3 = totalDist;
                    minTime3 = duration;
                } else if ((minDist3 == totalDist) && (minTime3 > duration)) {
                    minTime3 = duration;
                }

                if (readLine()) {
                    return 0;
                }
                if (minDist2 > totalDist) {
                    minDist2 = totalDist;
                    minTime2 = duration;
                } else if ((minDist2 == totalDist) && (minTime2 > duration)) {
                    minTime2 = duration;
                }

                if (readLine()) {
                    return 0;
                }
                if (minDist1 > totalDist) {
                    minDist1 = totalDist;
                    minTime1 = duration;
                } else if ((minDist1 == totalDist) && (minTime1 > duration)) {
                    minTime1 = duration;
                }
            }

            ofs << "pmedcap" << k << '(' << j << ')' << ", "
                << optima[k - 1][j - 1] << ", "
                << minDist1 << ", " << minTime1 << ", "
                << minDist2 << ", " << minTime2 << ", "
                << minDist3 << ", " << minTime3 << endl;
        }
    }

    return 0;
}
