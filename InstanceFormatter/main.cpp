#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

using namespace std;

void generateAllGroups();
void formatToCplexInput();


int main()
{
    //generateAllGroups();
    formatToCplexInput();
}


void generateAllGroups()
{
    int problemNum, optima;
    int vertexNum, medianNum, medianCap;

    char name[512];
    char buf[8192];
    FILE *pf1;
    FILE *pf2;
    int i = 1;
    int n = 0;
    for (; i <= 40; i++) {
        sprintf( name, "D:/workspace/cpp/CapacitatedPMedianProblem/Instances/pmedcap1(%d).txt", i );
        pf1 = fopen( name, "r" );
        sprintf( name, "D:/workspace/cpp/CapacitatedPMedianProblem/Instances/pmedcap2(%d).txt", i );
        pf2 = fopen( name, "w" );

        do {
            fscanf( pf1, "%d%d%d%d%d", &problemNum, &optima, &vertexNum, &medianNum, &medianCap );
            medianNum = vertexNum * 2 / 5;
            fprintf( pf2, " %d %d\n %d %d %d", problemNum, optima, vertexNum, vertexNum / 4, 12 * vertexNum / medianNum );

            n = fread( buf, 1, 4096, pf1 );
            fwrite( buf, n, 1, pf2 );
            printf( "%d\n", n );
        } while (n == 4096);

        fclose( pf2 );
        fclose( pf1 );
    }
}

void formatToCplexInput()
{
    int problemNum, optima;
    int vertexNum, medianNum, medianCap;

    char name[512];
    FILE *pf1;
    FILE *pf2;
    int i = 1;

    for (int j = 1; j <= 4; j++) {
        for (int i = 1; i <= 40; i++) {
            sprintf( name, "D:/workspace/cpp/CapacitatedPMedianProblem/Instances/pmedcap%d(%d).txt", j, i );
            pf1 = fopen( name, "r" );
            sprintf( name, "D:/workspace/cplex/cpmp/pmedcap%d(%d).dat", j, i );
            pf2 = fopen( name, "w" );

            int nodeSeqNum;
            int x, y;
            int demand;

            vector<int> px, py, d;

            fscanf( pf1, "%d%d%d%d%d", &problemNum, &optima, &vertexNum, &medianNum, &medianCap );

            for (int k = 0; k < vertexNum; k++) {
                fscanf( pf1, "%d%d%d%d", &nodeSeqNum, &x, &y, &demand );
                px.push_back( x );
                py.push_back( y );
                d.push_back( demand );
            }

            fprintf( pf2, "n = %d;\n", vertexNum );
            fprintf( pf2, "p = %d;\n", medianNum );
            fprintf( pf2, "Q = %d;\n", medianCap );

            fprintf( pf2, "q = [ %d", d[0] );
            for (int k = 1; k < vertexNum; k++) {
                fprintf( pf2, ", %d", d[k] );
            }
            fprintf( pf2, " ];\n" );

            fprintf( pf2, "px = [ %d", px[0] );
            for (int k = 1; k < vertexNum; k++) {
                fprintf( pf2, ", %d", px[k] );
            }
            fprintf( pf2, " ];\n" );


            fprintf( pf2, "py = [ %d", py[0] );
            for (int k = 1; k < vertexNum; k++) {
                fprintf( pf2, ", %d", py[k] );
            }
            fprintf( pf2, " ];\n" );

            fclose( pf1 );
            fclose( pf2 );
        }
    }
}