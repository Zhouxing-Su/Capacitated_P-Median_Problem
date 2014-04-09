/**
*   usage : display node distribution and assignment for CPMP.
*           this program can be only applied to instance with geometrical graph.
*
*   note :  1.
*/

#ifndef DISPLAY_NODE_DISTRIBUTION_H

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstring>

#include "../CPPutilibs/Graph.h"

namespace RGBcolor
{
    const cv::Scalar WHITE( 255, 255, 255 );
    const cv::Scalar BLACK( 0, 0, 0 );
    const cv::Scalar RED( 0, 0, 255 );
}

#define DISPLAY_NODE_DISTRIBUTION_H
#endif