#ifndef _M_QuickHull_H
#define _M_QuickHull_H

#include <maya/MIntArray.h>
#include <maya/MPointArray.h>

MStatus buildConvexHull(const MPointArray &points, MIntArray &triangleMesh);

#endif