#include "mQuickHull.h"

#include <maya/MGlobal.h>
#include <maya/MVectorArray.h>
#include <maya/MBoundingBox.h>

#include <iostream>
using namespace std;

#define HULL_TOLERANCE 1e-10

// Perpendicular Distance Above/Below Face
static double distanceOutsideFace(const MPoint &faceCenter, const MVector &faceNormal, const MPoint &pt) {
  return (MVector(pt - faceCenter) * faceNormal);
}

// Add Triangle to List
static void addTriangle(const int p0, const int p1, const int p2, const MPointArray &points, MIntArray &triangleVertices, MVectorArray &triangleNormals) {
  triangleVertices.append(p0);
  triangleVertices.append(p1);
  triangleVertices.append(p2);

  MVector v01 = points[p1] - points[p0];
  MVector v02 = points[p2] - points[p0];
  triangleNormals.append((v01 ^ v02).normal());
}

// Remove Triangles from List (Copy Triangles Not Being Removed)
static void removeTriangles(const MIntArray &remove, MIntArray &triangleVertices, MVectorArray &triangleNormals) {
  MIntArray vList;
  MVectorArray nList;
  for(unsigned t = 0; 3*t < triangleVertices.length(); t++) {
	bool removeTri = false;
	for(unsigned r = 0; r < remove.length(); r++) {
	  if(remove[r] == t) {
		removeTri = true;
		break;
	  }
	}
	if(removeTri == false) {
	  vList.append(triangleVertices[3*t+0]);
	  vList.append(triangleVertices[3*t+1]);
	  vList.append(triangleVertices[3*t+2]);
	  nList.append(triangleNormals[t]);
	}
  }
  triangleVertices = vList;
  triangleNormals = nList;
} 

static bool triangleContainsEdge(const MIntArray &triangleVertices, int tri, int e0, int e1, int &other_vertex) {
  int v0 = triangleVertices[3*tri+0];
  int v1 = triangleVertices[3*tri+1];
  int v2 = triangleVertices[3*tri+2];

  if(e0 == v0) {
	if(e1 == v1) {
	  other_vertex = v2;
	  return true;
	} else if(e1 == v2) {
	  other_vertex = v1;
	  return true;
	}
  } else if(e0 == v1) {
	if(e1 == v0) {
	  other_vertex = v2;
	  return true;
	} else if(e1 == v2) {
	  other_vertex = v0;
	  return true;
	}
  } else if(e0 == v2) {
	if(e1 == v0) {
	  other_vertex = v1;
	  return true;
	} else if(e1 == v1) {
	  other_vertex = v0;
	  return true;
	}
  }

  return false;
}

// Returns List of Triangle's Neighbour (Sharing 1 Edge with Triangle)
static MIntArray getTriangleNeighbours(const MIntArray &triangleVertices, int tri) {
  tri *= 3;
  int e0 = triangleVertices[tri+0];
  int e1 = triangleVertices[tri+1];
  int e2 = triangleVertices[tri+2];
  MIntArray neighbours;
  for(unsigned t = 0; t < triangleVertices.length(); t+=3) {
	if(t != tri) {
	  int count = 0;
	  count += (e0 == triangleVertices[t+0]) || (e0 == triangleVertices[t+1]) || (e0 == triangleVertices[t+2]);
	  count += (e1 == triangleVertices[t+0]) || (e1 == triangleVertices[t+1]) || (e1 == triangleVertices[t+2]);
	  count += (e2 == triangleVertices[t+0]) || (e2 == triangleVertices[t+1]) || (e2 == triangleVertices[t+2]);
	  if(count == 2) {
	    neighbours.append(t / 3);
	  }
	}
  }
  return neighbours;
}

// Determines The 'Ridge Edges' That Mark The Boundary of the Visible Face Set (To Be Replaced) & the Vertex from that edge being lost
static void getRidgeEdge(const MPointArray &points, const MIntArray &triangleVertices, const MVectorArray &triangleNormals, const MIntArray &visibleFaceSet, MIntArray &edges) {
  MStatus status;

  // Ridge Edges: Edge IDs That Are Used Exactly Once by visibleFaceSet
  for(unsigned i = 0; i < visibleFaceSet.length(); i++) {
	unsigned tri = 3*visibleFaceSet[i];
	int e0 = triangleVertices[tri+0];
	int e1 = triangleVertices[tri+1];
	int e2 = triangleVertices[tri+2];

	int count0 = 0, count1 = 0, count2 = 0;
	int other_vertex0, other_vertex1, other_vertex2;
	for(unsigned j = 0; j < visibleFaceSet.length(); j++) {
	  if(triangleContainsEdge(triangleVertices, visibleFaceSet[j], e0, e1, other_vertex0)) {
		count0++;
	  }
	  if(triangleContainsEdge(triangleVertices, visibleFaceSet[j], e1, e2, other_vertex1)) {
		count1++;
	  }
	  if(triangleContainsEdge(triangleVertices, visibleFaceSet[j], e2, e0, other_vertex2)) {
		count2++;
	  }
	}
	if(count0 == 1) {
	  edges.append(e0);
	  edges.append(e1);
	  edges.append(other_vertex0);
	}
	if(count1 == 1) {
	  edges.append(e1);
	  edges.append(e2);
	  edges.append(other_vertex1);
	}
	if(count2 == 1) {
	  edges.append(e2);
	  edges.append(e0);
	  edges.append(other_vertex2);
	}
  }
}

// Recursively Builds List of Neighbouring Faces That All Have Specified Point Above It
static void getNeighbouringOutsideFaces(const MPointArray &points, const MIntArray &triangleVertices, const MVectorArray &triangleNormals, MIntArray &visibleFaceSet, const MPoint &pt, int tri) {
  MPoint faceCenter = points[triangleVertices[3*tri]];
  MVector faceNormal = triangleNormals[tri];

  double len = distanceOutsideFace(faceCenter, faceNormal, pt);
  if(len > HULL_TOLERANCE) {
    for(unsigned f = 0; f < visibleFaceSet.length(); f++) {
	  if(visibleFaceSet[f] == tri) {
	    return;
	  }
    }
    visibleFaceSet.append(tri);

	MIntArray neighbours = getTriangleNeighbours(triangleVertices, tri);
	for(unsigned n = 0; n < neighbours.length(); n++) {
	  getNeighbouringOutsideFaces(points, triangleVertices, triangleNormals, visibleFaceSet, pt, neighbours[n]);
    }
  }
}

// Expands Hull to Include Point Furthest Above Any Hull Mesh Face
static MStatus expandHull(const MPointArray &points, MIntArray &unassignedPoints, MIntArray &triangleVertices, MVectorArray &triangleNormals) {
  MStatus status;

  // Determine Point Highest Above Any Face
  int maxFace = -1;
  int maxPoint = -1;
  double maxDistance = HULL_TOLERANCE;

  for(unsigned i = 0; i < unassignedPoints.length(); i++) {
	unsigned p = unassignedPoints[i];
	for(unsigned t = 0; 3*t < triangleVertices.length(); t++) {
	  MPoint faceCenter = points[triangleVertices[3*t]];
	  MVector faceNormal = triangleNormals[t];
		
	  double len = distanceOutsideFace(faceCenter, faceNormal, points[p]);
	  if(len > maxDistance) {
	    maxFace = t;
		maxPoint = p;
	    maxDistance = len;
	  }
	}
  }

  // All Points Inside Hull ?
  if(maxPoint < 0) {
	return MStatus::kFailure;
  }
  
  // Build Visible Face Set of Faces that Point is Above
  MIntArray visibleFaceSet;
  getNeighbouringOutsideFaces(points, triangleVertices, triangleNormals, visibleFaceSet, points[maxPoint], maxFace);

  // Determine 'Ridge' Edges Around Visible Set
  MIntArray edges;
  getRidgeEdge(points, triangleVertices, triangleNormals, visibleFaceSet, edges);

  // Remove Old Faces
  removeTriangles(visibleFaceSet, triangleVertices, triangleNormals);

  // Build New triangleVertices using Ridge Edges
  for(unsigned r = 0; 3*r < edges.length(); r++) {
	int p0 = edges[3*r];
	int p1 = edges[3*r+1];
	int oldp2 = edges[3*r+2];
	int p2 = maxPoint;

	MPoint pt0 = points[p0];
	MPoint pt1 = points[p1];
	MPoint pt2 = points[p2];
	MVector shift = points[p2] - points[oldp2];

    MVector vec01(pt1 - pt0);
    MVector vec02(pt2 - pt0);
    MVector vx = (vec01 ^ vec02).normal();

	if(vx*shift > 0.0) {
	  addTriangle(p0, p1, p2, points, triangleVertices, triangleNormals);
	} else {
	  addTriangle(p0, p2, p1, points, triangleVertices, triangleNormals);
	}
  }

  // Assign Points Now Inside Mesh
  for(unsigned i = 0; i < unassignedPoints.length(); ) {
	unsigned p = unassignedPoints[i];
	if(p == maxPoint) {
	  unassignedPoints.remove(i);
	} else {
	  bool outside = false;
	  for(unsigned t = 0; 3*t < triangleVertices.length(); t++) {
	    MPoint faceCenter = points[triangleVertices[3*t]];
	    MVector faceNormal = triangleNormals[t];
	    if(distanceOutsideFace(faceCenter, faceNormal, points[p]) > HULL_TOLERANCE) {
		  outside = true;
		  break;
	    }
	  }
	  if(outside) {
		i++;
	  } else {
		unassignedPoints.remove(i);
	  }
	}
  }

  return MStatus::kSuccess;
}

// Determine Point Cloud Bounding Box (inc. Low/High X/Y/Z Point Indexes)
static MBoundingBox rangePointCloud(const MPointArray &points, unsigned &loX, unsigned &hiX, unsigned &loY, unsigned &hiY, unsigned &loZ, unsigned &hiZ) {
  MBoundingBox bbox;
  bbox.expand(points[0]);

  loX = hiX = loY = hiY = loZ = hiZ = 0;
  for(unsigned p = 1; p < points.length(); p++) {
	MPoint pt = points[p];
	bbox.expand(pt);
	loX = (pt.x < points[loX].x) ? p : loX;
	hiX = (pt.x > points[hiX].x) ? p : hiX;
	loY = (pt.y < points[loY].y) ? p : loY;
	hiY = (pt.y > points[hiY].y) ? p : hiY;
	loZ = (pt.z < points[loZ].z) ? p : loZ;
	hiZ = (pt.z > points[hiZ].z) ? p : hiZ;
  }
  return bbox;
}

// Build Initial Simplex Hull (Tetrahedron) of Most Spread Out Points
static MStatus simplexHull(const MPointArray &points, MIntArray &unassignedPoints, MIntArray &triangleVertices, MVectorArray &triangleNormals) {
  MStatus status;

  unsigned p0, p1, p2, p3;

  // Determine Bounding Box of Points
  unsigned minX, maxX, minY, maxY, minZ, maxZ;
  MBoundingBox pointRange = rangePointCloud(points, minX, maxX, minY, maxY, minZ, maxZ);

  // Which Dimension has Greatest Spread? Set That As Simplex Vertices 0 & 1
  if(pointRange.width() > pointRange.height()) {
	if(pointRange.width() > pointRange.depth()) {
	  p0 = minX; p1 = maxX;
	} else {
	  p0 = minZ; p1 = maxZ;
	}
  } else {
	if(pointRange.height() > pointRange.depth()) {
	  p0 = minY; p1 = maxY;
	} else {
	  p0 = minZ; p1 = maxZ;
	}
  }

  // Check We Have A Decent Setup Line
  if(MVector(points[p1] - points[p0]).length() < HULL_TOLERANCE) {
	status = MStatus::kFailure;
	status.perror("Insufficient Range of Points for Hull Creation");
	return status;
  }

  // Find Vertex 2: Furthest from p0 -> p1 Line
  MVector dir01(points[p1] - points[p0]);
  dir01.normalize();

  int p2candidate = -1;
  double p2distance = HULL_TOLERANCE;
  for(unsigned i = 0; i < points.length(); i++) {
	if((i != p0) && (i != p1)) {
	  double u = dir01*MVector(points[i] - points[p0]);
	  double distance = MVector(points[i] - (points[p0]+u*dir01)).length();
	  if(distance > p2distance) {
		p2candidate = i;
		p2distance = distance;
	  }
	}
  }

  if(p2candidate < 0) {
	status = MStatus::kFailure;
	status.perror("Insufficient Range of Points for Hull Creation (Colinear?)");
	return status;
  }

  p2 = p2candidate;

  // Find Vertex 3: Furthest from Plane of p0 -> p1 & p0 -> p2 Lines (via Cross Product to Normal)
  MVector dir02(points[p2] - points[p0]);
  dir02.normalize();

  MVector dirx = (dir01 ^ dir02).normal();

  int p3candidate = -1;
  double p3distance = HULL_TOLERANCE;
  for(unsigned i = 0; i < points.length(); i++) {
	if((i != p0) && (i != p1) && (i != p2)) {
	  double distance = dirx*MVector(points[i] - points[p0]);
	  if(distance > p3distance) {
		p3candidate = i;
		p3distance = distance;
	  }
	}
  }

  if(p3candidate < 0) {
	status = MStatus::kFailure;
	status.perror("Insufficient Range of Points for Hull Creation (Coplanar?)");
	return status;
  }

  p3 = p3candidate;

  // Setup Simplex Mesh
  if(dirx * (points[p3] - points[p0]) < 0.0) {
	addTriangle(p0, p1, p2, points, triangleVertices, triangleNormals);
	addTriangle(p3, p1, p0, points, triangleVertices, triangleNormals);
	addTriangle(p3, p2, p1, points, triangleVertices, triangleNormals);
	addTriangle(p3, p0, p2, points, triangleVertices, triangleNormals);
  } else {
	addTriangle(p0, p2, p1, points, triangleVertices, triangleNormals);
	addTriangle(p3, p0, p1, points, triangleVertices, triangleNormals);
	addTriangle(p3, p1, p2, points, triangleVertices, triangleNormals);
	addTriangle(p3, p2, p0, points, triangleVertices, triangleNormals);
  }

  // Assign Points Now Inside Mesh
  for(unsigned i = 0; i < unassignedPoints.length(); ) {
	unsigned p = unassignedPoints[i];
	if((p == p0) || (p == p1) || (p == p2) || (p == p3)) {
	  unassignedPoints.remove(i);
	} else {
	  bool outside = false;
	  for(unsigned t = 0; 3*t < triangleVertices.length(); t++) {
	    MPoint faceCenter = points[triangleVertices[3*t]];
	    MVector faceNormal = triangleNormals[t];
	    if(distanceOutsideFace(faceCenter, faceNormal, points[p]) > HULL_TOLERANCE) {
		  outside = true;
		  break;
	    }
	  }
	  if(outside) {
		i++;
	  } else {
		unassignedPoints.remove(i);
	  }
	}
  }

  return MStatus::kSuccess;
}

// Builds Convex Hull to Wrap Points Cloud. Returns Triangle Mesh (VertexIDs: [Size = 3 * TriangleCount])
MStatus buildConvexHull(const MPointArray &points, MIntArray &triangleMesh) {
  MStatus status;

  if(points.length() < 4) {
    status = MStatus::kFailure;
	status.perror(MString("Insufficient Points for Hull Creation: ") + points.length());
    return status;
  }

  // Unassigned - List of VertexIDs Not On Hull & Not Inside Hull
  MIntArray unassignedPoints(points.length());
  for(unsigned p = 0; p < unassignedPoints.length(); p++) {
	unassignedPoints[p] = p;
  }

  // TriangleNormals - Triangle Normals [Size = TriangleCount]
  MVectorArray triangleNormals;

  // Build Simplex
  status = simplexHull(points, unassignedPoints, triangleMesh, triangleNormals);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  // Loop & Expand Hull Until All Points Inside (within Tolerance)
  do {
	status = expandHull(points, unassignedPoints, triangleMesh, triangleNormals);
  } while(status);

  return MStatus::kSuccess;
}