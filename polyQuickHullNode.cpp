#include "polyQuickHullNode.h"

#include <math.h>
#include <maya/MGlobal.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MArrayDataHandle.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MPointArray.h>

#include "mQuickHull.h"

#error Choose A Suitable Maya API ID
MTypeId polyQuickHull::id(0x00000);

MObject polyQuickHull::outMesh;
MObject polyQuickHull::inMesh;
MObject polyQuickHull::inPoints;
MObject polyQuickHull::inPointsX;
MObject polyQuickHull::inPointsY;
MObject polyQuickHull::inPointsZ;

#include <iostream>
using namespace std;

void *polyQuickHull::creator() {
  return new polyQuickHull();
}

MStatus polyQuickHull::initialize() {
  MStatus status;
  MFnNumericAttribute nattr;
  MFnTypedAttribute tattr;

  // Create outMesh Attribute - Polygon Mesh of Convex Hull
  outMesh = tattr.create("outMesh", "om", MFnData::Type::kMesh, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  CHECK_MSTATUS(tattr.setWritable(false));

  // Create inMesh Attribute - 3D Coordinates Array of Vertices to Create Convex Hull Around
  inMesh = tattr.create("inMesh", "im", MFnData::Type::kMesh, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  CHECK_MSTATUS(tattr.setArray(true));
  CHECK_MSTATUS(tattr.setHidden(true));
  CHECK_MSTATUS(tattr.setIndexMatters(false));
  CHECK_MSTATUS(tattr.setDisconnectBehavior(MFnAttribute::DisconnectBehavior::kDelete));

  // Create inPoints Attribute - 3D Coordinates Array of Points to Create Convex Hull Around
  inPointsX = nattr.create("inPointsX", "ipx", MFnNumericData::Type::kDouble, 0.0, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  inPointsY = nattr.create("inPointsY", "ipy", MFnNumericData::Type::kDouble, 0.0, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  inPointsZ = nattr.create("inPointsZ", "ipz", MFnNumericData::Type::kDouble, 0.0, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  inPoints = nattr.create("inPoints", "ip", inPointsX, inPointsY, inPointsZ, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  CHECK_MSTATUS(nattr.setArray(true));
  CHECK_MSTATUS(nattr.setHidden(true));
  CHECK_MSTATUS(nattr.setIndexMatters(false));
  CHECK_MSTATUS(nattr.setDisconnectBehavior(MFnAttribute::DisconnectBehavior::kDelete));

  status = addAttribute(polyQuickHull::outMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = addAttribute(polyQuickHull::inMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = addAttribute(polyQuickHull::inPoints);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  status = attributeAffects(polyQuickHull::inMesh, polyQuickHull::outMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = attributeAffects(polyQuickHull::inPoints, polyQuickHull::outMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  return MStatus::kSuccess;
}

polyQuickHull::polyQuickHull() {
}

polyQuickHull::~polyQuickHull() {
}

MStatus polyQuickHull::compute(const MPlug &plug, MDataBlock &data) {
  MStatus status;

  // Attribute Access
  MArrayDataHandle adhInMesh = data.inputArrayValue(inMesh, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  MArrayDataHandle adhInPoints = data.inputArrayValue(inPoints, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  MDataHandle dhOutMesh = data.outputValue(outMesh, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  // Build Point Cloud
  MPointArray points;

  for(unsigned i = 0; i < adhInMesh.elementCount(); i++) {
	MFnMesh mesh(adhInMesh.inputValue().asMesh());
	MPointArray meshPoints;
	mesh.getPoints(meshPoints);
	for(unsigned p = 0; p < meshPoints.length(); p++) {
	  points.append(meshPoints[p]);
	}
	adhInMesh.next();
  }

  for(unsigned i = 0; i < adhInPoints.elementCount(); i++) {
    points.append(MPoint(adhInPoints.inputValue().asVector()));
    adhInPoints.next();
  }

  if(points.length() < 4) {
    status = MStatus::kFailure;
	status.perror(MString("Insufficient Points for Hull Creation: ") + points.length());
    return status;
  }

  // Triangles - Triangle Vertices Indices [Size = 3*TriangleCount]
  MIntArray triangles;

  // Get Convex Hull
  status = buildConvexHull(points, triangles);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  // Construct Polygon Mesh from Triangle List
  MFnMeshData fnMeshData;
  MObject meshData = fnMeshData.create(&status);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  MFnMesh fnMesh;
  fnMesh.create(points.length(), triangles.length()/3, points, MIntArray(triangles.length()/3,3), triangles, meshData, &status);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  status = dhOutMesh.set(meshData);
  CHECK_MSTATUS_AND_RETURN_IT(status);
  status = data.setClean(outMesh);
  CHECK_MSTATUS_AND_RETURN_IT(status);

  return MStatus::kSuccess;
}