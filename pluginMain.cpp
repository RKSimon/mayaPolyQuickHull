#include <maya/MFnPlugin.h>

#include "polyQuickHullNode.h"

MStatus initializePlugin(MObject obj) { 
  MStatus status;
  MFnPlugin plugin(obj, "Simon Pilgrim", "6.0", "Any");

  status = plugin.registerNode("polyQuickHull", polyQuickHull::id, polyQuickHull::creator, polyQuickHull::initialize);
  if(!status) {
    status.perror("registerNode");
    return status;
  }

  return status;
}

MStatus uninitializePlugin(MObject obj) {
  MStatus status;
  MFnPlugin plugin(obj);

  status = plugin.deregisterNode(polyQuickHull::id);
  if(!status) {
    status.perror("deregisterNode");
    return status;
  }

  return status;
}