#ifndef _polyQuickHullNode
#define _polyQuickHullNode

#include <maya/MPxNode.h>

class polyQuickHull : public MPxNode {
 public:
   static void *creator();
   static MStatus initialize();

   polyQuickHull();
   virtual ~polyQuickHull(); 

   virtual MStatus compute(const MPlug &plug, MDataBlock &data);

   static MTypeId id;

   static MObject outMesh;
   static MObject inMesh;
   static MObject inPoints;
   static MObject inPointsX;
   static MObject inPointsY;
   static MObject inPointsZ;
};

#endif
