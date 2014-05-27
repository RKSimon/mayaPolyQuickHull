polyQuickHull
=============

Maya plugin that creates a polygonal mesh of the convex hull of irregular point cloud, using
"The Quickhull Algorithm for Convex Hulls (1996)" by Barber, Dobkin & Huhdanpaa.

http://citeseer.ist.psu.edu/barber96quickhull.html

Actual use is quite simple, just connect the meshes to the inMesh array attribute, any individual points to the inPoints array
attribute and connect the outMesh to a mesh node. Transform space is largely ignored so be consistant with what you connect!

testHull.mel
Script to demonstrate adding 100 random locator positions to the polyQuickHull node and outputing the resulting
convex hull. Moving / Deleting locators will update the mesh.
