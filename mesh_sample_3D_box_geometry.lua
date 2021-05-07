-- Copyright (c) 2010-2021:  G-CSC, Goethe University Frankfurt
-- Author: Martin Stepniewski
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.

--!--------------------------------------------------------------------------------
--! Script to create a sample 3D box geometry using plugin ProMesh consisting of
--! (1) a thin box on top of the larger box, 
--! (2) a tall box inside the large box that sits on the bottom of the large box, 
--! (3) a box inside a box sitting on the bottom of the large box, and 
--! (4) a sphere inside a sphere sitting on the bottom of the large box
--!
--! To execute type:
--! ugshell -ex mesh_sample_3D_box_geometry.lua
--!--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"ProMesh"})

print("Meshing custom 3D box geometry")

-- User parameter
outerSphereRadius = util.GetParamNumber("-outerSphereRadius", 0.2)
innerSphereRadius = util.GetParamNumber("-innerSphereRadius", 0.1)
outerBoxSubsetName = util.GetParam("-outerBoxSubsetName", "largeBox")

--	Create a mesh object first. All promesh-methods operate on a mesh
obj = Mesh()

-- Create nested boxes with parameters: 
-- Mesh* mesh, const Vec3d* min corner, const Vec3d* max corner, integer subset, bool fill 
CreateBox(obj, MakeVec(0.0, 0.0, 0.0), MakeVec(2.0, 1.0, 1.0), 0, false)
CreateBox(obj, MakeVec(0.0, 0.0, 1.0), MakeVec(2.0, 1.0, 1.3), 1, false)
CreateBox(obj, MakeVec(0.25, 0.25, 0.0), MakeVec(0.35, 0.35, 0.5), 2, false)
CreateBox(obj, MakeVec(0.75, 0.25, 0.0), MakeVec(1.25, 0.75, 0.4), 3, false)
CreateBox(obj, MakeVec(0.9, 0.4, 0.0), MakeVec(1.1, 0.6, 0.2), 4, false)
ClearSelection(obj)

-- Select intersecting bottom faces with parameters:
-- Mesh* mesh, const Vec3d* min, const Vec3d* max
SelectFacesInBox(obj, MakeVec(0.0, 0.0, -0.01), MakeVec(2.0, 1.0, 0.01))

-- Erase previously selected bottom faces with parameters:
-- Mesh* mesh, bool erase unused vertices, bool erase unused edges, bool erase unused faces
EraseSelectedElements(obj, false, false, false)
ClearSelection(obj)

-- Mark crease edges with parameters:
-- Mesh* mesh, number min angle, bool clear marks
MarkCreaseEdges(obj, 25, false)

-- Create nested spheres with parameters: 
-- mesh, center coordinate, radius, number of refinements, subset index
CreateSphere(obj, MakeVec(1.6, 0.5, 0.2), outerSphereRadius, 2, 5)
CreateSphere(obj, MakeVec(1.6, 0.5, 0.1), innerSphereRadius, 2, 6)

-- Select everything and remove double vertices with parameters:
-- mesh, threshold
SelectAll(obj)
RemoveDoubles(obj, 0.0001)
ClearSelection(obj)

-- Select intersecting sphere triangles 
SelectIntersectingTriangles(obj)

-- Resolve self intersecting sphere triangles with parameters:
-- Mesh* mesh, number snap threshold
ResolveSelfIntersections(obj, 0.000001)

-- Retriangulate previous selection with parameters:
-- Mesh* mesh, number min angle
Retriangulate (obj, 20)
ClearSelection(obj)

-- Select bottom edges with parameters:
-- Mesh* mesh, const Vec3d* min, const Vec3d* max
SelectEdgesInBox(obj, MakeVec(0.0, 0.0, -0.001), MakeVec(2.0, 1.0, 0.001))

-- Triangulate bottom edge selection s.t. crease constraints with parameters:
-- Mesh* mesh, bool quality generation, number min angle, integer subset
TriangleFill(obj, true, 20, 0)
ClearSelection(obj)

-- Resolve self intersection of bottom plane and sphere vertex with parameters:
-- Mesh* mesh, number snap threshold
SelectAll(obj)
ResolveSelfIntersections(obj, 0.000001)
ClearSelection(obj)

-- Select bottom faces with parameters:
-- Mesh* mesh, const Vec3d* min, const Vec3d* max
SelectFacesInBox(obj, MakeVec(0.0, 0.0, -0.001), MakeVec(2.0, 1.0, 0.001))

-- Retriangulate selection with parameters:
-- Mesh* mesh, number min angle
Retriangulate (obj, 20)
ClearSelection(obj)

-- Clear marks
ClearMarks(obj)

-- Tetrahedralize mesh with parameters: 
-- Mesh* mesh, number quality (min radius-edge ratio), bool preserve outer, bool preserve all, bool separate volumes, bool append subsets at end, integer verbosity
Tetrahedralize(obj, 1.414, false, false, true, true, 0)

-- Subset management
AssignSubsetColors(obj)
obj:subset_handler():set_subset_name("largeBox", 0)
obj:subset_handler():set_subset_name("thinTopBox", 1)
obj:subset_handler():set_subset_name("tallInnerBox", 2)
obj:subset_handler():set_subset_name("outerNestedBox", 3)
obj:subset_handler():set_subset_name("innerNestedBox", 4)
obj:subset_handler():set_subset_name("outerSphere", 5)
obj:subset_handler():set_subset_name("innerSphere", 6)

--	Finally save the resulting mesh to a file
SaveMesh(obj, "3D_box_geometry.ugx")