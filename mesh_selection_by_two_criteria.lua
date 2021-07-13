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
--! Script to create a sample 3D sphere geometry using plugin ProMesh to
--! illustrate the selection by two criteria.
--!
--! To execute type:
--! ugshell -ex mesh_selection_by_two_criteria.lua
--!--------------------------------------------------------------------------------

-----------------------------
-- Vector subtraction
-----------------------------
function Vec3dSubtract(vec1, vec2)
	result = Vec3d()
	result:set_coord(0, vec1:coord(0) - vec2:coord(0))
	result:set_coord(1, vec1:coord(1) - vec2:coord(1))
	result:set_coord(2, vec1:coord(2) - vec2:coord(2))
	
	return result
end

-----------------------------
-- Euklidian norm
-----------------------------
function Vec3dNorm(vec)
	normSq = vec:coord(0)*vec:coord(0) + vec:coord(1)*vec:coord(1) + vec:coord(2)*vec:coord(2)
	norm = math.sqrt(normSq)

	return norm
end


-----------------------------
-- Start
-----------------------------
print(" > generating test mesh")
mesh = Mesh()

-- mesh, center coordinate, radius, number of refinements, subset index
CreateSphere(mesh, MakeVec(0.0, 0.0, 0.0), 1, 2, 0)
ClearSelection(mesh)

-- Subset management
AssignSubsetColors(mesh)

-----------------------------
-- first selection criterion
-----------------------------
SelectLongEdges(mesh, 0.28)
SaveMesh(mesh, "1_criterion_selection.ugx")

-----------------------------
-- second selection criterion
-----------------------------
targetDirection = Vec3d(0.14655714625849, 0.14203952192021, -0.23713444394043)

-- edge iterator pointing to first of previously selected edges 
eIter = mesh:edge_selection_begin()

-- loop over previously selected edges as long as iterator does not point to end of selection
while eIter:unequal(mesh:edge_selection_end()) do 
	-- get edge from iterator
	edge = eIter:value() 

	-- increment iterator to next edge
	eIter:advance()

	-- get vertices from edge
	vrt1 = edge:vertex(0)
	vrt2 = edge:vertex(1)

	-- calculate edge direction
	dirVec = Vec3dSubtract(mesh:position(vrt2), mesh:position(vrt1))

	-- debug output
	--print(dirVec:coord(0) .. "; " .. dirVec:coord(1) .. "; " .. dirVec:coord(2))

	-- if direction of currently investigated edge does not meet the target direction, then deselect
	if Vec3dNorm(Vec3dSubtract(targetDirection, dirVec)) > 1e-10 then				
		mesh:selector():deselect(edge)
	end
end

SaveMesh(mesh, "2_criterion_selection.ugx")