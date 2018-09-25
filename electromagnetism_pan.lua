--------------------------------------------------------------------------------
--[[!
-- \file apps/electromagnetism/pan.lua
-- \author Dmitry Logashenko
-- \brief Simulation of electric field around a coil and a conducting disc.
--
-- The problem setting simplifies the induction cooking: The conducting disc
-- is heated by the eddy currents induced by the magnetic field of the coil.
]]--
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
-- Remark on the solvability of the problem:
-- As the matrix of the discretized problem is singular, the problem is
-- only solvable if the rhs is in the image of the matrix, whereas this
-- image is not the whole space. The image consists of the numerically
-- divergence-free functions which is computed (cf. the computation of
-- the loop source below). The computation of the rhs involves the solution
-- of the discretized Poisson equation with the same solver as used for
-- the projection. The accuracy of the computation of the rhs is essential
-- for the solvability of the eddy current problem, so if your solver for
-- the eddy current problem does not converge (in any way), try to set
-- a higher accuracy for the computation of the rhs.
--------------------------------------------------------------------------------

PluginRequired("Electromagnetism")

ug_load_script ("ug_util.lua")
ug_load_script ("util/load_balancing_util.lua")

-- constants
dim        = 3; -- the problem is formulated in 3d
numPreRefs = util.GetParamNumber ("-numPreRefs", 0, "number of refinements before parallel distribution")
numRefs    = util.GetParamNumber ("-numRefs",    2, "number of refinements")

-- Geometry
--geometry = "coil_and_pan" -- more elements in the coarse grid, very angular corners, no projectors
--geometry = "coil_and_pan_v3" -- less elements in the coarse grid, better resolution of the corners, but no projectors
geometry = "coil_and_pan_proj" -- the same as _v3, but with the projectors rounding the corners

-- Remark:
-- The order of the subsets in the geometry influences the convergence of the
-- linear solver for the discretization of the Maxwell equation (because it
-- determines the ordering of the DoFs and the problem suffers jumps of the coefficients).
-- The order that works looks like
-- box, pan, coilPos, coilNeg, boxBnd, panBnd, cut0, cut1, coilBnd

-- Grid file name
gridName = "grids/" .. geometry .. ".ugx"

-- Subsets used in the problem
neededSubsets = {"box", "pan", "coilPos", "coilNeg", "cut0", "cut1", "boxBnd"}

print (" Geometry: " .. geometry .. " (file " .. gridName .. ")")
print (" Choosen Parameter:")
print ("    numRefs    = " .. numRefs)
print ("    numPreRefs = " .. numPreRefs)

-- Init UG for dimension and algebra
InitUG (dim, AlgebraType("CPU", 2)); -- note: the block size should be 2

--------------------------------------------------------------------------------
-- Domain Setup
--------------------------------------------------------------------------------

-- Create the domain, load the grid and refine it
dom = util.CreateDomain (gridName, numPreRefs, neededSubsets)
balancer.RefineAndRebalanceDomain (dom, numRefs - numPreRefs)

print ("Domain-info:")
print (dom:domain_info():to_string())

-- Save the geometry of the grid hierarchy
--SaveGridHierarchyTransformed (dom:grid (), dom:subset_handler (), "coil_and_pan_refined.ugx", 2.5)

-- Set the electromagnetic parameters to the subdomains
em = EMaterial (dom)
em:add ("box", 1.0, 0.0)
em:add ("coilPos,coilNeg", 1.0, 0.0)
em:add ("pan", 1.0, 1.0)
em:close ()

-- Frequency of the current
omega = 1

-- Create the edge-centered approximation space for E
print ("--> Edge-centered DoF distribution")
edgeApproxSpace = ApproximationSpace (dom)
edgeApproxSpace:add_fct ("r", "Nedelec", 1) -- Re
edgeApproxSpace:add_fct ("i", "Nedelec", 1) -- Im

edgeApproxSpace:init_levels ()
edgeApproxSpace:init_top_surface ()
edgeApproxSpace:print_statistic ()

-- Create the vertex-centered approximation space for the potential
print ("--> Vertex-centered DoF distribution")
vertApproxSpace = ApproximationSpace (dom, AlgebraType("CPU", 1))
vertApproxSpace:add_fct ("u", "Lagrange", 1)

vertApproxSpace:init_levels ()
vertApproxSpace:init_top_surface ()
vertApproxSpace:print_statistic ()

--------------------------------------------------------------------------------
-- FE Disc setup
--------------------------------------------------------------------------------

-- Local discretization
elemDisc = EddyCurrent_E_Nedelec ("r,i", em, omega) -- last argument = frequency omega

-- Dirichlet BC
dirichletBC = NedelecDirichletBC ("r,i")
dirichletBC:add_0 ("boxBnd")

-- Global discretization
domainDisc = DomainDiscretization (edgeApproxSpace)
domainDisc:add (elemDisc)
domainDisc:add (dirichletBC)

--------------------------------------------------------------------------------
--  Linear solver for the Nedelec-element-based discretization
--------------------------------------------------------------------------------

-- matrix and vectors
A = MatrixOperator ()
u = GridFunction (edgeApproxSpace)
b = GridFunction (edgeApproxSpace)
JG = GridFunction (edgeApproxSpace)

-- edge-centered smoother in the coarse grid solver
edgeBaseSmoother = GaussSeidel ()

-- vertex-centered smoother in the coarse grid solver (should be based on a scalar algebra)
vertBaseSmoother = GaussSeidelCPU1 ()

-- hybrid smoother in the coarse grid solver
baseHybridSmoother = HiptmairHybridSmoother (vertApproxSpace, edgeBaseSmoother, vertBaseSmoother)
baseHybridSmoother:set_Dirichlet (dirichletBC)

-- convergence check for the coarse solver
baseConvCheck = ConvCheck ()
baseConvCheck:set_maximum_steps (1024)
baseConvCheck:set_minimum_defect (1e-10)
baseConvCheck:set_reduction (1e-10)
baseConvCheck:set_verbose (false)

-- coarse grid solver
baseSolver = BiCGStab ()
baseSolver:set_preconditioner (baseHybridSmoother)
baseSolver:set_convergence_check (baseConvCheck)

-- edge-centered smoother
edgeSmoother = GaussSeidel ()
edgeSmoother:enable_consistent_interfaces (true)

-- vertex-centered smoother (should be based on a scalar algebra)
vertSmoother = GaussSeidelCPU1 ()
vertSmoother:enable_consistent_interfaces (true)

-- hybrid smoother
hybridSmoother = HiptmairHybridSmoother (vertApproxSpace, edgeSmoother, vertSmoother)
hybridSmoother:set_Dirichlet (dirichletBC)

-- transfer operators
transferOp = NedelecTransfer (edgeApproxSpace)

-- geometric multigrid method
gMG = GeometricMultiGrid (edgeApproxSpace)
gMG:set_discretization (domainDisc)
gMG:set_base_level (0)
gMG:set_base_solver (baseSolver)
gMG:set_gathered_base_solver_if_ambiguous (false)
gMG:set_smoother (hybridSmoother)
gMG:set_transfer (transferOp)
gMG:set_cycle_type (1)
gMG:set_num_presmooth (4)
gMG:set_num_postsmooth (4)

-- convergence check
ConvCheck = ConvCheck ()
ConvCheck:set_maximum_steps (1024)
ConvCheck:set_minimum_defect (1e-9)
ConvCheck:set_reduction (1e-10)
ConvCheck:set_verbose (true)

-- solver for the discretization
--linSolver = LinearSolver ()
linSolver = BiCGStab ()
linSolver:set_preconditioner (gMG)
--linSolver:set_preconditioner (hybridSmoother)
linSolver:set_convergence_check (ConvCheck)

--------------------------------------------------------------------------------
--  Projection
--------------------------------------------------------------------------------

-- vertex-centered smoother for the coarse grid solver of the projection
projBaseSmoother = ILUCPU1 ()
--projBaseSmoother:enable_consistent_interfaces (true)

-- convergence check for the coarse grid solver of the projection
projBaseConvCheck = ConvCheckCPU1 ()
projBaseConvCheck:set_maximum_steps (1024)
projBaseConvCheck:set_minimum_defect (1e-14)
projBaseConvCheck:set_reduction (1e-14)
projBaseConvCheck:set_verbose (false)

-- coarse grid solver for the projection
projBaseSolver = CGCPU1 ()
projBaseSolver:set_preconditioner (projBaseSmoother)
projBaseSolver:set_convergence_check (projBaseConvCheck)

-- vertex-centered smoother for the projection
projSmoother = ILUCPU1 ()
--projSmoother:enable_consistent_interfaces (true)

-- geometric multigrid method for the projection
projGMG = _G ["GeometricMultiGrid"..dim.."dCPU1"] (vertApproxSpace)
projGMG:set_base_level (0)
projGMG:set_base_solver (projBaseSolver)
projGMG:set_smoother (projSmoother)
projGMG:set_cycle_type (1)
projGMG:set_num_presmooth (1)
projGMG:set_num_postsmooth (1)

-- convergence check for the projection
projConvCheck = ConvCheckCPU1 ()
projConvCheck:set_maximum_steps (1024)
projConvCheck:set_minimum_defect (1e-14)
projConvCheck:set_reduction (1e-14)
projConvCheck:set_verbose (true)

-- linear solver for the projection
projSolver = CGCPU1 ()
projSolver:set_preconditioner (projGMG)
projSolver:set_convergence_check (projConvCheck)

-- projection
projection = NedelecProject (em, vertApproxSpace, projSolver)
projection:set_Dirichlet (dirichletBC)

--------------------------------------------------------------------------------
--  Generator current
--------------------------------------------------------------------------------

-- Loop source in the coil
gen_current = NedelecLoopCurrent ("coilNeg", "coilPos", "cut0", vertApproxSpace, projSolver)
gen_current:set ("r", 1.0)

print ("--> Computation of the generator current")

-- compute source
projGMG:add_prolongation_post_process (gen_current:zero_average ())
JG:set (0.0)
gen_current:compute (JG)
elemDisc:set_generator_current (JG, "r,i", gen_current:subsets ())
projGMG:clear_transfer_post_process ()

--------------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------------

print ("--> Solving the discretized problem")

-- 1. assemble matrix and rhs
domainDisc:assemble_linear (A, b)

-- 2. set dirichlet values and start iterate
u:set (0.0)
domainDisc:adjust_solution (u)

-- 3. init solver for linear Operator
linSolver:init (A, u)

-- 4. apply solver
linSolver:apply_return_defect (u, b)

--------------------------------------------------------------------------------
--  Apply projection
--------------------------------------------------------------------------------

print ("--> Projection of the solution")

projection:apply (u, "r,i")

--------------------------------------------------------------------------------
--  Output
--------------------------------------------------------------------------------

print ("--> Output")

-- compose the VTU-file names
grid_level = numPreRefs + numRefs
entire_vtu_file_name = "PanSolution3d-".. geometry .. "-frq" .. omega .. "-lev" .. grid_level;
coil_vtu_file_name = "PanSolution3d-".. geometry .. "-Coil-frq" .. omega .. "-lev" .. grid_level;
pan_vtu_file_name = "PanSolution3d-".. geometry .. "-Pan-frq" .. omega .. "-lev" .. grid_level;

-- electric field
ReEData = NedelecGridFunctionData (u, "r")
ImEData = NedelecGridFunctionData (u, "i")

ReJData = NedelecGridFunctionData (JG, "r")
ImJData = NedelecGridFunctionData (JG, "i")

-- magnetic induction
ReBData = EddyCurrentReBofEUserData (u, "r,i", omega)
ImBData = EddyCurrentImBofEUserData (u, "r,i", omega)

-- heat source
heat = EddyCurrentHeat (u, "r,i", em)

-- VTK-objects

print ("  > Entire domain")
entire_out = VTKOutput ()
--entire_out:set_binary (false)
entire_out:clear_selection ()
entire_out:select_element (ReEData, "ReE")
entire_out:select_element (ImEData, "ImE")
entire_out:select_element (ReBData, "ReB")
entire_out:select_element (ImBData, "ImB")
entire_out:select_element (heat, "HeatSrc")
entire_out:print (entire_vtu_file_name, u)

print ("  > The coil subdomain")
coil_out = VTKOutput ()
--coil_out:set_binary (false)
coil_out:clear_selection ()
coil_out:select_element (ReJData, "ReJ")
coil_out:select_element (ImJData, "ImJ")
coil_out:select_element (ReEData, "ReE")
coil_out:select_element (ImEData, "ImE")
coil_out:select_element (ReBData, "ReB")
coil_out:select_element (ImBData, "ImB")
coil_out:select_element (SubsetIndicatorUserData (dom, "coilPos,coilNeg,coilBnd,cut0,cut1"), "Coil")
coil_out:print_subsets (coil_vtu_file_name, u, "coilPos,coilNeg,coilBnd,cut0,cut1")

print ("  > The pan subdomain")
pan_out = VTKOutput ()
--pan_out:set_binary (false)
pan_out:clear_selection ()
pan_out:select_element (ReEData, "ReE")
pan_out:select_element (ImEData, "ImE")
pan_out:select_element (ReBData, "ReB")
pan_out:select_element (ImBData, "ImB")
pan_out:select_element (heat, "HeatSrc")
pan_out:select_element (SubsetIndicatorUserData (dom, "pan,panBnd"), "Pan")
pan_out:print_subsets (pan_vtu_file_name, u, "pan,panBnd")

--------------------------------------------------------------------------------
--  Done
--------------------------------------------------------------------------------

print ("--> Done")

-- End of File
