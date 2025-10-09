F90 PREVIOUS VERSION
20230522 - Added Material Matrix to example
20230523 - Added nodal Boundary Conditions boolean flags
         - Added bc reinforcement implementation
         - Fixed full integration (only 2 points were active)
20230524 - Added element strain calcs
         - Fixed element dissasemble (pass global to elem) to calculate strains
         - Added steps for alternative solver (From Benson 1992).
         - Added fixity and accel bc 
20230529 - Fixed mass matrix assembly
20230530 - Begining to work with density calc, reducing indices
20230531 - Adding 3x3 3D matrix inversion subroutines for 3d elements
         - Working on 3D Domain Generation
         - Adding hourglass control
-----------------------------------------------------------------------
20230602 - Allocate internal forces
20230607 - Added vtu output file writer
20230608 - Corrected several things of Mass matrix redimension
         - Fixed hourglass forces calculation (was not initialized)
         - Summed hourglass correction
20230609 - Fixed density calc errors
20230608 - Added element shared by each node to average constant stresses and strains
20230612 - Fixed Element volume for full integration elements
         - Fixed error in strain rate calculation (indices of tensor diagonal were wrong)
20230613 - Fixed hourglass controls (corrected coefficients)
         - Added hourglass for 2D elements
         - Fixed 2D reduced integral volume calc
         - Fixed global ext force assembe 
         - Fixed element forces calculation error (in diagonal term was not incremented by node)
         - Fixed Sign of rotation rate
         - Fixed element strain rate calculation (derivative matrix indices were wrong)
         - Fixed nodal forces calc, relaced derivative matrix by the onw with determinant calculated
20230617 - Begining to add contact things
20230625 - Working on cylinder example, mesh reading
         - Working with several elements example (cylinder)
20230623 - Fixed domain allocation in plane mesh reading (CSV)
20230627 - Fixed fortran vtk element output 
         - Fixed addboxlength node arrangement
         - Added Verlet Solver
         - Added Equivalent stress calc
20230628 - Fixed ERROR in rot rate value
         - Fixed hourglass calculation (increment)
         - Fixed Hourglass forces magnitude
20230629 - Added sigma symm tensor output to vtk
-----------------------------------------------------------------------------------
20230702 - Added Contact Mesh
2023075  - Fixed Derived matrices, volume and jacobian calc for 2D full integration
202307XX - Added Chung Hulbert Integrator
20230725 - Added shock viscosity calculation
-----------------------------------------------------------------------------------
20230810 - Added Stress calculated by deformation (and not by def gradient)
         - Fixed calcs
         - Corrected BC for acceleration (is important IN Genrealized alpha method since u is calculated from a)
         - Corrected DEFORMATION GRADIENT CALCULATION (WAS USED THE TOTAL DISPLACEMENT)
         - Correct first step u_inc definition, now displacements after first step are ok.
20230829 - Finally working Verlet integration with strain rates (OpenRadioss, ABAQUS and LS-DYNA style).
         - This implies is not necessary to solve polar decomposition, neither F (Green Naghdi).
         - On the other side, results are not exactly as above, so is pending to check TS, and corrections.-
         - Hourglass and reduced integration is not working yet.
20230830 - Chung Hulbert with strain rate (without polar decomp) is working ok. BUT LEAPFROG AND VERLET ARE DIVERGING BY THE MOMENT.
         - Added Leapfrog, and Chung hulbert
-------------------------------------------------------------------------------------
20231108 - Fixed 2D plain strain stress decomposition calculation
         - Fixed output (With numbers like E-221)
20231113 - Added contact surfaces and things.
20231128 - Added parallel processing on jacobian calc
				 - Ommitted Shape matrix calculation (only derivs needed). 
				 - Prallelized deriv matrix calc
20231129 - Parallelized Stresses & Strains
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
20230207 - FIXED HOURGLASSING 
20230425 - Fixed BCs

20231128 - Added MemCopy From Host to Device, Nodes and Connectivity
----------------------------------------------------------------------------------
20231201 - Written Element matrices (created locally by the moment)
20231211 - Fixed temporally jacobian creation (not crtashing, defining matrices at the begingin)
         - Adding Element Strain calculations
20231227 - Completed derivatives for 3D
20231228 - Added strain computation, allocating element shape derivatives
20231229 - Changed velocity allocation
         - Working with strain rate
----------------------------------------------------------------------------------
----------------------------------------------------------------------------------
20240206 - Added to flat ptr
         - Fixed derivative saving
20240208 - Working with reduced integration
20240209 - Added BC assignment
         - Added Velocity and Disp Prediction
         - Added Volume Calculation
20240214 - Added Sigma Allocation
         - Added material allocation
         - Begining to add sigma calculations
20240222 - Making to be build both as CPU AND GPU 
         - Fixed Volume allocation
         - ADDED ASSEMBLY WITHOUT LOCKING
20240223 - Ended Node Elements generation
         - Added m_nodel_local in order to assemblying through nodal loop
20240224 - Begining to add Vec3D (from tiny and great DynELA) TODO: MAKING GPU/CPU MATH SUBMODULE (BASED ON DYNELA (CPU) AND (EVENTUALLY) TENSOR (GPUSPH)
         - DISCARD DOUBLE3 AND REPLACE IT WITH VECTORS/MATRICES
20240225 - Added memcpy
20240226 - Added Solver, Mechanical and Matrices files for CPU Solver
20240226 - Added Matrix Derivative dimensions and offset
         - Fixed another derivarive allocation
         - ADDED COPY CONSTRUCTOR TO MATRIX CLASS AND OPERATOR=
20240229 - Fixed some GPU building, matrix is working
         - Corrected Derivatives calculation
         ----------------------------------------
20240302 - Added velocity and accel corrections.
         - Added mass matrices
20240308 - Fixed Vol Calc (Gasuss weight)
         - Added mass matrices
         - Fixed nodal pos calculations in predictions
20240312 - Fixed Vector to Ptr allocation
20240315 - Fixed getVel fnc. Now strains are calculated ok.
         - Corrected velocity components calculation of strani rate
         - Corrected getDerivative function (several indices and offsets were wrong).
         - Added velocity corrections
20240315 - Added Shape matrix calculation (for element masses and lumped matrix)
20240320 - Corrected mass matrix dimension
         - Added AssemblyForces device call from Kernel.
         - Added mass matrix calculation (parallel and CPU)
20240321 - FIXED assembly (Two allocations were missed).
         - Changed Initial Density
         - Added first stress calc
20240326 - Corrected Shear Stress calc. NOW GIVES OK STRESSES AND PRESSURE 
           AT INTERMEDIATE CALCS!!
         - Fixed material props example. 
         - Fixed Forces calculation (Gauss integration Weight was lost, and never cleared)
         - DISPLACEMENTS ARE CORRECT FOR REDUCED INTEGRATION BRICK AFTER FIRST TIME STEP!
          (CHUNG HULBERT SOLVER)
         - Added ACCEL BOUNDARY CONDITIONS
20240327 - Added HOURGLASS FORCES! Attention: GPU Force Elem calc is wrong
-----------------------------------------------
20240403 - Fixed symm indices
20240405 - Added hourglass forces sum (GPU)
         - Found error in displacement correction.
         - Found error in acceleration (internal forces were inverted)
         - Fixed Strain and Rot rate calc
         - Added CS0 to hourglass calc
         
20240409 - Added device print funcions (vectors and tensors)
-----------------------------------------------
20240516 - Fixed for 3D some crashing thing.  
         - Seems not to be working still on CUDA
         - 2D Mesh generation is done
20240517 - Added Example 2D plain strain, FIRST Working with deriv calc.
         - Fixed 2D derivatives
20240520 - FInally fixed C++ version! Problem was in time integration
20240525 - Fixed initialization of Shared Elements by each node.
         - Also fixed shared element nodes in box mesh (zincrement was wrong)
20240527 - Fixing some bug related to several element domain.
         - Added getNodePos() to domain.
         - Fixed BCs of 4 element example.
         - Found a bug on strain rate calc
         - Begining to add vtk output
20240528 - Found errors in fast strain rate calac (determinant offset)
         - str rate were not equal to zero on each gp iter 
         - Fixed getVelem() function. 
20240529 - Fixed fast strain rate calc
         - Corrected 4 elem example .
         - Adding vtk writer
20240531 - Added a missing Syncronize after element mass matrix calc
----------------------------------------------------------------------
20240605 - Removed Arch type. Still with problems on GPU/CPU
20240606 - Begining to Add VTU export 
20240607 - Begining to add reading from LSDyna file.
20240610 - Added Beam example -Not working-
         - Added Plastic Strain Output for not hardening rule
20240627 - F90 Version, made area weighted by default
20240628 - Fixed axisymm example, still not working
----------------------------------------------------------------------
20240702 - F90 version - Fixed volumetric weight 
20240703 - F90 version. Fixed mass for axisymm area weighted
         - WORKING AXISYMMETRIC AREA WEIGHTED VERSION (F90).
20240724 - Added exporting V and F to displacements.
         - FIXED BOUNDARY CONDITIONS ON C++.
         - C++ 3D BEAM EXAMPLE IS WORKING
         - FIXED CUDA diagonal mass (added calculate density, vol and mass)
         - Still with error.
----------------------------------------------------------------------
20240906 - Fixed CPU parallel (and should be built as release)
----------------------------------------------------------------------
20241105 - FIXED 2D strain rate calc
         - Added Axisymmetric calculation (still not working)
20241107 - Converting vectors to 2D
20241122 -  FIXED 3D NO CONVENGECE! PROBLEMWAS NOT INITIALIZING V AND ACCEL
----------------------------------------------------------------------
20241220 - Fixed contact! mesh generation was wrong
20241221 - Added cross product for CPU version,calc normals
         - Added CPU functions dot, length and calcPlaneCoeff, updated mesh move.
20241222 - Fixed contact mesh nfar calculation, added setvel to example
         - Removed BC movements to beam contact example
         - Fixed beam contact example element normals (were inverted)
         - Allocated and calculated contforce
20241223 - Added VTK format (no html, to include contact meshes)
         - Fixed beam contact example 
20241225 - Fixed tetra jacobian calc
20241226 - Fixed connectivity of tetra, and shape function calc
         - Fixed also global derivative matrix
20241227 - Write cell values (stresses), failing passing to point data.
20241228 - Added plastic strain increment equations
20241230 - Fixed 
----------------------------------------------------------------------
20250102 - Begining to add node element after create elements. 
         - Begining to add json reader
         - Adding reading Domain block as .k file
20250104 - Fixed mem things. Domain was allocated twice.
         - FIXED ALLOCATION WIRH SET DIMENSION 
         - Fixed node location device allocation 
20250106 - Fixed VTK output, rewriting elnod_h on each remeshing thing.
20250107 - Fixed yield stress initialization (was up until element count)
         - Added plastic increment (only for bilinear)
         - Added Cylidner example (still missing bc applying from json)
20250116 - Fixed VTK wiriting, not copying to elnod_h. Remains remesh
         - Fixing search external nodes for tetras
         - Fixed vtk writer formatting
         - Added sim time and output time output reading from json
20250117 - Begining to add contact to input reading
         - Add destructor to cpu vector_t, still freezing
         - Fixed vtk writing for contact problems
         - Corrected dtout to const, addedoutput to contact beam
         - Added Formula CalcPressure Average Noodal Pressure (ANP).
20250129 - Added omega_h library. Adapted some signals to MINGW32
         - Fixing macro things for CUDA & CPU MSVC
         - Change __inline__ to inline since is not vendor specifig
         - Change "ssize_t" to "size_t" of omega_h/tpl/pss
         - Now omega_h compiling ok for MSVC         
20250131 - Fixing contact algorithm (wrong elnode index)
         - Add contact algorithm without preduction byut with current node position
-------------------------------------------------------------------------------
20250203 - Added ReMesher example 
20250204 - Adding Function class
20250205 - Fixed tetrahedral masses 
         - Fixed Friction Contact Forces!!
         - Changing GPU broken things.
20250506 - Adapted both tetras & hexas via pressure and mass calculation 
         - Fixed Huge mem leaking in CPU file  
         - Fixed some commented hourglass for hexa
20250510 - Fixed mass calculation from vol
         - FIXED CHUNG HULBERT ACCEL CORRECTION! NOW C++ version WORKING 
         - Cube results is exactly like F90 ver.
         - With change in mass friction is working
20250213 - FINALLY FOUND ERROR ON HEAX BEAM! Problem is Nonlock (nodal) assembly.
20250214 - The error was founf in fact in HMOD
         - CFL factor is now readed from json.
         - Found error in read some .k tetra meshes 
20250217 - Adding script to run solver 
20250218 - Changing Matrix to all device in CUDA version
         - Corrected openmp, still not tested
         - Begining to add new parallel format

20250219 - Fixed several paralllel things. GPU RESULTS ARE GIVING OK.
         - Fixed & simplified prediction and correction displacements, vels, and accels for cuda.
20250220 - Fixed several element displacements. Working ok BUT NOT FOR ALL THE PROBLEM.
         - Added openmp cores to json
20250221 - Added headers to Sources.
         - Added Clock Wall time.
--------------------------------------------------------------------------------------
20250311 - Added mesh refinement built in examples amr_test2.cpp
         - BEGINING WITH ADAPT BUT STILL NEEDS TO CORRECT METRIC
         - cylinder_adapt_test.cpp
         - r3d_test.cpp to pass mesh
         - FIXED PARALEL LOOP PRAGMA IN MSVC
20250312 - ADDED SEVERAL MESH ADAPTING METHODS (Warp and now angle and length)
         - Fixed compiling for mingw
20250326 - Added 2d test. Runs ok but fails inside WF code.
         - Fixed Omega_h on MSVC.
         - Added example flag.
---------------------------------------------------------------------------------------
20250401 - Added Domain::Free() for remesh old fields (mem leak)
         - Mapped materials, and shear stresses and strain rates.
         - REMESH IN Solver type works
         - FIXED REMESHGING: ELEMENT SAVED vars
         - FIXED: Calculate cacobian BEFORE Vol and dens and nodal mass.
         - REMESHING SEEMS TO BE WORKING
20250404 - Added mmg.
         - Added selective angle. 
         - FIXED Initial and Volumes, Volumes. Now are calculating via detF. 
         - FIXED mmg Remeshing generation. Working but still have to pass nodes and element vals. 
20250409 - FIXED Nodal Volume and Mass OUTSIDE Remesher
20250411 - FIXED REMAP OF Plain Strain, and Sig yield
         - ADDED REMAP BCs
         - FIXED MMG remeshing errors of bad indices.
---------------------------------------------------------------
20250514 - Removing some remeshing things. 
         - Begining to test parallel.
         - Checking ok several processors.-
         - Restorign CUDA.
20250523 - FOUND ERROR IN REMESHING CRASH: 
         - MUST CHANGE VELOCITY MAPPING FOLLOWING MOMENTUM.
         - CHECK MASS CONSERVATION.
20250523 - Found memory leaking (errors in Domain->Free()).
20250526 - Added Linear Matrices
         - Added Force Assembly
         - Computation of Deformation Gradient
20250529 - Fixed Bmat calculation
         - Added class Eigen Solver
         
20250530 - Adding EIGEN Elastic Simple Solver 
         - Fixing getElementNode()
         - Adding PETSC first functions.
         - Adding Diricket Eigen BCs
         - SINCE SOLVING ELASTIC COMPRESSION EXAMPLE SYSTEM DIVERGES:
         - Adding hourglass forces for tetras (#5)
         - FIXED Eigen solver for column major
         - Setting some defaults
         - Fixed Material matrix. 
         - Validated Single Steel Element with Python and EIGEN.
---------------------------------------------------------------------------------------------------------------------
20250601 - Fixed nonlinear material code.
         - Adding single element UL mat+geo Nonlinearities python example.
20250602 - Added MATRIX by tensor mult matrices.
         - Changed almasi to Green Lagrange 
20250605 - Fixed SLIDING CONTACT!
         - Using accum displacement algorithm.
20250606 - Added Contact Forces
20250609 - Fixed Tetra Element Reading from LSDyna!! Was reading hexas.
         - Fixed stress output (were wrongly calculated)
         - Added Hollomon reading.
20250610 - Added EA/l Geometric Contact stiffness instead of 
         - Added damping. 
20250612 - Added Contact Forces by integrated Pressures instead of Forces
         - Added velocity reading (constant)
20250613 - Added Contact Force output   
         - Fixed Critical Contact Force calculation.
         - Fixed BUG in nfar calc.
20250617 - Fixed bug on m_mesh_in_contact (was nodal, but is elemental)
         - Added element area (but not updated) 
20250618 - Added friction coeff reading from file. 
         - Restored penalty forces calc and output. 
         - Added Total Contact Force Pressure Nodal and Elemental 
         - Fixed Element Area calc!
20250619 - Fied Element area in outer elements. Assumed 1 areaper element.
20250625 - Fixed Contact Conduction (was doubled multiplied by dt)
         - Added Plastic strain heating equations 
20250627 - Fixed crash at external face & node search in Remeshing
         - Problem in element size. 
20250629 - Set remesh length from domain instead of auto.
20250630 - Make thread safe nodal mass summation.
         - Make log with commit and date. 
         - Added New Pressure Calculation with hourglass.
         - FIXED CRITICAL pressure. 
         - No more need of ANP.
--------------------------------------------------------------------------------------------------------------------------
20250701 - Added several pressure calc algorithms.
         - Added laplacian filter.
20250705 - SET VELOCITY AND ACCEL TO 0 at remesh. First working.
20250707 - Added contact force reset to remesh
         - FIXED NODAL AREA CALCULATION!! (Problem with .other_face, restored)
         - Added Implicit main loop.
20250711 - Fixed Nodal mass.
         - Fixed BC based problems (crashing without contact).
         - Fixed Transpositions on implicit solver.
20250714 - Added XYZ symmetry.
20250714 - Fixed small accelerarion export.
20250716 - FIXED Force Calc.
         - Added GMT Material
         - Set To Zero Stab Params.
         - Fixed critic cp calc (now is dependent on material but constant)
20250717 - Fixed Contact Area Calculation (Common)
20250718 - Added Some remeshing improvements and settings. 
         - IMPLICIT SOLVER
         - Fixed element volume mult on K
20250723 - Added Implicit ficticious velocities assembly
         - Working on NR Solver
         - Fixed Dirichlet (symmetry and simple problems) application         
         - Fixed Main NR solver loop vars.
20250725 - FIXED Element Forces (and hence contact forces) in 3D tetras.
         - Fixed Kmat to be affected by 1/detJ (derivatives are in fact BxdetJ)
         - Corrected also  Kgeo
         - Reset K and r after each NR iteration.
20250728 - Fixed assembly (wrong on residuals)
         - Added selective damping based on residuals.
         - Working on relative tolerances.
         - Added diagonal regularization
         - Fixed sigma entirely
20250728 - Fixed derivatives in Fcalculation (were affected by detJ)
         - Fixed Element Kgeo dimensions.
----------------------------------------------------------------------------------------------------
20250805 - Fixed axisymmetric vol weight calcs.
20250805 - Working on AxiSymmetric, added tri/tetra domain.
         - Corrected Axisymm weight in strain and forces calc to tri/tetras
20250812 - Rolling back to previous Pressure with J total calculation.
         - Working on remesh: set interval FROM LAST REMESH! not begining
         - Fix: Calc Faces Areas when remesh
20250818 - FIXED LAMBA TOLS TO TOLERANCES TO 1.0e-10
         - ADDED NODAL FROM ELEM VELOC
20250820 - Fixed AXISYMM Mass weighted eqns.
20250821 - Added Boundary Conditions Zone reading!!
20250826 - Begining to work with Contact to work also on 2D  
20250827 - Fixed Hourglass to woirk with AxiSymm
------------------------------------------------------------------------------------------------------
20250908 - Working on loading meshes from files
         - Fixed mesh, generated with msh = new TriMesh() instead of copy constructor
         - Fixed Assign mesh domain mesh
20250909 - Fixed Normals (orientation and nodes were both inverted with flipnormals)
         - Working Contact with external mesh.
         - Wirking with implicit solver. Now reads bcs (Problems w/o contact).
20250916 - Added tangent Matrix for Hollomon.
         - Solved first elastic large problem with implicit.   
20250917 - Fixed crashing when parallel on implicit version.
20250918 - Fixed Nonzero Dirichlet BC on Eigen Solver. 
         - Corrected Dirichlet BC, being the difference between applied and 
         - FIRST IMPLICIT ELASTIC WORKING.
-------------------------------------------------------------------------------------------------------
20251009  - Fixed Kinetic Energy
          - Fixed m_fe(which should be zero and is not )