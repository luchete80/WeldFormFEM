20250526 - Added Linear Matrices
         - Added Force Assembly
         - Computation of Deformation Gradient
20250529 - Fixed Bmat calculation
         - Added class Eigen Solver
         
<<<<<<< HEAD
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
=======
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

