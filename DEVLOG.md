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
