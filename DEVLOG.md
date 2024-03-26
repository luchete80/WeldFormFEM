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
         