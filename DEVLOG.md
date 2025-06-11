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
