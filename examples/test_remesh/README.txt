Current commit:
2c36ae0..5afeca0

Sometime it breaks here


  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   MODULE MMG3D: 5.8.0 (Oct. 30, 2024)
  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  -- MMG3DLIB: INPUT DATA
  --  INPUT DATA COMPLETED.     0.000s

  -- PHASE 1 : ANALYSIS

  -- MESH QUALITY  15685
     BEST   0.998334  AVRG.   0.807578  WRST.   0.335050 (5446)
  -- PHASE 1 COMPLETED.     0.000s

  -- PHASE 2 : ISOTROPIC MESHING

  -- GRADATION : 1.300000 (2.300000)
                               0 splitted,        1 collapsed,        0 swapped, 5 iter.
                                                                      0 swapped,      202 moved, 3 iter.
          631 filtered,       12 splitted,        9 collapsed,       51 swapped,       32 moved, 6 iter.
                                                                      0 swapped,       38 moved, 4 iter.
  -- PHASE 2 COMPLETED.     1.000s

  -- MESH QUALITY  15703
     BEST   0.998334  AVRG.   0.810886  WRST.   0.404064 (12402)

  -- MESH PACKED UP
     NUMBER OF VERTICES       3163   CORNERS        0
     NUMBER OF TETRAHEDRA    15703
     NUMBER OF EDGES            80   RIDGES        80
     NUMBER OF TRIANGLES      2028

   MMG3DLIB: ELAPSED TIME  1.000s

  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   END OF MODULE MMG3D
  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

New node count 3163
Done.
OVERALL tetrahedron count 15703
NEW MESH. Done mapping
Node count 3161, ELEM COUNT 15685
MESH CREATED.
WRITING DOMAIN 3163 NODES 15703ELEMS
MAP NODAL
MAP NODAL VECTOR RAW (MMMG)
FOUND new node 0 For node 0
DONE
New Volume: 1.51286e-05
Element 0 barycenter INSIDE old element 0 with barycentric coords: (0.25, 0.25, 0.25, 0.25)
Element 1 barycenter INSIDE old element 1 with barycentric coords: (0.253012, 0.285977, 0.259013, 0.201999)
Element 2 barycenter INSIDE old element 2 with barycentric coords: (0.25, 0.25, 0.25, 0.25)
Element 3 barycenter INSIDE old element 3 with barycentric coords: (0.25, 0.25, 0.25, 0.25)
Element 4 barycenter INSIDE old element 4 with barycentric coords: (0.25, 0.25, 0.25, 0.25)
foundelem closrsrt
map pressure
Map done
MAP NODAL VECTOR RAW (MMMG)
FOUND new node 0 For node 0
creating domain
Setting explicit dimensions
setting deriv allocation size of 188436
mdim3 4 15703 1
SETTING DENSITY TO 27000
G 2.65e+10
copying
Setting conn
Allocating 15703 and 4nodes x elem
COPYING 62812 element nodes
Done
Size of Nodal shared Elements vector 62812
Node Elemets set.
Calculating masses
recovering velocities
deleting
MESH CHANGED
done. Face count: 32420
Total External Nodal Area: 3.3640e-03
Writing VTK output out_remesh_7999.vtk
Writing VTK output out_remesh_smooth7999.vtk
s warmup: 0
New dt: 3.92677e-10
Max vel 5.53514
ERROR, 40 nodes with acceleration too large
Writing VTK output out_remesh_after1_7999.vtk
Step 8000, Time 0.001462, End Time: 1.0000e-01, Step Time 3.9268e-10
Overall Time391.901 seconds
ERROR: DT 1.97083e-07
s warmup: 0.05
New dt: 1.57071e-09
Max vel 5.53514
ERROR, 40 nodes with acceleration too large
ERROR: DT 1.97082e-07
s warmup: 0.1
New dt: 0
Max vel inf
ERROR, 306 nodes with acceleration too large
ERROR: DT 1.97082e-07
s warmup: 0.15
New dt: 6.28284e-09
Max vel 5.53484
ERROR: DT 2.15579e-07
s warmup: 0.2
New dt: 1.07488e-08
Max vel 4.08016
s warmup: 0.25
New dt: nan
Max vel 0.709606
MAX DISP 0 0 0
Writing output
Writing VTK output out.vtk
Done.
Writing VTK output out_remesh.vtk
Overall elapsed time: 394.786 seconds
Program ended.


ANOTHER (COMMIT)
240ace8..67e949f
ONLY WITH MINOR CHANGES:

  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   MODULE MMG3D: 5.8.0 (Oct. 30, 2024)
  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  -- MMG3DLIB: INPUT DATA
  --  INPUT DATA COMPLETED.     0.000s

  -- PHASE 1 : ANALYSIS

  -- MESH QUALITY  15685
     BEST   0.998334  AVRG.   0.807578  WRST.   0.335050 (5446)
  -- PHASE 1 COMPLETED.     0.000s

  -- PHASE 2 : ISOTROPIC MESHING

  -- GRADATION : 1.300000 (2.300000)
                               0 splitted,        1 collapsed,        0 swapped, 5 iter.
                                                                      0 swapped,      202 moved, 3 iter.
          631 filtered,       12 splitted,        9 collapsed,       51 swapped,       32 moved, 6 iter.
                                                                      0 swapped,       38 moved, 4 iter.
  -- PHASE 2 COMPLETED.     1.000s

  -- MESH QUALITY  15703
     BEST   0.998334  AVRG.   0.810886  WRST.   0.404064 (12402)

  -- MESH PACKED UP
     NUMBER OF VERTICES       3163   CORNERS        0
     NUMBER OF TETRAHEDRA    15703
     NUMBER OF EDGES            80   RIDGES        80
     NUMBER OF TRIANGLES      2028

   MMG3DLIB: ELAPSED TIME  1.000s

  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   END OF MODULE MMG3D
  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

New node count 3163
Done.
OVERALL tetrahedron count 15703
NEW MESH. Done mapping
Node count 3161, ELEM COUNT 15685
MESH CREATED.
WRITING DOMAIN 3163 NODES 15703ELEMS
MAP NODAL
MAP NODAL VECTOR RAW (MMMG)
FOUND new node 0 For node 0
DONE
New Volume: 1.51286e-05
Element 0 barycenter INSIDE old element 0 with barycentric coords: (0.25, 0.25, 0.25, 0.25)
Element 1 barycenter INSIDE old element 1 with barycentric coords: (0.253012, 0.285977, 0.259013, 0.201999)
Element 2 barycenter INSIDE old element 2 with barycentric coords: (0.25, 0.25, 0.25, 0.25)
Element 3 barycenter INSIDE old element 3 with barycentric coords: (0.25, 0.25, 0.25, 0.25)
Element 4 barycenter INSIDE old element 4 with barycentric coords: (0.25, 0.25, 0.25, 0.25)
foundelem closrsrt
map pressure
Map done
MAP NODAL VECTOR RAW (MMMG)
FOUND new node 0 For node 0
creating domain
Setting explicit dimensions
setting deriv allocation size of 188436
mdim3 4 15703 1
SETTING DENSITY TO 27000
G 2.65e+10
copying
Setting conn
Allocating 15703 and 4nodes x elem
COPYING 62812 element nodes
Done
Size of Nodal shared Elements vector 62812
Node Elemets set.
Calculating masses
recovering velocities
deleting
MESH CHANGED
done. Face count: 32420
Total External Nodal Area: 3.3640e-03
Writing VTK output out_remesh_7999.vtk
Writing VTK output out_remesh_smooth7999.vtk
s warmup: 0
New dt: 3.92677e-10
Max vel 5.53514
ERROR, 40 nodes with acceleration too large
Writing VTK output out_remesh_after1_7999.vtk
Step 8000, Time 0.001462, End Time: 1.0000e-01, Step Time 3.9268e-10
Overall Time391.901 seconds
ERROR: DT 1.97083e-07
s warmup: 0.05
New dt: 1.57071e-09
Max vel 5.53514
ERROR, 40 nodes with acceleration too large
ERROR: DT 1.97082e-07
s warmup: 0.1
New dt: 0
Max vel inf
ERROR, 306 nodes with acceleration too large
ERROR: DT 1.97082e-07
s warmup: 0.15
New dt: 6.28284e-09
Max vel 5.53484
ERROR: DT 2.15579e-07
s warmup: 0.2
New dt: 1.07488e-08
Max vel 4.08016
s warmup: 0.25
New dt: nan
Max vel 0.709606
MAX DISP 0 0 0
Writing output
Writing VTK output out.vtk
Done.
Writing VTK output out_remesh.vtk
Overall elapsed time: 394.786 seconds
Program ended.

--------------------- THIS RUNS OK 

Writing VTK output out_remesh_7999.vtk
Writing VTK output out_remesh_smooth7999.vtk
s warmup: 0
New dt: 3.92677e-10
Max vel 5.53514
Writing VTK output out_remesh_after1_7999.vtk
Step 8000, Time 0.001462, End Time: 1.0000e-01, Step Time 3.9268e-10
Overall Time181.131 seconds
ERROR: DT 1.97083e-07
s warmup: 0.05
New dt: 1.57071e-09
Max vel 5.53514
ERROR: DT 1.97082e-07
s warmup: 0.1
New dt: 3.5341e-09
Max vel 5.53484
ERROR: DT 1.9708e-07
s warmup: 0.15
New dt: 6.28286e-09
Max vel 5.53321
ERROR: DT 1.97074e-07
s warmup: 0.2
New dt: 9.81702e-09
Max vel 5.52835
ERROR: DT 1.97066e-07
s warmup: 0.25
New dt: 1.41367e-08
Max vel 5.51771
ERROR: DT 1.97053e-07
s warmup: 0.3
New dt: 1.92419e-08
Max vel 5.49819
ERROR: DT 1.97036e-07
s warmup: 0.35
New dt: 2.51331e-08
Max vel 5.46631
s warmup: 0.4
New dt: 3.18105e-08
Max vel 5.41832
s warmup: 0.45
New dt: 3.92747e-08
Max vel 5.35054
s warmup: 0.5
New dt: 4.75264e-08
Max vel 5.25967
s warmup: 0.55
New dt: 5.65666e-08
Max vel 5.14344
s warmup: 0.6
New dt: 6.63967e-08
Max vel 5.00133
s warmup: 0.65
New dt: 7.70182e-08
Max vel 4.83544
s warmup: 0.7
New dt: 8.84332e-08
Max vel 4.65107
s warmup: 0.75
New dt: 1.00644e-07
Max vel 4.45593
s warmup: 0.8
New dt: 1.13653e-07
Max vel 4.25662
s warmup: 0.85
New dt: 1.27465e-07
Max vel 4.0509
s warmup: 0.9
New dt: 1.42086e-07
Max vel 3.81734
s warmup: 0.95
New dt: 1.57464e-07
Max vel 4.07734
Step 8100, Time 0.001480, End Time: 1.0000e-01, Step Time 1.9730e-07
Overall Time184.156 seconds
Step 8200, Time 0.001499, End Time: 1.0000e-01, Step Time 1.9742e-07