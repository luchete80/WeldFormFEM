// TODO: extend to all dirs
//NOTE: DENSITY IS OF ELEMENTS
//This also will be passed to device
#ifndef CUDA_BUILD
#include "double3_c.h"
#endif
#include "Mesh.h"

#define max(a,b) a>b?a:b;
#define min(a,b) a<b?a:b;

#define PRINT_V(v) printf("%f %f %f\n",v.x,v.y,v.z);
namespace MetFEM{
// ORIGINAL CPU version
// inline void TriMesh::Move(const double &dt){
	// //Seems to be More accurate to do this by node vel
	// //This is used by normals
  // Vec3_t min = 1000.;
  // Vec3_t max = -1000.;
	// for (int n=0;n<node.Size();n++){
    // Vec3_t vr 	= cross(m_w, *node[n]);
    // *node_v[n] = m_v + vr;
    // for (int i=0;i<3;i++) {
      // if      ((*node[n])(i) < min(i)) min[i] = (*node[n])(i);
      // else if ((*node[n])(i) > max(i)) max[i] = (*node[n])(i);
    // } 
		// *node[n] += (*node_v[n])*dt;
	// }
  
  // //printf( "Min Max Node pos" << min<< "; " <<max<<endl;
  
  // CalcCentroids();
  // CalcNormals();        //From node positions
  // UpdatePlaneCoeff();   //pplane
// }

/*
__global__ inline void MeshUpdateKernel(TriMesh_d *mesh_d, double dt) {
 	mesh_d->Move(dt);
  mesh_d->CalcCentroids();
  mesh_d->CalcNormals();
  mesh_d->UpdatePlaneCoeff(); 
}


*/
  TriMesh_d::TriMesh_d() : 
    node(nullptr), node_v(nullptr), elnode(nullptr),
    centroid(nullptr), normal(nullptr), pplane(nullptr),
    nfar(nullptr), 
    //nod_mesh_id(nullptr), 
    ele_mesh_id(nullptr),
    T(nullptr),
    nodecount(0), elemcount(0), mesh_count(0),
    current_node_capacity(0), current_elem_capacity(0),
    allocated_meshes(0) {
    m_v = m_w = make_double3(0.,0.,0.);
}

//NOW THIS IS ZORIENTED, CHANGE TO EVERY PLANE
void TriMesh_d::AxisPlaneMesh(const int &id, const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens){
	#ifndef CUDA_BUILD
	double x1,x2,x3;
	double l1,l2;
	double3 p = p2-p1;
	int dir[3];
	if 			(axis == 0 )	{dir[0] = 1; dir[1] = 2;}
	else if (axis == 1 )	{dir[0] = 0; dir[1] = 2;}
	else									{dir[0] = 0; dir[1] = 1;}
	
	dir [2] = axis; //dir2 is which remains constant
	
	x3 = p1.z;
	x2 = p1.y; 
  
  //TODO: CORRECT
  //x3 = p1.y;
  //x2 = p1.y;
  
	//double dl = p(dir[0])/dens;	//Could be allowed 2 diff densities
  double dl = p.x/dens;
  nodecount = (dens+1)*(dens+1);
  
  //Is it necessary to paralellize mesh nodes??
  //cudaMalloc((void **)&node   , 	nodecount * sizeof (double3));
  //cudaMalloc((void **)&node_v , 	nodecount * sizeof (double3));
  
  malloc_t (node,      double3,nodecount);
  malloc_t (node_v,    double3,nodecount);
  
  double3 *node_h, *node_vh;
  node_h  =  new double3 [nodecount];
  node_vh =  new double3 [nodecount];
  
	//printf("dens: "<<dens<<endl;
	//Plane is in 0 and 1 dirs
  printf("Creating nodes..");
	int vi=0;
	int test =dens+1;
	for (int j=0; j<test; j++) {
		//x1 = p1(dir[0]);
    x1 = p1.x;
		for (int i=0; i<test; i++){
			double3 v;
			v.x=x1;v.y=x2;v.z=x3;
			//printf( "i,j" << i << ", " << j<<endl; 
			//node.Push(new double3(x1,x2,x3));
			node_h[vi]		=make_double3(v.x,v.y,v.z);
			node_vh[vi]	=make_double3(0.,0.,0.);
      vi++;
			// node.Push(new double3(v(0),v(1),v(2)));
			// node_v.Push(new double3(0.,0.,0.));
			//printf( "xyz: "<<x1 << ", "<<x2<<", "<<x3<<endl;
			x1+=dl;
		}
		x2+=dl;
	}
  //cudaMemcpy(node, node_h,    nodecount * sizeof (double3), cudaMemcpyHostToDevice);
  //cudaMemcpy(node_v, node_vh, nodecount * sizeof (double3), cudaMemcpyHostToDevice);

  memcpy_t(node,       node_h,  nodecount * sizeof(double3));
  memcpy_t(node_v,    node_vh,  nodecount * sizeof(double3));
    

	int n[4];
	int el =0;
	int i;
	
	elemcount = dens * dens * 2;
  printf( "Element count: %d",elemcount );  
  printf( "done. Creating elements... ");
  malloc_t (centroid,      double3,elemcount);
  malloc_t (normal,        double3,elemcount);
  malloc_t (elnode,        int, 3 * elemcount);
  
  malloc_t (react_force,   	double3, 1);
  malloc_t (react_p_force,  double, 1);
	//cudaMalloc((void **)&centroid , 	elemcount * sizeof (double3));
	//cudaMalloc((void **)&normal 	, 	elemcount * sizeof (double3));
	//cudaMalloc((void **)&elnode 	, 	3 * elemcount * sizeof (int));	
  
  int *elnode_h = new int[3*elemcount];
  double3 *centroid_h = new double3[elemcount];
  double3 *normal_h   = new double3[elemcount];
	
  if (!node_h || !node_vh || !elnode_h || !centroid_h || !normal_h) {
    printf("Memory allocation failed!\n");
    exit(1);
  } else {
    printf("Allocation success.\n");
  }

	for (size_t j = 0 ;j  < dens; j++ ) {
				// printf("j, dens" <<j<<", "<<dens<<endl;
				// printf("j<dens"<< (j  < dens)<<endl;
		for ( i = 0; i < dens; i++ ){
				// printf("i, dens" <<i<<", "<<dens<<endl;
				// printf("i <dens"<< (i  < dens)<<endl;
				n[0] = (dens + 1)* j + i; 		n[1] = n[0] + 1; 
				n[2] = (dens + 1)* (j+1) + i; n[3] = n[2] + 1;
			//printf(" jj" << jj<<endl;
			int elcon[2][3];	// TODO: check x, y and z normals and node direction 
												// For all plane orientations
			//If connectivity  is anticlockwise normal is outwards
			if (positaxisorent) {
				elcon[0][0] = n[0];elcon[0][1] = n[1];elcon[0][2] = n[2];
				elcon[1][0] = n[1];elcon[1][1] = n[3];elcon[1][2] = n[2];
			} else {
				elcon[0][0] = n[0];elcon[0][1] = n[2];elcon[0][2] = n[1];
				elcon[1][0] = n[1];elcon[1][1] = n[2];elcon[1][2] = n[3];				
			}
			for ( int e= 0; e<2;e++) { // 2 triangles
				int elnodeid = 3*el;
				//element.Push(new Element(elcon[e][0],elcon[e][1],elcon[e][2]));		
				elnode_h[elnodeid + 0] = elcon[e][0]; 
				elnode_h[elnodeid + 1] = elcon[e][1]; 
				elnode_h[elnodeid + 2] = elcon[e][2];
				//printf( "Element "<< el <<": ";
				// for (int en = 0 ; en<3; en++) printf( elcon[e][en]<<", ";
				// printf(endl;
				
				double3 v = ( node_h[elcon[e][0]] + node_h[elcon[e][1]] + node_h[elcon[e][2]] ) / 3. ;
				//element[el] -> centroid = v; 
				centroid_h[el] = v;
				//printf( "Centroid element %d = %f %f %f\n" ,el,v.x,v.y,v.z);
				el++;
			}
		}// i for
		
	}

	///////////////////////////////////////////
	//// MESH GENERATION END
	printf( "Done. Creating normals\n");
  printf( "elem count %d\n",elemcount);
	for (int e = 0; e < elemcount; e++){ 
		double f=-1.;
    normal_h[e].x = normal_h[e].y = normal_h[e].z = 0.0;
		if (positaxisorent) f= 1.;
		//element[e] -> normal (axis) = f;
		if (axis == 0)			normal_h[e].x = f;
		else if (axis == 1)	normal_h[e].y = f;
		else 								normal_h[e].z = f;
    
    // if (length(normal_h[e])<1.0e-3) printf( "ERROR. ZERO NORMAL"<<endl;
    // if (normal_h[e].y > 1.0e-10) printf( "ERROR. NORMAL Y NOT ZERO"<<endl;    
    
    //printf( "normal_h[e] "<<normal_h[e].x << ", " << normal_h[e].y << ", " <<normal_h[e].z<<endl;
	}
  //printf("Done\n");
  m_v = m_w = make_double3(0.,0.,0.);
  
  
  malloc_t (pplane,      double,elemcount);  
  malloc_t (nfar,        int   ,elemcount);  

  malloc_t (ele_mesh_id,      int,elemcount);  

  malloc_t (mu_sta,      double,1);  
  malloc_t (mu_dyn,      double,1);
  //malloc_t (nod_mesh_id,      int,nodecount);  

  if (!pplane || !nfar ) {
    printf("Memory allocation failed!\n");
    exit(1);
  } else {
    printf("Allocation success.\n");
  }
  
  int *ele_mesh_id_h   = new int[elemcount];
  for (int n=0;n<elemcount;n++){ele_mesh_id_h[n]=id;}
    
  //cudaMalloc((void **)&pplane , 	elemcount * sizeof (double));
  //cudaMalloc((void **)&nfar   , 	elemcount * sizeof (int));
  
  //cudaMemcpy(elnode, elnode_h, 3 * elemcount * sizeof(int), cudaMemcpyHostToDevice);
  //cudaMemcpy(centroid, centroid_h, elemcount * sizeof(double3), cudaMemcpyHostToDevice);
  //cudaMemcpy(normal, normal_h, elemcount * sizeof(double3), cudaMemcpyHostToDevice);
  
  memcpy_t(elnode,    elnode_h, 3 * elemcount * sizeof(int));
  memcpy_t(centroid,  centroid_h,    elemcount * sizeof(double3));
  memcpy_t(normal,    normal_h,     elemcount * sizeof(double3));


  memcpy_t(ele_mesh_id,    ele_mesh_id_h,     elemcount * sizeof(int));  
  
  

  delete[] node_h; //WHY THIS CRASHES
  delete[] node_vh; //WHY THIS CRASHES
  
  printf("End mesh gen\n");
  //////CRASHES IN CPU
  /////--------------
  delete[] elnode_h;
  delete[] centroid_h;
  delete[] normal_h;  
  delete[] ele_mesh_id_h;  
  #endif
  mesh_count = 1;
}

///// TODO: REPLACE WITH memcpy_t, malloc_t

// //// RESIZING
void TriMesh_d::ResizeStorage(int new_mesh_capacity) {
    // Allocate with malloc_t
    double3 *new_react_force;
    double *new_react_p_force;
    // malloc_t(new_mesh_vel, double3, new_mesh_capacity);
    // malloc_t(new_mesh_rot, double3, new_mesh_capacity);
    printf("Newmesh capacity: %d\n", new_mesh_capacity);
    printf("NEW MESH CAPACIDY %d\n",new_mesh_capacity);
    malloc_t(new_react_force, 	double3, new_mesh_capacity);
    malloc_t(new_react_p_force, double, new_mesh_capacity);
	
    // Copy existing data
    if (mesh_count > 0) {
        //memcpy(new_node_offsets, node_offsets, (mesh_count + 1) * sizeof(int));
        //memcpy(new_elem_offsets, elem_offsets, (mesh_count + 1) * sizeof(int));
        memcpy(new_react_force, react_force, mesh_count * sizeof(double3));
        memcpy(new_react_p_force, react_p_force, mesh_count * sizeof(double));
        // memcpy(new_mesh_vel, mesh_velocities, mesh_count * sizeof(double3));
        // memcpy(new_mesh_rot, mesh_rotations, mesh_count * sizeof(double3));
    }

    // Initialize new slots
    // if (new_mesh_capacity > allocated_meshes) {
        // for (int i = mesh_count + 1; i <= new_mesh_capacity; i++) {
            // new_node_offsets[i] = new_node_offsets[mesh_count];
            // new_elem_offsets[i] = new_elem_offsets[mesh_count];
        // }
    // }

    // Cleanup old memory
    // free_t(node_offsets);
    // free_t(elem_offsets);
    free_t(react_force);
    free_t(react_p_force);
    // free_t(mesh_velocities);
    // free_t(mesh_rotations);

    // Update pointers
    // node_offsets = new_node_offsets;
    // elem_offsets = new_elem_offsets;
    react_force 	= new_react_force;
    react_p_force 	= new_react_p_force;
    // mesh_velocities = new_mesh_vel;
    // mesh_rotations = new_mesh_rot;

    allocated_meshes = new_mesh_capacity;
}


int TriMesh_d::ResizeNodeData(int new_capacity) {
    printf("New capacity %d\n", new_capacity);
    if (new_capacity <= nodecount) return false;

    // Allocate new arrays using malloc_t
    double3* new_nodes = nullptr;
    double3* new_node_v = nullptr;

    
    malloc_t(new_nodes, double3, new_capacity);
    malloc_t(new_node_v, double3, new_capacity);


    if (!new_nodes || !new_node_v 
    //|| !new_ele_mesh_id
    ) {
        printf("Memory allocation failed!\n");
        free_t(new_nodes);
        free_t(new_node_v);

        return false;
    }
    cout << "copying " <<nodecount<<" nodes "<<endl;
    if (node) {
        // Copy existing data
        memcpy_t(new_nodes, node, nodecount * sizeof(double3));
        memcpy_t(new_node_v, node_v, nodecount * sizeof(double3));
        
        // Free old memory
        free_t(node);
        free_t(node_v);
    }
    cout << "assigning "<<endl;
    // Update pointers
    node = new_nodes;
    node_v = new_node_v;
    //ele_mesh_id = new_ele_mesh_id;
    current_node_capacity = new_capacity;
    
    return true;
}

int TriMesh_d::ResizeElementData(int new_capacity) {
    if(new_capacity <= elemcount) return 0;
    
    // Allocate new arrays using malloc_t
    int* new_elnode = nullptr;
    double3* new_centroid = nullptr;
    double3* new_normal = nullptr;
    double* new_pplane = nullptr;
    int* new_nfar = nullptr;
    int* new_ele_mesh_id = nullptr;
        
    malloc_t(new_elnode, int, new_capacity * 3);
    malloc_t(new_centroid, double3, new_capacity);
    malloc_t(new_normal, double3, new_capacity);
    malloc_t(new_pplane, double, new_capacity);
    malloc_t(new_nfar, int, new_capacity);
    malloc_t(new_ele_mesh_id, int, new_capacity);

    if(!new_elnode || !new_centroid || !new_normal || !new_pplane || !new_nfar) {
        printf("Element memory allocation failed!\n");
        free_t(new_elnode);
        free_t(new_centroid);
        free_t(new_normal);
        free_t(new_pplane);
        free_t(new_nfar);
        free_t(new_ele_mesh_id);
        return 0;
    }

    if(elnode) {
        // Copy existing data
        memcpy_t(new_elnode, elnode, elemcount * 3 * sizeof(int));
        memcpy_t(new_centroid, centroid, elemcount * sizeof(double3));
        memcpy_t(new_normal, normal, elemcount * sizeof(double3));
        memcpy_t(new_pplane, pplane, elemcount * sizeof(double));
        memcpy_t(new_nfar, nfar, elemcount * sizeof(int));
        memcpy_t(new_ele_mesh_id, ele_mesh_id, elemcount * sizeof(int));
        
        // Free old memory
        free_t(elnode);
        free_t(centroid);
        free_t(normal);
        free_t(pplane);
        free_t(nfar);
        free_t(ele_mesh_id);
    }
    
    // Update pointers
    elnode = new_elnode;
    centroid = new_centroid;
    normal = new_normal;
    pplane = new_pplane;
    nfar = new_nfar;
    current_elem_capacity = new_capacity;
    ele_mesh_id= new_ele_mesh_id;
    
    return 1;
}

//AddMesh coordinates everything
void TriMesh_d::AddMesh(const TriMesh_d& new_mesh) {
    // 1. Validate input
    printf("Resizing storage\n");
    cout << "Allocated meshes: "<<allocated_meshes<<endl;
    if(mesh_count >= allocated_meshes) {
      cout << "assigning "<<endl;
        int new_capacity = allocated_meshes ? allocated_meshes * 2 : 4;
        cout << "done "<<endl;
        ResizeStorage(new_capacity);
    }
    printf("Resizing Node Data\n");
    // 2. Calculate new totals
    const int new_nodes = nodecount + new_mesh.nodecount;
    const int new_elems = elemcount + new_mesh.elemcount;

    printf("Old Node count %d, New node count %d, current_node_capacity %d\n", nodecount, new_nodes,current_node_capacity);
    //~ // 3. Resize data arrays with safety checks
    if(new_nodes > current_node_capacity) {
        int new_node_cap = max(new_nodes, current_node_capacity * 2);
        if(!ResizeNodeData(new_node_cap)) return;
    }else{
      printf("No allocated extra nodes.\n");
    }
    printf("Resizing Element Data\n");
    if(new_elems > current_elem_capacity) {
        int new_elem_cap = max(new_elems, current_elem_capacity * 2);
        if(!ResizeElementData(new_elem_cap)) return;
    }
    printf("Copying memory\n");
    // 4. Copy data with bounds checking
    if(node && new_mesh.node) {
        memcpy_t(node + nodecount, new_mesh.node, 
              new_mesh.nodecount * sizeof(double3));
    }
    if(node_v && new_mesh.node_v) {
        memcpy_t(node_v + nodecount, new_mesh.node_v,
              new_mesh.nodecount * sizeof(double3));
    }


    // Copy element data with node index offset
    if(new_mesh.elnode && new_mesh.elemcount > 0) {
        printf("Copying %d elements to offset %d\n", new_mesh.elemcount, elemcount);
        
        for(int e = 0; e < new_mesh.elemcount; e++) {
            // Calculate source and destination indices
            int src_idx = 3 * e;
            int dst_idx = 3 * (elemcount + e);
            
            // Copy and offset node indices
            elnode[dst_idx]     = new_mesh.elnode[src_idx]     + nodecount;
            elnode[dst_idx + 1] = new_mesh.elnode[src_idx + 1] + nodecount;
            elnode[dst_idx + 2] = new_mesh.elnode[src_idx + 2] + nodecount;
            
            // Copy other element data
            centroid[elemcount + e] = new_mesh.centroid[e];
            normal[elemcount + e] = new_mesh.normal[e];
            pplane[elemcount + e] = new_mesh.pplane[e];
            nfar[elemcount + e] = new_mesh.nfar[e];
        }
    }
    
    printf("Setting vel\n");
    //SETTING VEL
    // After copying data but before updating counters:
    int new_nodes_start = nodecount;
    int new_nodes_end = nodecount + new_mesh.nodecount;
    
    printf ("new_nodes_start %d, new_nodes_end %d\n",new_nodes_start, new_nodes_end);
    /////Set velocities only for the new nodes
    int no = 0;
    for (int n = new_nodes_start; n < new_nodes_end; n++) {
        node_v[n] = new_mesh.m_v;
        no++;
    }



    int new_elem_start = elemcount;
    int new_elem_end = elemcount + new_mesh.elemcount;
    
    printf ("new_elem_start %d, new_elem_end %d\n",new_nodes_start, new_nodes_end);
    /////Set velocities only for the new nodes
    no = 0;
    for (int n = new_elem_start; n < new_elem_end; n++) {
        ele_mesh_id[n] = new_mesh.ele_mesh_id[no];
        no++;
    }


    printf("Updating metadata..\n");
    // 5. Update metadata
    //node_offsets[mesh_count + 1] = new_nodes;
    //elem_offsets[mesh_count + 1] = new_elems;


    // 6. Update counters
    mesh_count++;
    nodecount = new_nodes;
    elemcount = new_elems;


}


  ///ONLY FOR NEW MESH ADDED
  void TriMesh_d::SetNodesVel(const double3 &v, int start_node, int end_node) {
      if (end_node == -1) end_node = nodecount; // default to all nodes
      end_node = min(end_node, nodecount); // ensure we don't go out of bounds
      
      m_v = v; // Still set the overall mesh velocity
      for (int n = start_node; n < end_node; n++) {
          node_v[n] = v;
      }
  }
  
  
//TODO: CHANGE TRIMESH NAME
TriMesh_d::TriMesh_d(NastranReader &nr, bool flipnormals){
  int dimension = nr.dim;
  //Insert nodes
  
  nodecount = nr.node_count;
  cout << "Nastran Mesh NodeCount: "<<nodecount<<endl;

  //Is it necessary to paralellize mesh nodes??
  //cudaMalloc((void **)&node   , 	nodecount * sizeof (double3));
  //cudaMalloc((void **)&node_v , 	nodecount * sizeof (double3));
  
  malloc_t (node,      double3,nodecount);
  malloc_t (node_v,    double3,nodecount);
  
  double3 *node_h, *node_vh;
  node_h  =  new double3 [nodecount];
  node_vh =  new double3 [nodecount];

  cout << "Copying "<< nodecount<<" nodes "<<endl;
  for (int n=0;n<nr.node_count;n++){

      //~ node.Push(new Vec3_t(nr.node[3*n],nr.node[3*n+1],nr.node[3*n+2]));
      node_h[n] = make_double3(nr.node[3*n],nr.node[3*n+1],nr.node[3*n+2]);
      cout << "Pos "<<nr.node[3*n]<<", "<<nr.node[3*n+1]<<", "<<nr.node[3*n+2]<<endl;

    node_vh[n] = make_double3(0.0,0.0,0.0);
		//~ node_v.Push(new Vec3_t(0.,0.,0.));
  }

  memcpy_t(node,       node_h,  nodecount * sizeof(double3));
  memcpy_t(node_v,    node_vh,  nodecount * sizeof(double3));
  cout << "endl"<<endl;
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //~ cout << "Generated "<<node.Size()<< " trimesh nodes. "<<endl;
  //~ //cout << "Normals"<<endl;
  elemcount = nr.elem_count;
  
  printf( "Element count: %d",elemcount );  
  printf( "done. Creating elements... ");
  
  malloc_t (centroid,      double3,elemcount);
  malloc_t (normal,        double3,elemcount);
  malloc_t (elnode,        int, 3 * elemcount);
  
  malloc_t (react_force,   	double3, 1);
  malloc_t (react_p_force,  double, 1);

  
  int *elnode_h = new int[3*elemcount];
  double3 *centroid_h = new double3[elemcount];
  double3 *normal_h   = new double3[elemcount];
	
  
  cout << "Writing elements..."<<endl;
  for (int e=0;e<nr.elem_count;e++){
    double cz = nr.elcon[3*e+2];
    if (dimension == 2) cz = 0;
    //~ if (!flipnormals)   element.Push(new Element(nr.elcon[3*e],nr.elcon[3*e+1],cz));		  
    //~ else                element.Push(new Element(nr.elcon[3*e+1],nr.elcon[3*e],cz));		        
    if (!flipnormals) {
      elnode_h[3*e  ] = nr.elcon[3*e];
      elnode_h[3*e+1] = nr.elcon[3*e+1];
      elnode_h[3*e+2] = cz;		  
    }else{   
      elnode_h[3*e  ] = nr.elcon[3*e+1];
      elnode_h[3*e+1] = nr.elcon[3*e  ];
      elnode_h[3*e+2] = cz;		                    
      //element.Push(new Element(nr.elcon[3*e+1],nr.elcon[3*e],cz));		  
    }
    //~ Vec3_t v;
		//~ if (dimension ==3) v = ( *node[nr.elcon[3*e]] + *node[nr.elcon[3*e+1]] + *node[nr.elcon[3*e+2]] ) / 3. ;
    //~ else               v = ( *node[nr.elcon[3*e]] + *node[nr.elcon[3*e+1]])  / 2. ;
    
    double3 v;
		//~ if (dimension ==3) v = ( *node[nr.elcon[3*e]] + *node[nr.elcon[3*e+1]] + *node[nr.elcon[3*e+2]] ) / 3. ;
    //~ else               v = ( *node[nr.elcon[3*e]] + *node[nr.elcon[3*e+1]])  / 2. ;

		if (dimension ==3) v = ( node_h[nr.elcon[3*e]] + node_h[nr.elcon[3*e+1]] + node_h[nr.elcon[3*e+2]] ) / 3. ;
    else               v = ( node_h[nr.elcon[3*e]] + node_h[nr.elcon[3*e+1]])  / 2. ;
    centroid_h[e] = v;
    
    //~ //TODO: CHANGE FOR CALCNORMALS
    if (dimension==3){
      //~ Vec3_t v1, v2;
      double3 v1,v2;      
      //~ //In COUNTERCLOCKWISE
      //~ v1 = *node[nr.elcon[3*e+1]] - *node[nr.elcon[3*e]];
      //~ v2 = *node[nr.elcon[3*e+2]] - *node[nr.elcon[3*e]];
      //~ element[e] ->normal = cross (v1,v2);
      //~ element[e] ->normal /= Norm(element[e] ->normal);
      v1 = node_h[nr.elcon[3*e+1]] - node_h[nr.elcon[3*e]];
      v2 = node_h[nr.elcon[3*e+2]] - node_h[nr.elcon[3*e]];
      normal_h[e] = cross (v1,v2);
      normal_h[e] = normal_h[e]/norm(normal_h[e]);
      
      cout << "normal "<<normal_h[e].x<<", "<<normal_h[e].y<<", "<<normal_h[e].z <<endl;
    } else { //See calc normals
        //~ Vec3_t u = *node [element[e]->node[1]] - *node [element[e]->node[0]];
        double3 u = node_h[nr.elcon[3*e+1]] - node_h[nr.elcon[3*e]];
        //~ v[0] = -u[1];
        //~ v[1] =  u[0];
        //~ v[2] =  0.0;
        v.x = -u.y;
        v.y =  u.x;
        v.z =  0.0;
        //~ element[e] -> normal = v/norm(v);
        normal_h[e] = v/norm(v);
    }
  }/////ELEMENT
  //~ cout << "Generated "<<element.Size()<< " trimesh elements. "<<endl;  
  
  m_v = m_w = make_double3(0.,0.,0.);
  
  
  malloc_t (pplane,      double,elemcount);  
  malloc_t (nfar,        int   ,elemcount);  
  malloc_t (ele_mesh_id,      int,elemcount);  
  malloc_t (mu_sta,      double,1);  
  malloc_t (mu_dyn,      double,1);

  if (!pplane || !nfar ) {
    printf("Memory allocation failed!\n");
    exit(1);
  } else {
    printf("Allocation success.\n");
  }
  
  int *ele_mesh_id_h   = new int[elemcount];
  for (int n=0;n<elemcount;n++){ele_mesh_id_h[n]=id;}

  
  memcpy_t(elnode,    elnode_h, 3 * elemcount * sizeof(int));
  memcpy_t(centroid,  centroid_h,    elemcount * sizeof(double3));
  memcpy_t(normal,    normal_h,     elemcount * sizeof(double3));
  memcpy_t(ele_mesh_id,    ele_mesh_id_h,     elemcount * sizeof(int));  
  
  mesh_count = 1;

  delete[] node_h; //WHY THIS CRASHES
  delete[] node_vh; //WHY THIS CRASHES
  delete[] elnode_h;
  delete[] centroid_h;
  delete[] normal_h;  
  delete[] ele_mesh_id_h;  
}
  

};
