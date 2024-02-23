// TODO: extend to all dirs
//NOTE: DENSITY IS OF ELEMENTS
//This also will be passed to device
#include "Mesh.cuh"

#define PRINT_V(v) printf("%f %f %f\n",v.x,v.y,v.z);
namespace SPH{
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
  
  // //cout << "Min Max Node pos" << min<< "; " <<max<<endl;
  
  // CalcCentroids();
  // CalcNormals();        //From node positions
  // UpdatePlaneCoeff();   //pplane
// }

__global__ inline void MeshUpdateKernel(TriMesh_d *mesh_d, double dt) {
 	mesh_d->Move(dt);
  mesh_d->CalcCentroids();
  mesh_d->CalcNormals();
  mesh_d->UpdatePlaneCoeff(); 
}

//NOW THIS IS ZORIENTED, CHANGE TO EVERY PLANE
inline void TriMesh_d::AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens){
	
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
  cudaMalloc((void **)&node   , 	nodecount * sizeof (double3));
  cudaMalloc((void **)&node_v , 	nodecount * sizeof (double3));
  
  double3 *node_h, *node_vh;
  node_h  =  new double3 [nodecount];
  node_vh =  new double3 [nodecount];
  
	//cout <<"dens: "<<dens<<endl;
	//Plane is in 0 and 1 dirs
  cout << "Creating nodes.."<<endl;
	int vi=0;
	int test =dens+1;
	for (int j=0; j<test; j++) {
		//x1 = p1(dir[0]);
    x1 = p1.x;
		for (int i=0; i<test; i++){
			double3 v;
			v.x=x1;v.y=x2;v.z=x3;
			//cout << "i,j" << i << ", " << j<<endl; 
			//node.Push(new double3(x1,x2,x3));
			node_h[vi]		=make_double3(v.x,v.y,v.z);
			node_vh[vi]	=make_double3(0.,0.,0.);
      vi++;
			// node.Push(new double3(v(0),v(1),v(2)));
			// node_v.Push(new double3(0.,0.,0.));
			//cout << "xyz: "<<x1 << ", "<<x2<<", "<<x3<<endl;
			x1+=dl;
		}
		x2+=dl;
	}
  cudaMemcpy(node, node_h,    nodecount * sizeof (double3), cudaMemcpyHostToDevice);
  cudaMemcpy(node_v, node_vh, nodecount * sizeof (double3), cudaMemcpyHostToDevice);

  cout << "Element count: "<<elemcount << endl;  
  cout << "done. Creating elements... ";
	int n[4];
	int el =0;
	int i;
	
	elemcount = dens * dens * 2;
	cudaMalloc((void **)&centroid , 	elemcount * sizeof (double3));
	cudaMalloc((void **)&normal 	, 	elemcount * sizeof (double3));
	cudaMalloc((void **)&elnode 	, 	3 * elemcount * sizeof (int));	
  int *elnode_h = new int[3*elemcount];
  double3 *centroid_h = new double3[elemcount];
  double3 *normal_h   = new double3[elemcount];
	
	for (size_t j = 0 ;j  < dens; j++ ) {
				// cout <<"j, dens" <<j<<", "<<dens<<endl;
				// cout <<"j<dens"<< (j  < dens)<<endl;
		for ( i = 0; i < dens; i++ ){
				// cout <<"i, dens" <<i<<", "<<dens<<endl;
				// cout <<"i <dens"<< (i  < dens)<<endl;
				n[0] = (dens + 1)* j + i; 		n[1] = n[0] + 1; 
				n[2] = (dens + 1)* (j+1) + i; n[3] = n[2] + 1;
			//cout <<" jj" << jj<<endl;
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
				//cout << "Element "<< el <<": ";
				// for (int en = 0 ; en<3; en++) cout << elcon[e][en]<<", ";
				// cout <<endl;
				
				double3 v = ( node_h[elcon[e][0]] + node_h[elcon[e][1]] + node_h[elcon[e][2]] ) / 3. ;
				//element[el] -> centroid = v; 
				centroid_h[el] = v;
				//cout << "Centroid" << element[el] -> centroid << endl;
				el++;
			}
		}// i for
		
	}

	///////////////////////////////////////////
	//// MESH GENERATION END
	cout << endl<<"Done. Creating normals"<<endl;
  cout << "elem count "<<elemcount<<endl;
	for (int e = 0; e < elemcount; e++){ 
		double f=-1.;
    normal_h[e].x = normal_h[e].y = normal_h[e].z = 0.0;
		if (positaxisorent) f= 1.;
		//element[e] -> normal (axis) = f;
		if (axis == 0)			normal_h[e].x = f;
		else if (axis == 1)	normal_h[e].y = f;
		else 								normal_h[e].z = f;
    
    // if (length(normal_h[e])<1.0e-3) cout << "ERROR. ZERO NORMAL"<<endl;
    // if (normal_h[e].y > 1.0e-10) cout << "ERROR. NORMAL Y NOT ZERO"<<endl;    
    
    //cout << "normal_h[e] "<<normal_h[e].x << ", " << normal_h[e].y << ", " <<normal_h[e].z<<endl;
	}
  
  m_v = m_w = make_double3(0.,0.,0.);
  
  cudaMalloc((void **)&pplane , 	elemcount * sizeof (double));
  cudaMalloc((void **)&nfar   , 	elemcount * sizeof (int));
  
  cudaMemcpy(elnode, elnode_h, 3 * elemcount * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(centroid, centroid_h, elemcount * sizeof(double3), cudaMemcpyHostToDevice);
  cudaMemcpy(normal, normal_h, elemcount * sizeof(double3), cudaMemcpyHostToDevice);

  delete node_h;
  delete node_vh;
  delete elnode_h;
  delete centroid_h;
  delete normal_h;  
}

//This is done once, Since mesh is rigid
//Calculate radius and plane coefficient
inline __device__ void TriMesh_d::CalcSpheres(){
	// double max;
  int e = threadIdx.x + blockDim.x*blockIdx.x;
  if (e < elemcount) {
    double max = 0.;
    double3 rv;
    for (int n = 0 ;n < 3; n++){
      rv = node[3*e+n] - centroid[e];
      if (length(rv) > max) max = length(rv);
      nfar[e] = n;
    }
    
    //element[e]-> radius[e] = max;	//Fraser Eq 3-136
    
    UpdatePlaneCoeff();
	}
}

__global__ inline void CalcSpheresKernel(TriMesh_d *mesh_d) {
  mesh_d->CalcSpheres();
  
}

inline __device__ void TriMesh_d::UpdatePlaneCoeff(){
	//Update pplan
  int i = threadIdx.x + blockDim.x*blockIdx.x;
  if (i < elemcount) { //parallelize by element
    //printf("elnode %f %f %f \n",elnode[3*i+nfar[i]].x,elnode[3*i+nfar[i]].y,elnode[3*i+nfar[i]].z);
    pplane[i] = dot(node[elnode[3*i]+nfar[i]],normal[i]);
    //printf("pplane %.8e \n",pplane[i]);
  }
}

inline __device__ void TriMesh_d::CalcNormals(){
	double3 u, v, w;
  int e = threadIdx.x + blockDim.x*blockIdx.x;
  if (e < elemcount) {
    u = node [elnode[3*e+1]] - node [elnode[3*e]];
    v = node [elnode[3*e+2]] - node [elnode[3*e]];
    w = cross(u,v);
    normal[e] = w/length(w);
    // if (length(normal[e])<1.0e-3)
      // printf("ERROR: ZERO normal. Calc error in element %d\n",e);
    // if (abs(normal[e].y) >1.0e-5 || abs(normal[e].x) > 1.0e-5)
      // printf("CalcNormal %d %.6e %.6e %.6e\n u %.6e %.6e %.6e \n v %.6e %.6e %.6e\n",e, normal[e].x,normal[e].y,normal[e].z,u.x,u.y,u.z,v.x,v.y,v.z);
    normal[e].x = normal[e].y = 0.0;
    normal[e].z = -1.0;
      // //printf("elnodes z coord %.6e %.6e %.6e\n", node[elnode[3*e]].z,node[elnode[3*e+1]].z,node[elnode[3*e+2]].z);
    // }
    //Fraser Eqn 3.34
    //Uj x Vj / |UjxVj|
	}
}

inline __device__ void TriMesh_d::CalcCentroids(){
  int e = threadIdx.x + blockDim.x*blockIdx.x;
  if (e < elemcount)
    centroid[e] = ( node[elnode[3*e]] + node[elnode[3*e+1]] + node[elnode[3*e+2]]) / 3.; 
}

inline __device__ void TriMesh_d::Move(double dt){

	int n = threadIdx.x + blockDim.x*blockIdx.x; //Parallelize by node 
  if ( n < nodecount ){
    double3 vr 	= cross(m_w, node[n]);
    node_v[n] = m_v + vr;
    // for (int i=0;i<3;i++) {
      // if      ((*node[n])(i) < min(i)) min[i] = (*node[n])(i);
      // else if ((*node[n])(i) > max(i)) max[i] = (*node[n])(i);
    // } 

    node[n] += (node_v[n])*dt;
    //printf("after \n");
    //PRINT_V(node[n]); 
  }//n<nodecount
}

__global__ inline void CheckNormalsKernel(TriMesh_d *mesh_d){
  mesh_d->CheckNormals();
}

inline __device__ void TriMesh_d::CheckNormals(){
  int e = threadIdx.x + blockDim.x*blockIdx.x;
  //printf("CheckNormals: %d, elemcount %d\n", e, elemcount);
  if (e < elemcount){
    printf("%d %f %f %f\n", e, normal[e].x,normal[e].y,normal[e].z);
  }  
}

};