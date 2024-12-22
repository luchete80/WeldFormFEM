#ifndef _MESH_CUH_
#define _MESH_CUH_

//namespace SPH{

#include "defs.h"

#include <stdio.h>

#ifdef  CUDA_BUILD
#include <cuda.h>
#endif

#include "utils.h"

namespace MetFEM{

struct elem_data{
    double3             *centroid,*normal;
};
//Element is not anymore here, is everything flattened in mesh_d class
class TriMesh_d{
	
	public:

	//Element 						elem_data;
	double3 						*node,*node_v; //Positions and veloc, 
	int									*elnode;			//3 per element
	double 							*pplane; //Only for testing
  //Element data, TODO: PASS TO ELEMDATA
  double3             *centroid,*normal;
  int                 *nfar;
  double3             m_v,m_w;
  int nodecount, elemcount;
  int                 id;
  
	
	//double							v;						//Constant Uniform v
	TriMesh_d(){}
	void AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens);
  void SetVel(const double3 &v) {m_v = v;} //Like in WeldForm CPU version
	inline void ApplyConstVel(const double3 &v);
	inline void CalcCentroidVelFromNodes();
	inline dev_t void UpdatePlaneCoeff();
	inline void UpdatePos(const double &dt);
	inline dev_t void Move(double dt);
  inline dev_t void CalcNormals();
	inline dev_t void CalcSpheres();
	inline dev_t void CalcCentroids();
  inline dev_t void CheckNormals();
  
  
};


inline dev_t void   TriMesh_d::CalcCentroids(){
  //int e = threadIdx.x + blockDim.x*blockIdx.x;
  //if (e < elemcount)
    for (int e=0; e<elemcount;e++)
      centroid[e] = ( node[elnode[3*e]] + node[elnode[3*e+1]] + node[elnode[3*e+2]]) / 3.; 
}


inline dev_t void   TriMesh_d::CalcNormals(){
	double3 u, v, w;
  //int e = threadIdx.x + blockDim.x*blockIdx.x;
  //if (e < elemcount) {
  for (int e=0; e<elemcount;e++){
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


inline dev_t void  TriMesh_d::UpdatePlaneCoeff(){
	//Update pplan
  //int i = threadIdx.x + blockDim.x*blockIdx.x;
  //if (i < elemcount) { //parallelize by element
  for (int e=0; e<elemcount;e++){
    //printf("elnode %f %f %f \n",elnode[3*i+nfar[i]].x,elnode[3*i+nfar[i]].y,elnode[3*i+nfar[i]].z);
    pplane[e] = dot(node[elnode[3*e]+nfar[e]],normal[e]);
    //printf("pplane %.8e \n",pplane[e]);
  }
}


inline dev_t void  TriMesh_d::Move(double dt){

	//int n = threadIdx.x + blockDim.x*blockIdx.x; //Parallelize by node 
  for (int n=0;n<nodecount;n++)
  if ( n < nodecount ){
    //double3 vr 	= cross(m_w, node[n]);
    //node_v[n] = m_v + vr;
    node_v[n] = m_v;
    // for (int i=0;i<3;i++) {
      // if      ((*node[n])(i) < min(i)) min[i] = (*node[n])(i);
      // else if ((*node[n])(i) > max(i)) max[i] = (*node[n])(i);
    // } 

    node[n] = node[n] + (node_v[n])*dt;
    //printf("after \n");
    //PRINT_V(node[n]); 
  }//n<nodecount

  //CalcCentroids();
  //CalcNormals();        //From node positions
  //UpdatePlaneCoeff();   //pplane
  
}

inline  dev_t void   TriMesh_d::CalcSpheres(){
	// double max;
  //int e = threadIdx.x + blockDim.x*blockIdx.x;
  //if (e < elemcount) {
  for (int e=0; e<elemcount;e++){
    double max = 0.;
    double3 rv;
    for (int n = 0 ;n < 3; n++){
      rv = node[3*e+n] - centroid[e];
      if (length(rv) > max) max = length(rv);
      nfar[e] = n;
    }
    printf("Element %d nfar %d\n",e,nfar[e]);
    //element[e]-> radius[e] = max;	//Fraser Eq 3-136
    
	}
  UpdatePlaneCoeff();
}

//__global__ inline void MeshUpdateKernel(TriMesh_d *mesh_d, double dt) {
inline  dev_t void   MeshUpdate/*Kernel*/(TriMesh_d *mesh_d, double dt) {
 	mesh_d->Move(dt);
  mesh_d->CalcCentroids();
  mesh_d->CalcNormals();
  mesh_d->UpdatePlaneCoeff(); 
}

#ifdef  CUDA_BUILD

__global__ inline void MeshUpdateKernel(TriMesh_d *mesh_d, double dt);
__global__ inline void CalcSpheresKernel(TriMesh_d *mesh_d);

__global__ inline void CheckNormalsKernel(TriMesh_d *mesh_d);
#endif

};//SPH


#endif
