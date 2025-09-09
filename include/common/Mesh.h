#ifndef _MESH_CUH_
#define _MESH_CUH_

//namespace SPH{

#include "defs.h"

#include <stdio.h>

#ifdef  CUDA_BUILD
#include <cuda.h>
#include "vector_math.h"
#endif

#include "utils.h"
#include "NastranReader.h"

namespace MetFEM{

//////////////// FUTURE 


//~ class TriMeshData {  // Pure SoA for performance
//~ public:
    //~ // Geometry
    //~ double3* nodes;
    //~ double3* node_velocities;
    //~ int* elnodes;          // 3 indices per element
    //~ double* pplane;
    
    //~ // Element data
    //~ double3* centroids;
    //~ double3* normals;
    //~ int* nfar;
    
    //~ // Mesh properties
    //~ double3 m_v, m_w;
    //~ int nodecount, elemcount;
    //~ int id;

    //~ int* node_offsets;           // Where each mesh's nodes start [mesh_count+1]
    //~ int* elem_offsets;           // Where each mesh's elements start [mesh_count+1]

    //~ int mesh_count;
    
    //~ // Allocate/deallocate methods
    //~ void Allocate(int node_count, int elem_count);
    //~ void Free();
//~ };

///// TRIMESH DATA
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
  
  // TODO: CHANGE TO NEW TriMesh class
  double3             m_v,m_w;
  int nodecount, elemcount;
  int                 id;
  
  int *ele_mesh_id;
  double *T;        //if constant
  double T_const;   //FOR ALL MATS
  //int *nod_mesh_id;

  double3 *react_force;   //Per mesh
  double  *react_p_force= nullptr; //Per mesh
  double  cont_area;
  double  cfn;
  
  //// DATA FOR SEVERAL MESHES
  int* node_offsets = nullptr;
  int* elem_offsets = nullptr;

  //double3* mesh_velocities;  // Per-mesh velocity [mesh_count]
  //double3* mesh_rotations;   // Per-mesh rotation [mesh_count]

  // Counters
  int total_nodes;
  int total_elems;
  int mesh_count = 0;
  int allocated_meshes;
  // Capacity tracking
  int current_node_capacity = 0;
  int current_elem_capacity = 0;

  double *mu_sta, *mu_dyn;
  double heat_cond;
  
  ///////
	
	//double							v;						//Constant Uniform v
	TriMesh_d();
  TriMesh_d(NastranReader &nr, bool flipnormals);
	void AxisPlaneMesh(const int &id, const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens);
	void AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens);
  
  void SetMeshVel(const double3 &v) {m_v = v;} //Like in WeldForm CPU version
  void SetVel(const double3 &v) {m_v = v;} //Like in WeldForm CPU version
    
  ///ONLY FOR NEW MESH ADDED
  void SetNodesVel(const double3 &v, int start_node = 0, int end_node = -1);
  
	inline void ApplyConstVel(const double3 &v);
	inline void CalcCentroidVelFromNodes();
	inline dev_t void UpdatePlaneCoeff();
	inline void UpdatePos(const double &dt);
	inline dev_t void Move(double dt);
  inline dev_t void CalcNormals();
	inline dev_t void CalcSpheres();
	inline dev_t void CalcCentroids();
  inline dev_t void CheckNormals();



  int GetEleMeshId(int elem_id) const;
      
  
  ///// RESIZING FUNCTIONS
  void ResizeStorage(int new_mesh_capacity);
  int ResizeNodeData(int new_capacity);
  int ResizeElementData(int new_capacity);
  void CopyAndOffsetMeshData(const TriMesh_d& new_mesh, 
                                    int node_offset, int elem_offset);
  void AddMesh(const TriMesh_d& new_mesh); //ResizeStorage should only handle mesh-level metadata (offsets, velocities, etc.)
  
  // Destructor
  ~TriMesh_d() {
      free_t(node);
      free_t(node_v);
      free_t(elnode);
      free_t(centroid);
      free_t(normal);
      free_t(pplane);
      free_t(nfar);
      //free_t(nod_mesh_id);
      free_t(ele_mesh_id);
	  free_t(react_force);
	  free_t(react_p_force);
      //free_t(node_offsets;
      //free_t(elem_offsets;
  }
  
};


/// TO CHCK TO WHICH MESH BELONG AN ELEMENT
//~ int mesh_id = -1;
//~ for (int i = 0; i < mesh_count; i++) {
    //~ if (e >= elem_offsets[i] && e < elem_offsets[i+1]) {
        //~ mesh_id = i;
        //~ break;
    //~ }
//~ }

///// INDIVIDUAL MESH INTERFACE /////
class TriMesh {
private:
    TriMesh_d* storage;  // Shared storage
    int mesh_id;         // Index in storage
    
public:
    // Constructor creates new mesh in storage
    TriMesh(TriMesh_d* storage) : storage(storage) {
        if(storage->mesh_count >= storage->allocated_meshes) {
            storage->ResizeStorage(storage->allocated_meshes * 2 + 1);
        }
        mesh_id = storage->mesh_count++;
    }
    
    // Mesh manipulation
    void SetVelocity(double3 v) {
        //storage->mesh_velocities[mesh_id] = v;
    }
    
    // void AddNodes(const std::vector<double3>& nodes) {
        // // [Implementation would resize storage and add nodes]
        // // Updates node_offsets[mesh_id+1]
    // }
    
    // void AddElements(const std::vector<int>& elements) {
        // // [Implementation would add elements with proper offsets]
        // // Updates elem_offsets[mesh_id+1]
    // }
    
    // // Accessors
    // int GetNodeCount() const {
        // return storage->node_offsets[mesh_id+1] - storage->node_offsets[mesh_id];
    // }
    
    // [Other mesh-specific methods]
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
    //~ normal[e].x = normal[e].y = 0.0;
    //~ normal[e].z = -1.0;
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
    //printf("Element %d pplane %.8e, normalz %.4e \n",e, pplane[e],normal[e].z);
  }
}


inline dev_t void  TriMesh_d::Move(double dt){

	//int n = threadIdx.x + blockDim.x*blockIdx.x; //Parallelize by node 
  for (int n=0;n<nodecount;n++)
  if ( n < nodecount ){
    //double3 vr 	= cross(m_w, node[n]);
    //node_v[n] = m_v + vr;
    //node_v[n] = m_v;
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
      rv = node[elnode[3*e+n]] - centroid[e];
      if (length(rv) > max) max = length(rv);
      nfar[e] = n;
    }
    //printf("Element %d nfar %d\n",e,nfar[e]);
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





//~ class MultiMeshManager {
//~ private:
    //~ TriMeshData* device_mesh_data;  // Single allocation for all mesh data
    //~ int mesh_count;
    
//~ public:
      //~ #ifdef CUDA_BUILD
    //~ // Batch operations
    //~ __global__ void UpdateAllMeshes(double dt) {

        //~ int mesh_idx = blockIdx.x;
        //~ if (mesh_idx < mesh_count) {
            //~ TriMeshData& mesh = device_mesh_data[mesh_idx];
            //~ // Implement movement kernel directly on SoA data
            //~ for (int n=0; n<mesh.nodecount; n++) {
                //~ mesh.node_velocities[n] = mesh.m_v;
                //~ mesh.nodes[n] += mesh.node_velocities[n] * dt;
            //~ }
        //~ }
    //~ }
     //~ #endif
//~ };


};//SPH


#endif
