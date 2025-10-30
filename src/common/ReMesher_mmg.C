/*************************************************************************/
/*  ReMesher_mmg.C                                                     */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                              */
/*  https://www.opensourcemech.com                                                                */
/*                                                                       */
/*  Copyright (c) 2025-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/

#include <iostream>
#include <vector>
#include <array>



#include <vector>
#include <array>
#include <Vec3D.h>

#include "mmg/mmg3d/libmmg3d.h"

#include "ReMesher.h"
#include "Domain_d.h"
#include <iostream>

using namespace std;

// Function to compute barycentric coordinates

// Function to compute barycentric coordinates for a 3D tetrahedron
std::array<double, 4> barycentric_coordinates(const std::array<double, 3>& p,
                                              const std::array<double, 3>& p0,
                                              const std::array<double, 3>& p1,
                                              const std::array<double, 3>& p2,
                                              const std::array<double, 3>& p3) {
    // Compute volume of the tetrahedron
    std::array<double, 3> v0 = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    std::array<double, 3> v1 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
    std::array<double, 3> v2 = {p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]};
    
    double detT = v0[0] * (v1[1] * v2[2] - v1[2] * v2[1])
                - v0[1] * (v1[0] * v2[2] - v1[2] * v2[0])
                + v0[2] * (v1[0] * v2[1] - v1[1] * v2[0]);

    if (detT == 0.0) {
        std::cerr << "Degenerate tetrahedron encountered!" << std::endl;
        return {-1.0, -1.0, -1.0, -1.0}; // Invalid barycentric coordinates
    }


 if (std::abs(detT) < 1e-12) {  // Tolerance to avoid floating-point issues
        std::cerr << "Degenerate tetrahedron encountered!" << std::endl;
        return {-1.0, -1.0, -1.0, -1.0}; // Invalid coordinates
    }

    // Compute the barycentric coordinates
    std::array<double, 3> vp = {p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]};

    double d0 = vp[0] * (v1[1] * v2[2] - v1[2] * v2[1])
              - vp[1] * (v1[0] * v2[2] - v1[2] * v2[0])
              + vp[2] * (v1[0] * v2[1] - v1[1] * v2[0]);

    double d1 = v0[0] * (vp[1] * v2[2] - vp[2] * v2[1])
              - v0[1] * (vp[0] * v2[2] - vp[2] * v2[0])
              + v0[2] * (vp[0] * v2[1] - vp[1] * v2[0]);

    double d2 = v0[0] * (v1[1] * vp[2] - v1[2] * vp[1])
              - v0[1] * (v1[0] * vp[2] - v1[2] * vp[0])
              + v0[2] * (v1[0] * vp[1] - v1[1] * vp[0]);

    double d3 = v0[0] * (v1[1] * v2[2] - v1[2] * v2[1])
              - v0[1] * (v1[0] * v2[2] - v1[2] * v2[0])
              + v0[2] * (v1[0] * v2[1] - v1[1] * v2[0]);

    double w0 = d0 / detT;
    double w1 = d1 / detT;
    double w2 = d2 / detT;
    double w3 = 1.0 - (w0 + w1 + w2);  // Enforce sum-to-one constraint

    return {w0, w1, w2, w3};
}

// Function to interpolate scalar values at the nodes of a tetrahedron
double interpolate_scalar(const std::array<double, 3>& p,
                          const std::array<double, 3>& p0, const std::array<double, 3>& p1, 
                          const std::array<double, 3>& p2, const std::array<double, 3>& p3,
                          double scalar0, double scalar1, double scalar2, double scalar3) {
    auto lambdas = barycentric_coordinates(p, p0, p1, p2, p3);
    return lambdas[0] * scalar0 + lambdas[1] * scalar1 + lambdas[2] * scalar2 + lambdas[3] * scalar3;
}
// Function to interpolate scalar values
std::array<double, 3> interpolate_vector(const std::array<double, 3>& p,
                                         const std::array<double, 3>& p0, const std::array<double, 3>& p1,
                                         const std::array<double, 3>& p2, const std::array<double, 3>& p3,
                                         std::array<double, 3> v0, std::array<double, 3> v1,
                                         std::array<double, 3> v2, std::array<double, 3> v3) {
    auto lambdas = barycentric_coordinates(p, p0, p1, p2, p3);
    return {
        lambdas[0] * v0[0] + lambdas[1] * v1[0] + lambdas[2] * v2[0] + lambdas[3] * v3[0],
        lambdas[0] * v0[1] + lambdas[1] * v1[1] + lambdas[2] * v2[1] + lambdas[3] * v3[1],
        lambdas[0] * v0[2] + lambdas[1] * v1[2] + lambdas[2] * v2[2] + lambdas[3] * v3[2]
    };
}

/*

//WITHOUT CALCULATING AGAIN INTERNAL COORDS
double interp_scalar(std::array<double, 3> &lambdas,
                          double scalar0, double scalar1, double scalar2 ) {
    return lambdas[0] * scalar0 + lambdas[1] * scalar1 + lambdas[2] * scalar2;
}

inline Vec3D interp_vector (std::array<double, 3> &lambdas, //3 coordinates
                          Vec3D v0, Vec3D v1, Vec3D v2 ) {
    
    Vec3D ret = lambdas[0]*v0 + lambdas[1]*v1 + lambdas[2]*v2;
    return ret;
    //return lambdas[0] * scalar0 + lambdas[1] * scalar1 + lambdas[2] * scalar2;
}

*/

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MAX4(a,b,c,d) (MAX(MAX(a,b), MAX(c,d)))

namespace MetFEM{
void ReMesher::Generate_mmg(){

  cout << "Generate mmg"<<endl;
  
  int nreq, ref, nr, nc, *corner, *required, *ridge;  
  MMG5_int Tetra[4], Edge[6], k;
  double Point[3];

  int np, nt, na, nquad;

  MMG5_pMesh mmgMesh;
  MMG5_pSol mmgSol;

  mmgMesh = NULL;
  mmgSol = NULL;
  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
                  MMG5_ARG_end);

  //np = this->getNodesNumber();
  np = m_dom->m_node_count;


  // CHANGE DIM for 3D Tetra
  nt = nquad = 0;
  na = 0;
  for (int e = 0; e < m_dom->m_elem_count; e++) {
      //if (this->getElement(e)->getNumberOfNodes() == 4) 
          nt++; // count tetrahedral elements
      //else 
      //    nquad++;
  }
  //nt += 2 * nquad; // If needed, handle splits for other elements
  cout << "Number of tetrahedrons: " << nt << endl;

  cout << "Structure Node count " << endl;

// int MMG3D_Set_meshSize(MMG5_pMesh mesh, int np, int nt, int nprisms, int nquads, int nhexs, int nedges);

// mesh	MMG5_pMesh	Pointer to your MMG mesh structure.
// np	int32_t	Number of points (vertices) in the mesh.
// nt	int32_t	Number of tetrahedra.
// nprisms	int32_t	Number of prisms (wedge-shaped 3D elements).
// nquads	int32_t	Number of quadrilateral faces (used for boundary representation).
// nhexs	int32_t	Number of hexahedra (not often used unless you have hex meshes).
// nedges	int32_t	Number of edges (only needed if you define edges explicitly, e.g., for BCs).


  if (MMG3D_Set_meshSize(mmgMesh, np, nt, 0,nquad, 0,na) != 1)
      cout << "ERROR ALLOCATING MESH" << endl;
  else 
      cout << "MESH CREATED OK" << endl;
  cout << "Number of points: " << mmgMesh->na << endl;

  // In API_functions
  //cout << "struct nodecount " << Global_Structure->getNodesNumber() << endl;
  
  cout << "Node count " << mmgMesh->np << endl;
  cout << "Mesh node 0 " << mmgMesh->point[0].c[0] << endl;

  if (MMG3D_Chk_meshData(mmgMesh, mmgSol) != 1) 
      exit(EXIT_FAILURE);
  else 
      cout << "Initial Mesh check succeed " << endl;

  // Edge initialization


  // Set vertices for 3D mesh
  for (int n = 0; n < np; n++) {
      if (!MMG3D_Set_vertex(mmgMesh, m_dom->x[3*n],
                                     m_dom->x[3*n+1],
                                     m_dom->x[3*n+2],
                                          NULL, n + 1))
          cout << "ERROR ALLOCATING NODE " << n << endl;
  }
  cout << "Vertices allocated" << endl;

  // Set tetrahedrons
  cout << "SETTING TETRAHEDRONS " << endl;
  for (int e = 0; e < m_dom->m_elem_count; e++) {

          MMG3D_Set_tetrahedron(mmgMesh,      m_dom->m_elnod[4*e]+ 1,
                                              m_dom->m_elnod[4*e+1] + 1,//ISSUMMING 1 OR NOT?
                                              m_dom->m_elnod[4*e+2] + 1,
                                              m_dom->m_elnod[4*e+3] + 1,
                                              NULL, e + 1);
      
  }

  // Solution setup
  if (MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, np, MMG5_Scalar) != 1)
      exit(EXIT_FAILURE);
    
  // for (int k = 1; k <= np; k++) {
      // if (MMG3D_Set_scalarSol(mmgSol, 0.8 - m_dom->pl_strain[k-1], k) != 1) 
          // exit(EXIT_FAILURE);
  // }

  // for (int k = 1; k <= np; k++) {
      // double refinement_factor = 1.0 - m_dom->pl_strain[k-1];  // more strain → smaller element
      // double h_target = m_dom->m_remesh_length * refinement_factor;
      // if (h_target < 0.6*m_dom->m_remesh_length) h_target = 0.6*m_dom->m_remesh_length; // optional min size
      // MMG3D_Set_scalarSol(mmgSol, h_target, k);
  // }
    
  for (int k = 1; k <= np; k++) {
      if (MMG3D_Set_scalarSol(mmgSol, 0.5, k) != 1) 
          exit(EXIT_FAILURE);
  }

  // Set parameters (e.g., max edge size)
  MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmax, 0.5);

  // Higher verbosity level
  // MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_verbose, 5);

  MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, np, MMG5_Scalar);


  // for (int k = 1; k <= np; k++)
    // MMG3D_Set_scalarSol(mmgSol, 2.0*m_dom->m_remesh_length, k);  // uniform sizing via scalar field

  double max_ps = 0;
  int max_ps_elem;
  double min_el_size=1000;
  
  double min_size = m_dom->m_remesh_min_frac * m_dom->m_remesh_length;  //Tamaño mínimo en zonas de alta deformación
  double max_size = m_dom->m_remesh_max_frac * m_dom->m_remesh_length;  //Tamaño máximo en zonas sin deformación
  
  double *ps_nod = new double[m_dom->m_node_count];
  
  //TODO: IN GPU SHOULD BE WITH HOST ARRAY INSTEAD OF DEVICE vars 
  //avgScalarDom(m_dom->pl_strain,ps_nod,6,m_dom)

  cout << "BETA: "<<m_dom->m_remesh_beta<<endl;
  cout << "EPS REF: "<<m_dom->m_remesh_eps_ref<<endl;
   
  if (m_dom->m_remesh_type == 0) {
    cout << "REMESH TYPE IS LINEAR"<<endl;
    //////////////////////// HYPERBOLIC
    //~ cout << "Setting size "<<endl;
    for (int k = 1; k <= np; k++) { //REMEMBER: np = m_dom->m_node_count;
      
        //Calcular factor de refinamiento basado en deformación plástica
        //Más deformación → elementos más pequeños (refinamiento)
        double ps=0;
        for (int e=0; e<m_dom->m_nodel_count[k-1];e++) {
          int eglob   = m_dom->m_nodel     [m_dom->m_nodel_offset[k-1]+e];
          ps+=m_dom->pl_strain[eglob];
        }
        ps/=m_dom->m_nodel_count[k-1];


        //double plastic_strain = m_dom->pl_strain[k-1];
        double plastic_strain = ps;
        
        if (plastic_strain>max_ps){
          max_ps = plastic_strain;
          max_ps_elem = k-1;
          }
        
        //Definir rangos de refinamiento

        
        //Mapear la deformación plástica al tamaño de elemento
        //Puedes ajustar esta función según tus necesidades
        double refinement_factor = 1.0 / (1.0 + 2.0 * plastic_strain);  //Más strain → factor menor
        double h_target = min_size + (max_size - min_size) * refinement_factor;
        
        //Limitar el tamaño mínimo y máximo
        if (h_target < min_size) h_target = min_size;
        if (h_target > max_size) h_target = max_size;
        
        if (h_target<min_el_size) min_el_size = h_target;
        
        if (MMG3D_Set_scalarSol(mmgSol, h_target, k) != 1) 
            exit(EXIT_FAILURE);
    }
  } else if (m_dom->m_remesh_type == 1) {
    cout << "REMESH TYPE IS EXPONENTIAL"<<endl;
  //~ ///// EXPONENTIAL
  for (int k = 1; k <= np; k++) { //REMEMBER: np = m_dom->m_node_count;
    double h0        = m_dom->m_remesh_length;

    int    np_nodes  = np;


      //double plastic_strain = m_dom->pl_strain[k-1];

      double plastic_strain=0;
      for (int e=0; e<m_dom->m_nodel_count[k-1];e++) {
        int eglob   = m_dom->m_nodel     [m_dom->m_nodel_offset[k-1]+e];
        plastic_strain+=m_dom->pl_strain[eglob];
      }
      plastic_strain/=m_dom->m_nodel_count[k-1];
      

      // Normalizar y clamp entre 0 y 1
      double s = plastic_strain / m_dom->m_remesh_eps_ref;
      if (s < 0.0) s = 0.0;
      if (s > 1.0) s = 1.0;

      // Exponencial: factor ∈ (0,1], más s → factor → pequeño
      double refinement_factor = std::exp(-m_dom->m_remesh_beta * s); // cuando s=0 -> 1, s=1 -> exp(-beta)

      double h_target = min_size + (max_size - min_size) * refinement_factor;

      // Limitar por seguridad (por si acaso)
      if (h_target < min_size) h_target = min_size;
      if (h_target > max_size) h_target = max_size;

      if (MMG3D_Set_scalarSol(mmgSol, h_target, k) != 1) 
          exit(EXIT_FAILURE);
                    
    } 
      
    } else {
      
      cout << "NO REMESH TYPE!!!!"<<endl;
      exit(1);
      }    
    //~ // parámetros
    //~ double h0               = m_dom->m_remesh_length;
    //~ double min_frac         = 0.35;
    //~ double max_frac         = 2.0;
    //~ double eps_ref          = 0.1;
    //~ double beta             = 4.0;
    //~ double max_resize_ratio = 1.25; // no cambiar más del 25% por remesh
    //~ double quality_thresh   = 0.4;
    //~ double ps_activate      = 0.02;

    //~ for (int k = 1; k <= np; k++) {
      //~ // promedio de plastic strain alrededor del nodo
      //~ double plastic_strain = 0.0;
      //~ for (int e=0; e<m_dom->m_nodel_count[k-1]; e++) {
        //~ int eglob = m_dom->m_nodel[m_dom->m_nodel_offset[k-1]+e];
        //~ plastic_strain += m_dom->pl_strain[eglob];
      //~ }
      //~ plastic_strain /= m_dom->m_nodel_count[k-1];

      //~ // calcula calidad local (0..1) según tu criterio (angles, J, aspect)
      //~ double quality = computeLocalQualityForNode(k); // implementar: 1=buena,0=mala

      //~ // normaliza y limita s
      //~ double s = plastic_strain / eps_ref;
      //~ if (s < 0.0) s = 0.0;
      //~ if (s > 1.0) s = 1.0;

      //~ double min_size = min_frac * h0;
      //~ double max_size = max_frac * h0;

      //~ // mezcla strain y distorsión: si quality baja -> aumentar refinamiento
      //~ double q_factor = 1.0 - quality; // 0 si buena calidad, 1 si mala
      //~ double combined = 0.7 * s + 0.3 * q_factor; // pesos: ajustar (0.7,0.3)

      //~ double refinement_factor = std::exp(-beta * combined);
      //~ double h_target = min_size + (max_size - min_size) * refinement_factor;

      //~ // si plastic strain muy bajo y calidad buena, no remalles (previene ruido)
      //~ if (plastic_strain < ps_activate && quality > 0.7) {
        //~ h_target = max_size; // mantener malla gruesa
      //~ }

      //~ // limitar cambio respecto h anterior si existe (prev_h[k])
      //~ double prev_h = getPrevHTargetForNode(k); // implementa caché/array
      //~ if (prev_h > 0.0) {
        //~ double min_allowed = prev_h / max_resize_ratio;
        //~ double max_allowed = prev_h * max_resize_ratio;
        //~ if (h_target < min_allowed) h_target = min_allowed;
        //~ if (h_target > max_allowed) h_target = max_allowed;
      //~ }

      //~ // set scalar for MMG
      //~ if (MMG3D_Set_scalarSol(mmgSol, h_target, k) != 1)
        //~ exit(EXIT_FAILURE);
    //~ }



  printf("Max PS %.4e at element %d Min Element size %.4e.Size Limits: min %.4e , max %.4e\n", max_ps,max_ps_elem,min_el_size,min_size,max_size);
  
  

  //~ // --- Deterministic remesh parameters ---

  MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_verbose, 5);
  

  //~ MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmax, m_dom->m_remesh_length);
  //~ MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmin, m_dom->m_remesh_length);
  //~ MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hgrad, 1.0);

  //~ MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_noswap, 1);     // No swaps
  //~ MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_IPARAM_noinsert, 1);
  //~ MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_IPARAM_nomove, 1);
  
   // MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_nb_iter, 1);    // Min iter
        

  // Remesh
  int ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);

  // Get mesh size after remesh
// int MMG3D_Get_meshSize(MMG5_pMesh mesh,
                       // int32_t* np,
                       // int32_t* nt,
                       // int32_t* nprisms,// int32_t* nquads,// int32_t* nhexs,// int32_t* nedges);

// mesh	MMG5_pMesh	Pointer to the MMG mesh structure.
// np	int32_t*	Pointer to store number of points (nodes).
// nt	int32_t*	Pointer to store number of tetrahedra.
// nprisms	int32_t*	Pointer to store number of prisms (wedge elements).
// nquads	int32_t*	Pointer to store number of quadrilateral faces (typically for BCs).
// nhexs	int32_t*	Pointer to store number of hexahedra.
// nedges	int32_t*	Pointer to store number of edges (optional).
  if (MMG3D_Get_meshSize(mmgMesh, &np, &nt, NULL, NULL, NULL, &na) != 1)  
      exit(EXIT_FAILURE); 
  cout << "New node count " << np << endl;

  // Set up for corner and required elements
  corner = (int*)calloc(np + 1, sizeof(int));
  if (!corner) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
  }

  //CORRECT FROM HERE
  
  required = (int*)calloc(MAX4(np, 0, nt, na) + 1, sizeof(int));
  if (!required) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
  }

  int *required_nodes = (int*)calloc(np + 1, sizeof(int));
  int *required_tets  = (int*)calloc(nt + 1, sizeof(int));


  // Copy to destination
  std::vector<std::array<double, 3>> tgt_nodes(np);
  std::vector<std::array<int, 4>> tgt_tetras(nt);
  std::vector<double> tgt_scalar(np);


  //m_dom->Free();
  
  //m_dom->m_node_count = np;
  //m_dom->m_elem_count = nt;
  

  //m_dom->SetDimension(m_dom->m_node_count,m_dom->m_elem_count);	 //AFTER CREATING DOMAIN
  
  m_x = new double [3*np];
    
  nreq = 0;
  nc = 0;
  //cout << "Getting vertices "<<endl;
  for (k = 1; k <= np; k++) {
      //cout << "Vertex "<<k<<endl;
      if (MMG3D_Get_vertex(mmgMesh, &(Point[0]), &(Point[1]), &(Point[2]),
                          &ref, &(corner[k]), &(required[k])) != 1)
          exit(EXIT_FAILURE);

      std::array<double, 3> p0 = {Point[0], Point[1], Point[2]};
      tgt_nodes[k - 1] = p0;
      for (int d=0;d<3;d++) m_x[3*(k-1)+d] = Point[d];
      
      if (corner[k])  
          nc++;
      if (required[k])  
          nreq++;
  }
  cout << "Done. "<<endl;
  // Set up for tetrahedron data
  corner = (int*)calloc(np + 1, sizeof(int));
  if (!corner) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
  }


  //malloc_t(m_dom->m_elnod, unsigned int,m_dom->m_elem_count * m_dom->m_nodxelem);
  
  m_elnod = new int[nt * 4];
  
  cout << "OVERALL tetrahedron count " << nt << endl;
  int nt_corr = 0;
  for (int tetra = 0; tetra < nt; tetra++) {
      bool error = false;
      int ierror, terr;

      if (!error) {
          MMG5_int Tetra[4];
          int ref;

          MMG3D_Get_tetrahedron(mmgMesh, &(Tetra[0]), &(Tetra[1]), &(Tetra[2]), &(Tetra[3]), &ref, &(required_tets[tetra + 1]));

          std::array<int, 4> ta = {Tetra[0] - 1, Tetra[1] - 1, Tetra[2] - 1, Tetra[3] - 1};
          tgt_tetras[tetra] = ta;                                    
          for (int en=0;en<4;en++) m_elnod[tetra*4+en] = ta[en];
          //cout << ta[0]<<", "<<ta[1]<<"; "<<ta[2]<<", "<<ta[3]<<endl;
          nt_corr++;
      }
  }

  
  
  cout << "NEW MESH. Done mapping "<<endl;
  cout << "OLD Node count "<<m_dom->m_node_count<<", ELEM COUNT "<<m_dom->m_elem_count<<endl;
  cout << "NEW Node count "<<np<<", ELEM COUNT "<<nt<<endl;
  //memcpy_t(m_dom->x,      m_x , 3*sizeof(double) * m_dom->m_node_count);       
  //memcpy_t(m_dom->m_elnod,  m_elnod, 4*sizeof(int) * m_dom->m_elem_count);  
  
  //m_dom->setNodElem(m_elnod); 

  
    //NEW MESH; SO CAN BE USED OMEGA_H FUNCTIONS
    
    ///// IF USED OMEGA_H
  
  //cout << "Creating NEW mesh "<<endl;
  //create_mesh(m_mesh, x_h, np, (int *)m_elnod, nt);
  //cout << "done"<<endl;
  
  //m_x = new double [3*np];
  //m_elnod = new int[nt * 4];
  


  m_node_count = np;
  m_elem_count = nt;
  cout << "MESH CREATED."<<endl;
 
  malloc_t(m_closest_elem, int,m_elem_count);

  free(corner);
  free(required);
  free(required_tets);


  free(required_nodes);
  required_nodes = nullptr;  // optional but recommended

}


double tet_volume(const double v0[3], const double v1[3], const double v2[3], const double v3[3]) {
    double a[3], b[3], c[3];
    for (int i = 0; i < 3; ++i) {
        a[i] = v1[i] - v0[i];
        b[i] = v2[i] - v0[i];
        c[i] = v3[i] - v0[i];
    }
    double volume = (a[0]*(b[1]*c[2] - b[2]*c[1])
                   - a[1]*(b[0]*c[2] - b[2]*c[0])
                   + a[2]*(b[0]*c[1] - b[1]*c[0])) / 6.0;
    return std::abs(volume);
}

void calcVol(){
    
  int nelts;
  MMG5_pMesh mmgMesh;
  MMG5_pSol mmgSol;
/*

  MMG3D_Get_numberOfElements(mmgMesh, &nelts);

  // Then loop over elements:
  for (int i = 1; i <= nelts; ++i) {
      int verts[4]; // tetrahedron
      double coords[4][3];

      MMG3D_Get_tetrahedron(mmgMesh, &verts[0], &verts[1], &verts[2], &verts[3], NULL);

      for (int j = 0; j < 4; ++j) {
          MMG3D_Get_vertex(mmgMesh, &coords[j][0], &coords[j][1], &coords[j][2], NULL, NULL, verts[j]);
      }
      
      double[3] p1,p2,p3,p4;
      // Compute volume of this tetrahedron
      double v = std::abs(
          Omega_h::tet_volume(
          p1,p2,p3,p4
          
          
          )
    );
  } 
  */ 
  
}




};
