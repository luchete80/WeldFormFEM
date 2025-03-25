#include <iostream>
#include <vector>
#include <array>



#include <vector>
#include <array>
#include <Vec3D.h>

// Function to compute barycentric coordinates
std::array<double, 3> barycentric_coordinates(const std::array<double, 2>& p,
                                              const std::array<double, 2>& p0,
                                              const std::array<double, 2>& p1,
                                              const std::array<double, 2>& p2) {
    double denominator = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
    double lambda1 = ((p1[0] - p[0]) * (p2[1] - p[1]) - (p2[0] - p[0]) * (p1[1] - p[1])) / denominator;
    double lambda2 = ((p2[0] - p[0]) * (p0[1] - p[1]) - (p0[0] - p[0]) * (p2[1] - p[1])) / denominator;
    double lambda3 = 1.0 - lambda1 - lambda2;
    return {lambda1, lambda2, lambda3};
}

// Function to interpolate scalar values
double interpolate_scalar(const std::array<double, 2>& p,
                          const std::array<double, 2>& p0, const std::array<double, 2>& p1, const std::array<double, 2>& p2,
                          double scalar0, double scalar1, double scalar2) {
    auto lambdas = barycentric_coordinates(p, p0, p1, p2);
    return lambdas[0] * scalar0 + lambdas[1] * scalar1 + lambdas[2] * scalar2;
}

// Function to interpolate scalar values
Vec3D interpolate_vector (const std::array<double, 2>& p,
                          const std::array<double, 2>& p0, const std::array<double, 2>& p1, const std::array<double, 2>& p2,
                          Vec3D v0, Vec3D v1, Vec3D v2 ) {
    auto lambdas = barycentric_coordinates(p, p0, p1, p2);
    
    Vec3D ret = lambdas[0]*v0 + lambdas[1]*v1 + lambdas[2]*v2;
    return ret;
    //return lambdas[0] * scalar0 + lambdas[1] * scalar1 + lambdas[2] * scalar2;
}

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



#define COPY_NEW (VAR, n)     fnew[n].VAR = interpolate_vector(tgt_nodes[n], pp[0], pp[1], pp[2], \
                                                                                nnpoint[2]->New->VAR,\ 
                                                                                nnpoint[1]->New->VAR, \
                                                                                nnpoint[2]->New->VAR);    


inline void Interp_NodalField(NodalField *fnew,
                              std::array<double, 3> &lambdas,
                              NodalField *o0, NodalField *o1, NodalField *o2 ){
  
  
  fnew->disp = interp_vector(lambdas, o0->disp, o1->disp, o2->disp);    
  fnew->delta_disp = interp_vector(lambdas, o0->delta_disp, o1->delta_disp, o2->delta_disp);    
  
  fnew->ro =   interp_scalar(lambdas, o0->ro, o1->ro, o2->ro);   
  fnew->dro =   interp_scalar(lambdas, o0->dro, o1->ro, o2->dro);   
  
  fnew->mat_speed =   interp_vector(lambdas, o0->mat_speed, o1->mat_speed, o2->mat_speed);   
  
  fnew->dmat_speed = interp_vector(lambdas, o0->dmat_speed, o1->dmat_speed, o2->dmat_speed);  
  fnew->fe = interp_vector(lambdas, o0->fe, o1->fe, o2->fe);  
  
  fnew->e  =   interp_scalar(lambdas, o0->e, o1->e, o2->e);  
  fnew->de =   interp_scalar(lambdas, o0->de, o1->de, o2->de);  
  
  fnew->T    =   interp_scalar(lambdas, o0->T, o1->T, o2->T); 
  fnew->flux =   interp_vector(lambdas, o0->flux, o1->flux, o2->flux); 
  
}


void ReMesh::Generate_MMG(){


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

  np = this->getNodesNumber();



  // CHANGE DIM for 3D Tetra
  nt = nquad = 0;
  na = 0;
  for (int e = 0; e < this->getElementsNumber(); e++) {
      if (this->getElement(e)->getNumberOfNodes() == 4) 
          nt++; // count tetrahedral elements
      else 
          nquad++;
  }
  nt += 2 * nquad; // If needed, handle splits for other elements
  cout << "Number of tetrahedrons: " << nt << ", quads: " << nquad << endl;

  cout << "Structure Node count " << endl;

  if (MMG3D_Set_meshSize(mmgMesh, np, nt, nquad, na) != 1)
      cout << "ERROR ALLOCATING MESH" << endl;
  else 
      cout << "MESH CREATED OK" << endl;
  cout << "Number of points: " << mmgMesh->na << endl;

  // In API_functions
  cout << "struct nodecount " << Global_Structure->getNodesNumber() << endl;
  cout << "Node count " << mmgMesh->np << endl;
  cout << "Mesh node 0 " << mmgMesh->point[0].c[0] << endl;

  if (MMG3D_Chk_meshData(mmgMesh, mmgSol) != 1) 
      exit(EXIT_FAILURE);
  else 
      cout << "Initial Mesh check succeed " << endl;

  // Edge initialization
  int *edges = new int[3 * na];

  // Set vertices for 3D mesh
  for (int n = 0; n < np; n++) {
      if (!MMG3D_Set_vertex(mmgMesh, Global_Structure->getNode(n)->coords(0),
                                          Global_Structure->getNode(n)->coords(1),
                                          Global_Structure->getNode(n)->coords(2),
                                          NULL, n + 1))
          cout << "ERROR ALLOCATING NODE " << n << endl;
  }
  cout << "Vertices allocated" << endl;

  // Set tetrahedrons
  cout << "SETTING TETRAHEDRONS " << endl;
  for (int e = 0; e < this->getElementsNumber(); e++) {
      if (this->getElement(e)->getNumberOfNodes() == 4) {
          MMG3D_Set_tetrahedron(mmgMesh, Global_Structure->getElement(e)->nodes(0)->Id + 1,
                                              Global_Structure->getElement(e)->nodes(1)->Id + 1,
                                              Global_Structure->getElement(e)->nodes(2)->Id + 1,
                                              Global_Structure->getElement(e)->nodes(3)->Id + 1,
                                              NULL, e + 1);
      }
  }

  // Solution setup
  if (MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, np, MMG5_Scalar) != 1)
      exit(EXIT_FAILURE);
  for (int k = 1; k <= np; k++) {
      if (MMG3D_Set_scalarSol(mmgSol, 0.8 - Global_Structure->getNode(k - 1)->getNodalValue("plasticStrain", 0), k) != 1) 
          exit(EXIT_FAILURE);
  }

  // Set parameters (e.g., max edge size)
  MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmax, 0.1);

  // Higher verbosity level
  // MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_verbose, 5);

  // Remesh
  int ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);

  // Get mesh size after remesh
  if (MMG3D_Get_meshSize(mmgMesh, &np, &nt, NULL, &na) != 1)  
      exit(EXIT_FAILURE); 
  cout << "New node count " << np << endl;

  // Set up for corner and required elements
  corner = (int*)calloc(np + 1, sizeof(int));
  if (!corner) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
  }

  required = (int*)calloc(MAX4(np, 0, nt, na) + 1, sizeof(int));
  if (!required) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
  }

  ridge = (int*)calloc(na + 1, sizeof(int));
  if (!ridge) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
  }

  // Copy to destination
  std::vector<std::array<double, 3>> tgt_nodes(np);
  std::vector<std::array<int, 4>> tgt_tetras(nt);
  std::vector<double> tgt_scalar(np);

  nreq = 0;
  nc = 0;

  for (k = 1; k <= np; k++) {
      if (MMG3D_Get_vertex(mmgMesh, &(Point[0]), &(Point[1]), &(Point[2]),
                          &ref, &(corner[k]), &(required[k])) != 1)
          exit(EXIT_FAILURE);

      std::array<double, 3> p0 = {Point[0], Point[1], Point[2]};
      tgt_nodes[k - 1] = p0;
      
      if (corner[k])  
          nc++;
      if (required[k])  
          nreq++;
  }

  // Set up for tetrahedron data
  corner = (int*)calloc(np + 1, sizeof(int));
  if (!corner) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
  }

  Element* el4 = new ElTetra4n3D();
  Global_Structure->setDefaultElement(el4);

  cout << "NODE NUMBER " << Global_Structure->getNodesNumber() << endl;

  cout << "OVERALL tetrahedron count " << nt << endl;
  int nt_corr = 0;
  for (int tetra = 0; tetra < mmgMesh->nt; tetra++) {
      bool error = false;
      int ierror, terr;

      if (!error) {
          MMG5_int Tetra[4];
          int ref;

          MMG3D_Get_tetrahedron(mmgMesh, &(Tetra[0]), &(Tetra[1]), &(Tetra[2]), &(Tetra[3]), &ref, &(required[tetra + 1]));

          std::array<int, 4> ta = {Tetra[0] - 1, Tetra[1] - 1, Tetra[2] - 1, Tetra[3] - 1};
          tgt_tetras[tetra] = ta;                                    

          nt_corr++;
      }
  }

  // Mapping step
  std::vector<NodalField> fnew(np);
  std::vector<NodalField> fcur(np);

  int nf_nodes = 0;
  for (int n = 0; n < np; n++) {
      bool found = false;
      int i = 0;
      while (i < Global_Structure->getElementsNumber() && !found) {
          std::vector<std::array<int, 4>> conn = {{0, 1, 2, 3}};
          int pass = 1;
          
          for (int cp = 0; cp < pass; cp++) {
              Node* nnpoint[4];
              std::vector<std::array<double, 3>> pp(4);
              
              for (int p = 0; p < 4; p++) {
                  nnpoint[p] = Global_Structure->getElement(i)->nodes(conn[cp][p]);
                  pp[p] = {nnpoint[p]->coords(0), nnpoint[p]->coords(1), nnpoint[p]->coords(2)};
              }

              std::array<double, 4> lambdas = barycentric_coordinates(tgt_nodes[n], pp[0], pp[1], pp[2], pp[3]);

              if (lambdas[0] >= -5.0e-2 && lambdas[1] >= -5.0e-2 && lambdas[2] >= -5.0e-2 && lambdas[3] >= -5.0e-2) {
                  double scalar0 = nnpoint[0]->New->disp(1);
                  double scalar1 = nnpoint[1]->New->disp(1);
                  double scalar2 = nnpoint[2]->New->disp(1);
                  double scalar3 = nnpoint[3]->New->disp(1);
                  tgt_scalar[n] = interpolate_scalar(tgt_nodes[n], pp[0], pp[1], pp[2], pp[3], scalar0, scalar1, scalar2, scalar3);
                  found = true;
              }
          }
          i++;
      }
  }

  delete[] corner;
  delete[] required;
  delete[] ridge;
  delete[] edges;

}
