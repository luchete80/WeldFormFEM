/*************************************************************************/
/*  Contact.C                                                    */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                */
/*  ('https://www.opensourcemech.com',)                                    */
/*                                                                       */
/*  Copyright (c) 2023-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. */ 
/*  See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/



#include "Domain_d.h"
#include "Mesh.h"

namespace MetFEM{
#define MAX_NB_COUNT    20



#define HTOL 1.0e-6
#define DFAC 0.0


void dev_t Domain_d::CalcContactForces(){
  
  #ifndef CUDA_BUILD
	//int i = threadIdx.x + blockDim.x*blockIdx.x;	
  double pxa[m_node_count];
	double min_force_ts_=1000.;
// https://stackoverflow.com/questions/10850155/whats-the-difference-between-static-and-dynamic-schedule-in-openmp

  int nen = (trimesh->dimension == 3) ? 3 : 2;
			
  //max_contact_force = 0.;
	double min_contact_force = 1000.;
	int inside_pairs = 0;
  //printf("test\n");
  /*
  #ifdef BUILD_GPU
  par_loop(i,ext_nodes_count)
  #else
  for  (int i=0; i < ext_nodes_count;i++ )   //i particle is from SOLID domain, j are always rigid 
  #endif
  {
  */
  //~ par_loop(i,m_dim*m_node_count)
    //~ contforce[i]=0.0;
    m_dt_gap_min = 1.0;
    
  double cfn[m_node_count];
  par_loop(i,m_node_count){
    m_mesh_in_contact[i]=-1;
    pxa[i]=0.0;
  }
  par_loop(i,m_node_count){
    //printf("Node %d\n",i);
    //printf("Element count %d\n", trimesh->elemcount);
    double min_dist = 1e10;
    double delta;
    
    if (ext_nodes[i]){ 
      //Search NEGATIVE DISTANCE node ON MESH
      int j=0;
      bool end = false;
      while (!end){//----------------------------------------------------------
      //CHCK isNodeinElement WITH getShapeFunctionAtPoint
        //printf("Node i %d Element  j %d\n", i,j);      

        double d = dot(trimesh->normal[j],getNodePos3(i))  - trimesh->pplane[j];

        double3 x_pred = getPosVec3(i) + getVelVec(i) * dt + getAccVec(i) * dt * dt/2.0;
        
        double dalt =  dot (trimesh->normal[j], x_pred ) - trimesh-> pplane[j];

        //printf("PPLANE %.4e \n", trimesh-> pplane[j]);
        
        //double d = norm2(dist);
        //Original (wrong
        //delta = dot(dist,trimesh->normal[j]);
        //delta = dot(getNodePos3(i) - trimesh->centroid[j],trimesh->normal[j]);

        //~ x_pred = Particles[P1]->x + Particles[P1]->v * deltat + Particles[P1]->a * deltat * deltat/2.0;
        //~ vr_pred = Particles[P1]->v + Particles[P1]->a * deltat - Particles[P2]->v;
        
        //~ dist =  dot (Particles[P2]->normal, x_pred ) - trimesh[m]-> element[Particles[P2]->element] -> pplane;
        
        //ACCORDING TO WANG
			//double3 vr = Particles[P1]->v - Particles[P2]->v;		//Fraser 3-137
			//delta = - dot(trimesh->normal[j], vr);	//Penetration rate, Fraser 3-138        
        //printf("Node %d d %f\n",i, d);
        if (d < 0 /*&& dist < CERTAIN ELEMENT DISTANCE*/){
          //printf ("NODE %d DELTA <0------, d %.4e Z POS %.4e , PPLANE %.4e \n", i, d, getNodePos3(i).z,trimesh->pplane[j]);
         
          //ORIGINAL
          double3 Qj = getPosVec3(i) - d * trimesh->normal[j];
          //~ if (i==40)
            //~ printf("Element %d node 40 qj %f %f %f \n", j,Qj.x,Qj.y,Qj.z);

          bool inside = true;
          if (trimesh->dimension == 3){
            int l=0,n;		   
            //printf("Entering while \n");
            while (l<3 && inside){
              n = l+1;	if (n>2) n = 0;
              //double crit;
              double crit = dot (cross ( trimesh->node[trimesh->elnode[nen*j+n]] - trimesh->node[trimesh->elnode[nen*j+l]],
                                                                  Qj  - trimesh->node[trimesh->elnode[nen*j+l]]),
                                 trimesh->normal[j]);
              double3 test = cross ( trimesh->node[trimesh->elnode[nen*j+n]] - trimesh->node[trimesh->elnode[nen*j+l]],
                                                                  Qj  - trimesh->node[trimesh->elnode[nen*j+l]]);
                                        
              if (crit < 0.0) inside = false;
              l++;
            }
          } else { //dimension == 2
            int l=0, n;
            while (l<2 && inside){
              n = l+1;	if (n>1) n = 0;
              double crit = dot ( trimesh->node[trimesh->elnode[nen*j+n]] 
                                          - trimesh->node[trimesh->elnode[nen*j+l]],
                                          Qj  - trimesh->node[trimesh->elnode[nen*j+l]]);
                                          
              //double crit = dot ( *trimesh[m]->node[e -> node[j]] 
                                          //~ - *trimesh[m]->node[e -> node[i]],
                                          //~ qj  - *trimesh[m]->node[e -> node[i]]);
              if (crit < 0.0) inside = false;
              l++;
            }
                        
            
          }

          if (inside ){
 
              double nodlen = 0.0;
              double nodevol = 0.0;
              // TODO: MAKE FUNCTION FOR THIS
              for (int e=0; e<m_nodel_count[i];e++) { 
                int eglob   = m_nodel     [m_nodel_offset[i]+e]; //Element
                nodlen +=m_elem_length[eglob];
                nodevol+=vol[eglob];
              }
              nodlen /= m_nodel_count[i]; 
              nodevol /=m_nodel_count[i];
              
              //~ printf("Element %d, delta: %.4e\n",j, d);
              //~ printf("PPLANE: %.4e\n",trimesh->pplane[j]);
              //~ printf("MESH NORMAL: %.4e %.4e %.4e\n",trimesh->normal[j].x,trimesh->normal[j].y,trimesh->normal[j].z);
              //~ printf("ELNODES: %d %d %d \n", trimesh->elnode[3*j],trimesh->elnode[3*j+1],trimesh->elnode[3*j+2]);
              //~ printf("MESH POINT: %.4e %.4e %.4e\n",trimesh->node[j].x,trimesh->node[j].y,trimesh->node[j].z);
              //~ printf("BODY POINT : %.4e %.4e %.4e\n",getNodePos3(i).x,getNodePos3(i).y,getNodePos3(i).z);

                            
              double beta_d = 0.1;
              //printf("NODLEN %.4e\n", nodevol);
              //double kcont = 0.8 * mat[0]->Elastic().E() * node_area[i]/nodlen;       
              //double kcont = 0.2 * mat[0]->Elastic().E() * pow(nodevol,0.3333);           
              //double kcont = 0.2 * m_mdiag[i] /(dt*dt);
                            
              // THIS IS USED ALSO FOR CONTACT DAMPING
              double3 v_rel = getVelVec(i); // master assumed static or you can use relative motion     
              double v_reln = dot(v_rel, trimesh->normal[j]);         
              //double F_damp = 0.0;

              //time

              // double delta = -d; // penetración
              // double h_eff = pow(m_voln[n],1./3.); // o radio de contacto, o cbrt(nodevol)

              // double delta_eff = max(delta, 0.2 * h_eff);

              // double vn = fabs(v_reln);

              // if (vn > 1e-8) {
                // double dt_gap = 0.2 * delta_eff / vn;
                // m_dt_gap_min = min(m_dt_gap_min, dt_gap);
                 // printf("delta_eff %.3e, vn %.3e , dt gap %.3e\n", delta_eff, vn, dt_gap);
              // }

              double kcont_geo = mat[0]->Elastic().E() * node_area[i] / nodlen;
              //double kcont_geo = mat[0]->Elastic().BulkMod()*nodevol;
              double kcont_mass = 0.2*m_mdiag[i] / (dt * dt);
              double kcont = std::min(kcont_geo, kcont_mass);
              //printf("KGEO:%.4e KMASS: %.4e\n",kcont_geo,kcont_mass);
              
              //double kcont = kcont_mass;
              double damping_ratio = 0.2;
              double omega = sqrt(kcont / m_mdiag[i]);
              double ccrit = 2.0 * m_mdiag[i] * omega;
              double F_damp = damping_ratio * ccrit * v_reln;

              double F_normal = m_contPF * kcont * d;

              //double F_normal = m_contPF * kcont * eff_penet;

              //double F_damp = beta_d * sqrt(kcont * m_mdiag[i]) * v_reln;
              //printf("KCON geo %f KCON mass %f: \n", kcont, kcont2;
              double3 cf =  - (F_normal  + F_damp) * trimesh->normal[j];
              //~ if (getPosVec3(i).z<0.0001)
                //~ printf("Node %d CF %f %f %f, dist %f mass %f\n",i, cf.x,cf.y,cf.z, d,m_mdiag[i]);
              contforce[m_dim*i] = cf.x;contforce[m_dim*i+1] = cf.y;contforce[m_dim*i+2] = cf.z;
              pxa[i] = p_node[i] * node_area[i];
              //printf("Nodal pressure %.4e\n",p_node[i]);
              //printf("Cont Force %f %f %f \n",cf.x,cf.y,cf.z);
              cfn[i] = norm(cf);
              
              //printf("MESHIN CONTACT %d\n",trimesh->ele_mesh_id[j]);
              m_mesh_in_contact[i]=trimesh->ele_mesh_id[j];
              //m_mesh_in_contact[i]=0;
              //~ if (m_mesh_in_contact[i]>trimesh->mesh_count)
                //~ printf("ERROR!!! \n");
              ////

              double damping = 1.0;  // Tune this (0.01–0.5)
              ////////WANG 2013     
              double3 du = x_pred - getPosVec3(i) - getVelVec(i) * dt ;
              double3 delta_tg = du - dot(du, trimesh->normal[j])*trimesh->normal[j];
              double3 tgforce = kcont * delta_tg;

              
              
               // 1. Relative velocity at contact in tangential direction
              //double3 v_rel = getVelVec(i); // master assumed static or you can use relative motion
              double3 v_tan = v_rel - dot(v_rel, trimesh->normal[j]) * trimesh->normal[j];

              // 2. Incremental tangential displacement (Δut = vt * dt)
              double3 du_tangent = v_tan * dt;

              // 3. Accumulate tangential displacement
              double3 ut_acc = make_vector_t(ut_prev[3*i],ut_prev[3*i+1],ut_prev[3*i+2]);
              ut_acc = ut_acc + du_tangent;
              ut_prev[m_dim*i] += du_tangent.x; ut_prev[m_dim*i+1] += du_tangent.y;ut_prev[m_dim*i+2] += du_tangent.z;
              // 4. Trial tangential force from accumulated displacement
              double3 Ft_trial = -kcont * ut_acc;

              // 5. Coulomb cap

              double Ft_mag = length(Ft_trial);
              double3 Fn = dot(cf,trimesh->normal[j])*trimesh->normal[j];


              ////// THIS IS TO SMOOTH STICK-SLIP TRANSITION
              //~ double v_crit = 1e-3; // Umbral de velocidad (ajustable)  
              //~ double mu_eff = trimesh->mu_dyn[0] + (trimesh->mu_sta[0] - trimesh->mu_dyn[0]) * exp(-norm(v_tan) / v_crit);  
              //~ double Ft_max = mu_eff * norm(Fn);  

              double Ft_max = trimesh->mu_sta[0] * norm(Fn); // norm(Fn) is magnitude of normal force
              

              double3 Ft;

              double Ft_max_static = trimesh->mu_sta[0] * norm(Fn);
              if (Ft_mag <= Ft_max_static) {
                  // 6. Cap force, rescale
                  // Still sticking: Use trial force
                   //~ if (i==0){
                      //~ printf("STATIC FORCE %f %f %f\n",tgforce.x,tgforce.y,tgforce.z);
                     //~ printf("DISP NODE 0 %f %f %f\n", u[3*i],u[3*i+1],u[3*i+2]);
                     //~ print
                     //~ }
                  Ft = Ft_trial;
              }
              else {
                  double Ft_max_dynamic =  trimesh->mu_dyn[0] * norm(Fn);
                  Ft = -Ft_max_dynamic * v_tan/norm(v_tan);  // Use velocity direction now

                  //~ double decay_factor = 0.5; // Ajustable (0 = reset completo, 1 = sin cambio)  
                  //~ for (int d = 0; d < m_dim; d++) ut_prev[m_dim*i + d] *= decay_factor;  
    
                  for (int d=0;d<3;d++)ut_prev[m_dim*i+d] = 0;  // Optional: reset accumulation when sliding
              }              


              // 7. Apply friction force
              contforce[m_dim * i + 0] += Ft.x;
              contforce[m_dim * i + 1] += Ft.y;
              contforce[m_dim * i + 2] += Ft.z;                     


              q_cont_conv[i] = trimesh->heat_cond * node_area[i]*(trimesh->T_const - T[i]);

              
              end = true;//JUST ONE MASTER ELEMENT PER SLAVE NODE
            
            } //INSIDE

          
        } //if d< 0
        //cout << "MESH ELEMENT "<<j <<endl;
        j++;
        if (j==trimesh->elemcount)
          end = true;
      }//while !end
  
    }//external nodes
  } //i<first fem index
    
    trimesh->react_force[0].z = 0.0;
   for (int i=0;i<m_node_count;i++)
    if (m_mesh_in_contact[i]==0)
      trimesh->react_force[0].z+=cfn[i];
    
  #endif 
	//Correct time step!
//	std::min(deltat,dt_fext)
} //Contact Forces

////// BUG DETECTED

//~ void Domain_d::calcContactForceFromPressure(){
  
   //~ bool is_elem_sum[m_elem_count];
   //~ bool is_node_sum[m_node_count];
   //~ //double pxa_el[m_elem_count];
   //~ double cfsum=0.0;
   //~ double area = 0.0;
    //~ double cfnsum = 0.0;
    
   //~ for (int e=0;e<m_elem_count;e++)is_elem_sum[e]=false;
   
    //~ for (int i=0;i<m_node_count;i++){
      //~ if(m_mesh_in_contact[i]>-1 && m_mesh_in_contact[i]<1){
        //~ cfnsum += p_node[i]*node_area[i];
        //~ for (int ne=0; ne<m_nodel_count[i];ne++) {
          //~ int e   = m_nodel     [m_nodel_offset[i]+ne]; //Element
          //~ if (!is_elem_sum[e]){

            //~ //pxa_el[e]+=p[e]*m_elem_area[e];
            //~ is_elem_sum[e]=true;
            //~ //TODO: CRITICAL: ASSUMING ZZ
            //~ cfsum += m_sigma[6*e+2]*m_elem_area[e];
            //~ area+=m_elem_area[e];
          //~ }
            //~ if (!is_node_sum[i]){
              //~ area+=node_area[i];
              //~ }
        //~ }//nodel
      //~ }//mesh in contact
    //~ }
    
  
  //~ printf("Area %.3e\n", area);
  //~ trimesh->react_p_force[0] = cfsum;
  //~ trimesh->react_force[0].z = cfnsum;
   //~ trimesh->cont_area = area;
//~ }

void Domain_d::calcContactForceFromPressure(){
  
  double cfsum = 0.0;
  double area = 0.0;
  int facecount = 0;
  for (int i = 0; i < m_faceCount; i++) {
      if (faceList[i].count != 1) continue;

      int n0 = faceList[i].nodes[0];
      int n1 = faceList[i].nodes[1];
      int n2 = faceList[i].nodes[2];
      
      //APPROACH #1: Contact Force
      //~ // Verificar si los tres nodos tienen fuerza de contacto acumulada
      //~ if (norm(m_mesh_in_contact[n0]) == -1 ||
          //~ norm(m_mesh_in_contact[n1]) == -1 ||
          //~ norm(m_mesh_in_contact[n2]) == -1)
          //~ continue; // no hay contacto
      // Verificar si los tres nodos tienen fuerza de contacto acumulada
      if (m_mesh_in_contact[n0] == 0 &&
          m_mesh_in_contact[n1] == 0 &&
          m_mesh_in_contact[n2] == 0)
      {
          
        int e = faceList[i].elem_id;
        cfsum += m_sigma[6*e+2]*m_elem_area[e];
        area +=m_elem_area[e];
        facecount++;
        
        //~ // (opcional) Face normal
        //~ double3 p0 = getPosVec3(n0);
        //~ double3 p1 = getPosVec3(n1);
        //~ double3 p2 = getPosVec3(n2);

        //~ double3 v1 = p1 - p0;
        //~ double3 v2 = p2 - p0;
        //~ double3 face_normal = normalize(cross(v1, v2));

        //~ // (opcional) Promedio de normales de fuerza en nodos
        //~ double3 avg_contact_normal = normalize(
            //~ normalize(contact_force[n0]) +
            //~ normalize(contact_force[n1]) +
            //~ normalize(contact_force[n2])
        //~ );

        //~ // (opcional) Filter by angle
        //~ if (dot(face_normal, avg_contact_normal) < 0.8) continue;

        // Sumar fuerza: presión * área
      }
  }
  //printf("Contact Area %.3e\n", area);
  //printf("Contact Face Count: %d\n",facecount);
  trimesh->react_p_force[0] = cfsum;
  trimesh->cont_area = area;

}

}; //Namespace
//};//SPH
