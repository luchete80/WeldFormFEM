#include "Domain_d.h"
#include <iostream>
#include <vector>

// #include "tensor.cu"
#include "Matrix_temp.h"

#if CUDA_BUILD
#include "tensor.cuh"
#else
//#include "Matrix.h"
#endif

#include "tensor3.C"

using namespace std;

namespace MetFEM {
  
 

dev_t void Domain_d::calcElemStrainRates(){

  par_loop(e,m_elem_count){
    Matrix str_rate(3,3); //TODO: MAKE SYMM MATRIX
    Matrix rot_rate(3,3); //TODO: MAKE SYMM MATRIX);

//  if (e < m_elem_count) {

    for (int gp=0;gp<m_gp_count;gp++){
      int offset = e * m_gp_count * 6 + gp;
	    int offset_det = e * m_gp_count;
		  str_rate.SetZero();
      rot_rate.SetZero();
	//printf ("offset, %f , det %f\n", offset, m_detJ[offset + gp]);
      double f = 1.0 / m_detJ[offset_det + gp];
      //printf("det J: %f\n",m_detJ[offset_det + gp]);
      //double test = 0.0;
      for (int n=0; n<m_nodxelem;n++) {

        // double vele[3];
        // vele[0] = vele3.x;        vele[1] = vele3.y;        vele[2] = vele3.z;
        // do d=1, dim
          // !print *, "node dim dHxy vele", n,d,temp(d,n) , elem%vele (e,dim*(n-1)+d,1) 
          // elem%str_rate(e,gp, d,d) = elem%str_rate(e,gp, d,d) + temp(d,n) * elem%vele (e,dim*(n-1)+d,1) 
          // elem%rot_rate(e,gp, d,d) = 0.0d0
        // end do
        //test += getDerivative(e,gp,2,n) * f * getVElem(e,n,2);
        // printf("n %d deriv %f vele %f\n",n, getDerivative(e,gp,2,n),  getVElem(e,n,2));
        // printf ("Nod %d, vel %.6e  %.6e  %.6e \n", n, getVElem(e,n,0),getVElem(e,n,1),getVElem(e,n,2));
        for (int d=0;d<m_dim;d++){
          //printf("d %d n %d deriv %f vele %f\n",d, n, getDerivative(e,gp,d,n)*f,  getVElem(e,n,d));
          
          str_rate.Set(d,d, str_rate.getVal(d,d) + getDerivative(e,gp,d,n) * f * getVElem(e,n,d));
          rot_rate.Set(d,d, 0.0);

          // elem%str_rate(e,gp, d,d) = elem%str_rate(e,gp, d,d) + temp(d,n) * elem%vele (e,dim*(n-1)+d,1) 
          // elem%rot_rate(e,gp, d,d) = 0.0d0
          
        }//dim
        // !!!! TO AVOID ALL MATMULT
        str_rate.Set(0,1, str_rate.getVal(0,1) + f *(getDerivative(e,gp,1,n) * getVElem(e,n,0) +
                                                       getDerivative(e,gp,0,n) * getVElem(e,n,1)));
        rot_rate.Set(0,1, rot_rate.getVal(0,1) + f* (getDerivative(e,gp,1,n) * getVElem(e,n,0) -
                                                       getDerivative(e,gp,0,n) * getVElem(e,n,1)));
                                                       
        if (m_domtype == _Axi_Symm_){
          str_rate.Set(2,2,str_rate.getVal(0,1) + 0.25*getVElem(e,n,0)/getRadius(e,0));//0.25 is shapemat - Vr/r
          rot_rate.Set(2,2,0.0);
            // elem%str_rate(e,gp, 3,3) = elem%str_rate(e,gp, 3,3) + 0.25d0*elem%vele (e,dim*(n-1)+1,1) / elem%radius(e,gp) !! 0.25 is shapemat
          // ! print *, "hoop er", elem%str_rate(e,gp, 3,3) 
          // elem%rot_rate(e,gp, 3,3) = 0.0d0
        }
                                                       
        if (m_dim == 3) {
          //printf("elem %d velem %f %f %f\n", e, getVElem(e,n,1),getVElem(e,n,1),getVElem(e,n,2));
          //printf("deriv %f\n",getDerivative(e,gp,2,n));
          str_rate.Set(1,2, str_rate.getVal(1,2) + f *(getDerivative(e,gp,2,n) * getVElem(e,n,1) +
                                                         getDerivative(e,gp,1,n) * getVElem(e,n,2)));
          str_rate.Set(0,2, str_rate.getVal(0,2) + f *(getDerivative(e,gp,2,n) * getVElem(e,n,0) +
                                                         getDerivative(e,gp,0,n) * getVElem(e,n,2)));

          rot_rate.Set(1,2, rot_rate.getVal(1,2) + f *(getDerivative(e,gp,2,n) * getVElem(e,n,1) -
                                                         getDerivative(e,gp,1,n) * getVElem(e,n,2)));
          rot_rate.Set(0,2, rot_rate.getVal(0,2) + f *(getDerivative(e,gp,2,n) * getVElem(e,n,0) -
                                                         getDerivative(e,gp,0,n) * getVElem(e,n,2)));
        

        }// end if     
      }// end do !Nod x elem
      //printf ("test %fn", test);
      str_rate.Set(0,1, str_rate.getVal(0,1) *0.5);
      str_rate.Set(0,2, str_rate.getVal(0,2) *0.5);  
      str_rate.Set(1,2, str_rate.getVal(1,2) *0.5);
      
      rot_rate.Set(0,1, rot_rate.getVal(0,1) *0.5);
      rot_rate.Set(0,2, rot_rate.getVal(0,2) *0.5);  
      rot_rate.Set(1,2, rot_rate.getVal(1,2) *0.5);
      
      str_rate.Set(1,0,  str_rate.getVal(0,1));   
      str_rate.Set(2,1,  str_rate.getVal(1,2));  
      str_rate.Set(2,0,  str_rate.getVal(0,2));

      rot_rate.Set(1,0, -rot_rate.getVal(0,1));      
      rot_rate.Set(2,1, -rot_rate.getVal(1,2));  
      rot_rate.Set(2,0, -rot_rate.getVal(0,2)); 

      str_rate.ToFlatSymPtr(m_str_rate, offset);
      rot_rate.ToFlatSymPtr(m_rot_rate, offset); //UPPER PART

       //printf("Strain Rate\n");
       //str_rate.Print();

      // printf("Rot Rate\n");
      // str_rate.Print();
      
      // !elem%str_rate(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%vele (e,:,:)) 
      // !print *, "simlpified strain rate "

      //Inverse test
      // Matrix *test = new Matrix(3,3);
      // Matrix *invtest = new Matrix(3,3);
      // //printf("A\n");
      // test.Set(0,0,1);test.Set(0,1,1);test.Set(0,2,1);
      // test.Set(1,0,1);test.Set(1,1,2);test.Set(1,2,2);      
      // test.Set(2,0,1);test.Set(2,1,2);test.Set(2,2,3);
      // InvMat(*test,invtest);
      // ////printf("inv A\n");
      // test.Print();
      // invtest.Print();
        // delete test, invtest;
      } // Gauss Point

    }//if e<elem_count
    

} //calcElemStrains
  


////// STRAIN RATE CALC: SLOWER VARIANT 

// dev_t void Domain_d::calcElemStrainRates(){
    // double rot_rate[1][3][3];
    // double str_rate[1][3][3];
    // double tau[1][3][3];
    // double grad_v[m_nodxelem][m_dim][m_dim];

    // par_loop(e,m_elem_count){
    // Matrix str_rate_(3,3);
    // Matrix rot_rate_(3,3);    
    // for (int gp = 0; gp < m_gp_count; gp++) {
        // for (int I = 0; I < m_dim; I++) {
            // for (int J = 0; J < m_dim; J++){ 
                // grad_v[gp][I][J] = 0.0;
                // for (int k = 0; k < m_nodxelem; k++) {
                    // //grad_v[gp][I][J] += dNdX[gp][J][k] * vel[k][I];
                    // grad_v[gp][I][J] += getDerivative(0,gp,J,k) * getVElem(e,k,I)/m_detJ[e*m_gp_count+gp];
                    // //printf ("deriv %e " , getDerivative(0,gp,J,k)/m_detJ[gp]);
                    // //printf ("elem %d node %d, velem %f\n", e, k, getVElem(e,k,I));
                // }

            // }
        // }
    // }
    // cout << "Elem "<<e<<endl;
    // for (int gp = 0; gp < m_gp_count; gp++) {
      // int offset = e * m_gp_count * 6 + gp;
        // for (int i = 0; i < m_dim; i++) {
            // for (int j = 0; j < m_dim; j++) {
                // str_rate_.Set(i,j, 0.5*(grad_v[gp][i][j] + grad_v[gp][j][i]));
                // str_rate[gp][i][j] = 0.5 * (grad_v[gp][i][j] + grad_v[gp][j][i]);
                // rot_rate[gp][i][j] = 0.5 * (grad_v[gp][i][j] - grad_v[gp][j][i]);
                // //printf("str rate %e", str_rate[gp][i][j]);
            // }
        // }
        // // str_rate[gp][2][0]=rot_rate[gp][2][0]=0.0;                str_rate[gp][0][2]=rot_rate[gp][0][2]=0.0;        
        // // str_rate[gp][2][2]=rot_rate[gp][2][2]=0.0;
        // str_rate_.ToFlatSymPtr(m_str_rate, offset);
        // rot_rate_.ToFlatSymPtr(m_rot_rate, offset);
      // printf("Element %e Strain Rate\n", e);
      // str_rate_.Print();
      
    // }

    // // int o = 6*e;
    // // m_str_rate[o+0]=str_rate[0][0][0];     m_str_rate[o+1]=str_rate[0][1][1];    m_str_rate[o+2]=str_rate[0][2][2];  
    // // m_str_rate[o+3]=str_rate[0][0][1];     m_str_rate[o+4]=str_rate[0][1][2];    m_str_rate[o+5]=str_rate[0][0][2];  

    // // m_rot_rate[o+0]=rot_rate[0][0][0];     m_rot_rate[o+1]=rot_rate[0][1][1];    m_rot_rate[o+2]=rot_rate[0][2][2];  
    // // m_rot_rate[o+3]=rot_rate[0][0][1];     m_rot_rate[o+4]=rot_rate[0][1][2];    m_rot_rate[o+5]=rot_rate[0][0][2];  
    
    // } //Elem e

// }

__global__ void calcElemStrainRatesKernel(Domain_d *dom_d){
		
		dom_d->calcElemStrainRates();
}

  __global__ void calcElemStrainsKernel(Domain_d *dom_d){
		
		dom_d->calcElemStrainRates();
}

//To calculate before elemet Jacobian calc
dev_t void Domain_d::CalcElemVol(){
  par_loop(e,m_elem_count){
    double w;
    //TODO: CHANGE WEIGHT TO ARRAY
    if (m_gp_count == 1) {
      if (m_dim == 2) w = 4;//w = pow(2.0, m_dim);
      if (m_dim == 3)     
        if      (m_nodxelem == 4)  w = 1.0/6.0;
        else if (m_nodxelem == 8)  w = 8.0;
    } else                  w = 1.0;
    
    int offset = m_gp_count * e;
    vol[e] = 0.0;
    for (int gp=0;gp<m_gp_count;gp++){
      vol[e] += m_detJ[offset] * w;
    }  
  
    //if (e<10)
    //printf("Element %d Vol %f, det %f\n",e,vol[e],m_detJ[offset]);  
  }//el
}

dev_t void Domain_d::calcElemDensity(){ 
  par_loop(e,m_elem_count){
//!!!!! ASSUME VOLUME IS ALREADY CALCULATED
//subroutine calc_elem_density ()
  // implicit none
  // real(fp_kind), dimension(dim,dim) :: F !def gradient
  // real(fp_kind), dimension(nodxelem,dim) :: x !def gradient
  
  // integer :: e, n, gp
    for (int gp=0;gp<m_gp_count;gp++){
      int offset = e * m_gp_count ;
      rho[offset + gp] = rho_0[offset + gp] * vol_0[e] / vol[e];
      //printf("rho %.6e", rho[offset + gp]);    
    }

  // do e = 1, elem_count
    // do gp=1, elem%gausspc(e)
    // !if (elem%gausspc(e) .eq. 1) then
      // elem%rho(e,gp) = elem%rho_0(e,gp)*elem%vol_0(e)/elem%vol(e) !IS THE SAME
    // end do
  // end do
// end subroutine
  
  }
}

dev_t void Domain_d::calcAccel(){
  par_loop(n, m_node_count){
    int i = n*m_dim;
    for (int d=0;d<m_dim;d++){
      a[i+d] = (m_fe[i+d]-m_fi[i+d])/m_mdiag[n]; //TODO: REMAIN EXTERNAL FORCES
    }
    // if (n ==0){
      // printf("accel node 40 %f %f %f\n",a[i],a[i+1],a[i+2]);
    // }
    if (contact){
    for (int d=0;d<m_dim;d++){
      //printf("ADDED ACCEL: %f\n", contforce[i+d]);
      a[i+d] += contforce[i+d]/m_mdiag[n]; //TODO: REMAIN EXTERNAL FORCES
    }      
    }//contact
    //printf("mass %f\n",m_mdiag[n]);
    //printf("a %f %f %f \n",a[0],a[1],a[2]);
    //printf("f %f %f %f \n",m_fi[0],m_fi[1],m_fi[2]);
  }    
  
}

dev_t void Domain_d::CalcElemInitialVol(){
  CalcElemVol();
  par_loop(e,m_elem_count){
    int offset = m_gp_count * e;
    for (int gp=0;gp<m_gp_count;gp++){
      vol_0[e] = vol[e];
    }      
  }
  #ifndef CUDA_BUILD
  double vv = 0.0; 
  for (int e=0;e<m_elem_count;e++)
    vv+= vol_0[e];
  printf("Initial Vol: %.4e\n",vv);
  #endif
}

dev_t void Domain_d::calcTotMass(){
  // par_loop(e,m_elem_count){
// //!!!!! ASSUME VOLUME IS ALREADY CALCULATED
// //subroutine calc_elem_density ()
  // // implicit none
  // // real(fp_kind), dimension(dim,dim) :: F !def gradient
  // // real(fp_kind), dimension(nodxelem,dim) :: x !def gradient
  // int offset = e * m_gp_count ;
  // // integer :: e, n, gp
    // for (int gp=0;gp<m_gp_count;gp++){
      // rho[offset + gp] = rho[offset + gp] * vol_0[e] / vol[e];
    // }    

  // }  
}

dev_t void Domain_d::calcElemForces(){

  par_loop(e,m_elem_count){
    int offset = e*m_nodxelem*m_dim;
    double w = 1.;
    if (m_gp_count == 1) w = pow(2.0,m_dim);

    
    for (int n=0; n<m_nodxelem;n++) 
      for (int d=0;d<m_dim;d++)
        m_f_elem[offset + n*m_dim + d] = 0.0;

    int offset_det = m_gp_count * e;
    
    for (int gp=0;gp<m_gp_count;gp++){
      
      // tensor3 sigma     = FromFlatSym(m_sigma, e*m_gp_count+gp);
  // // integer :: e, i,j,k, gp,n, d
      // //printf("SIGMA\n");print(sigma);
  // real(fp_kind), dimension(dim*nodxelem,1) ::f
  // real(fp_kind) :: w
  // !TESTING
  // real (fp_kind) :: sigma_test(6,1) !ORDERED
  // real(fp_kind) :: test(24,1) !ifwanted to test in tensor form
  // elem%f_int = 0.0d0
  // w = 1.0d0 !!! Full integration
	
	// ! !$omp parallel do num_threads(Nproc) private (e,gp,d, w,n) 
  // do e=1, elem_count
    // if (elem%gausspc(e) .eq. 1) then
      // w = 2.0d0**dim
    // end if
    // do gp = 1, elem%gausspc(e)
      // !print *, "elem%dHxy_detJ(e,gp,1", elem%dHxy_detJ(e,gp,1,:)
      // !print *, "elem%dHxy_detJ(e,gp,2", elem%dHxy_detJ(e,gp,2,:)
      // ! sigma_test (:,1)=[elem%sigma (e,gp, 1,1),elem%sigma (e,gp, 2,2),elem%sigma (e,gp, 3,3),&
                        // ! elem%sigma (e,gp, 1,2),elem%sigma (e,gp, 2,3),elem%sigma (e,gp, 3,1)]
      // ! test = w*matmul(transpose(elem%bl(e,gp,:,:)),sigma_test)  ! (24x6)(6x1)
      // !print *, "test force", test
      
      // !print *, "dHdxy, 1", elem%dHxy_detJ(e,gp,1,:)
      // !print *, "dHdxy, 2", elem%dHxy_detJ(e,gp,2,:)
      // !print *, "dHdxy, 3", elem%dHxy_detJ(e,gp,1,:)
      
      
      // do n=1, nodxelem
      // !Is only linear matrix?    
      // !elem%f_int(e,n,d) =  
      // !f (:,:) = matmul(transpose(elem%bl(e,gp,:,:)),elem%sigma (e,:,:))
      // !!!! TO AVOID MATRIX MULTIPLICATIONS (8x6 = 48 in bathe notation with several nonzeros)
      // !!!!! F = BT x sigma = [dh1/dx dh1/dy ] x [ sxx sxy]
      // !!!!!                = [dh2/dx dh2/dy ]   [ syx syy]
      // !!!!! 
      //for (int i=0;i<3;i++)
      //  for (int j=0;j<3;j++){
          //printf("SIGMA %d %d %.6e\n",i,j,getSigma(e,gp,i,j));
      //  }
      for (int n=0; n<m_nodxelem;n++) {
        for (int d=0;d<m_dim;d++){
          m_f_elem[offset + n*m_dim + d] += getDerivative(e,gp,d,n) * getSigma(e,gp,d,d);
        }
      
        if (m_dim == 2){
          if (m_domtype != _Axi_Symm_){
          m_f_elem[offset + n*m_dim    ] +=  getDerivative(e,gp,1,n) * getSigma(e,gp,0,1);
          m_f_elem[offset + n*m_dim + 1] +=  getDerivative(e,gp,0,n) * getSigma(e,gp,0,1);
          } else {//2D AXISYMM, VOLUMETIRC WEIGHT (ALWAYS)
                if (m_axisymm_vol_weight) {
                    double r_gp = m_radius[e];           // r en el punto de integración
                    double detJ_gp = m_detJ[offset_det + gp];
                    
                    double sigma_rr  = getSigma(e, gp, 0, 0);
                    double sigma_zz  = getSigma(e, gp, 1, 1);
                    double sigma_tt  = getSigma(e, gp, 2, 2);  // sigma_theta_theta
                    double sigma_rz  = getSigma(e, gp, 0, 1);

                    // Fuerza en r
                    m_f_elem[offset + n*m_dim    ] += 
                        getDerivative(e, gp, 1, n) * sigma_rz * r_gp + 
                        0.25 * (sigma_rr - sigma_tt) * detJ_gp;

                    // Fuerza en z
                    m_f_elem[offset + n*m_dim + 1] += 
                        getDerivative(e, gp, 0, n) * sigma_rz * r_gp + 
                        0.25 * sigma_rz * detJ_gp;
                }else {
                    double r_gp = m_radius[e];
                    double detJ_gp = m_detJ[offset_det + gp];
                    double fa = 0.25 / r_gp * detJ_gp;

                    double sigma_rr  = getSigma(e, gp, 0, 0);
                    double sigma_tt  = getSigma(e, gp, 2, 2);
                    double sigma_rz  = getSigma(e, gp, 0, 1);

                    // Fuerza en r
                    m_f_elem[offset + n*m_dim    ] += 
                        getDerivative(e, gp, 1, n) * sigma_rz - 
                        (sigma_rr - sigma_tt) * fa;

                    // Fuerza en z
                    m_f_elem[offset + n*m_dim + 1] += 
                        getDerivative(e, gp, 0, n) * sigma_rz - 
                        sigma_rz * fa;
                }
            
          }
        } else { //3D
          //printf("offset %d\n", offset + n*m_dim    );
          //printf ("sigma 0 1 %f\n", getSigma(e,gp,0,1));
          m_f_elem[offset + n*m_dim    ] +=  getDerivative(e,gp,1,n) * getSigma(e,gp,0,1) +
                                             getDerivative(e,gp,2,n) * getSigma(e,gp,0,2);
          m_f_elem[offset + n*m_dim + 1] +=  getDerivative(e,gp,0,n) * getSigma(e,gp,0,1) + 
                                             getDerivative(e,gp,2,n) * getSigma(e,gp,1,2);        
          m_f_elem[offset + n*m_dim + 2] +=  getDerivative(e,gp,1,n) * getSigma(e,gp,1,2) + 
                                             getDerivative(e,gp,0,n) * getSigma(e,gp,0,2);     
        }

      
 
        
      }// nod x elem


    } // Gauss Point
    
    for (int n=0; n<m_nodxelem;n++) {
      for (int d=0;d<m_dim;d++){
        m_f_elem[offset + n*m_dim + d] *= w;
      }
    }  
    
    //TO CHECK
    
    //for (int n=0; n<m_nodxelem;n++) {
    //  printf("Element %d Node %d forces %.3e %.3e %.3e\n",e, n, m_f_elem[offset + n*m_dim ],m_f_elem[offset + n*m_dim + 1] ,m_f_elem[offset + n*m_dim + 2] );
    //} 
     
    
  }//if e<elem_count
}

// ORIGINAL
dev_t void Domain_d::calcElemPressure_Hex(){

  //~ par_loop(e,m_elem_count){
    //~ //printf("calc pressure \n");
    //~ int offset_t = e * m_gp_count *6;

     //~ double trace;
     //~ double press_inc = 0.0;
     //~ for (int gp=0;gp<m_gp_count;gp++){
       //~ trace = 0.0;
       //~ tensor3 str_inc     = FromFlatSym(m_str_rate,     offset_t +gp)*dt;
       //~ //printf("str inc, dt %f\n", dt);print(str_inc);
       //~ press_inc += Trace(str_inc);
      //~ }//gauss point
       //~ press_inc = -press_inc/m_gp_count;
       //~ //  //printf("trace %f\n",trace);
       //~ int offset = e * m_gp_count ;
       //~ for (int gp=0;gp<m_gp_count;gp++){
        //~ //  //printf("bulk mod:%f, press inc%f\n", mat[e]->Elastic().BulkMod(),press_inc);
       //~ trace = 0.0;
       //~ for (int d = 0; d<3;d++) trace += getSigma(e,gp,d,d);
        
       //~ p[offset + gp] = -1.0/3.0 * trace + mat[e]->Elastic().BulkMod() * press_inc;
       //~ //printf("pressure %f\n",p[offset + gp]);
     //~ }

   //~ } // e< elem_count
}

////ALT PRESSURE CALC (Previous, Low Press, no mix Nodal/Loc FBar)
//~ dev_t void Domain_d::calcElemPressure() {
  //~ // Compute nodal volumes (reused for hourglass control)
  //~ double *voln_0 = new double[m_node_count];
  //~ double *voln = new double[m_node_count];
  
  //~ par_loop(n, m_node_count) {
    //~ voln_0[n] = voln[n] = 0.0;
    //~ for (int i = 0; i < m_nodel_count[n]; ++i) {
      //~ int e = m_nodel[m_nodel_offset[n] + i];
      //~ voln_0[n] += vol_0[e];  // Sum elemental ref volumes
      //~ voln[n]   += vol[e];    // Sum elemental current volumes
    //~ }
  //~ }

  //~ for (int e=0;e<m_elem_count;e++){
  //~ //par_loop(e, m_elem_count) {
    //~ double K = mat[e]->Elastic().BulkMod();
    //~ double rho_e = rho[e];  // Densidad actual del elemento
    //~ double vol0 = vol_0[e];
    //~ double vol1 = vol[e];
    //~ double J =vol1/vol0;


    //~ // 1. Div Velocity 
    //~ double div_v = 0.0;

    //~ for (int a = 0; a < m_nodxelem; ++a) {
      //~ int nid = m_elnod[e * m_nodxelem + a];
      //~ double3 va = getVelVec(nid); // Velocidad nodal
      //~ double3 gradNa =make_double3(getDerivative(e,0,0,a),getDerivative(e,0,1,a),getDerivative(e,0,2,a));
      //~ div_v += dot(gradNa, va);
    //~ }

    //~ // 2. PSPG: Estabilización de presión
    //~ double h = pow(vol[e], 1.0/3.0); // Longitud característica
    //~ double mu_eff = 0.1 * sigma_y[e]; // Ej: σ_y = 100e6 Pa → mu_eff = 10e6 Pa·s
    //~ double tau = (h*h) / (4.0 * mu_eff);
    //~ double p_pspg = tau * div_v; // Término clave!

//~ // 3. F-bar: Corregir J (evitar locking)
    //~ double J_physical = vol[e] / vol_0[e]; // J sin corregir
    //~ double J_avg = 0.0;
    //~ for (int a = 0; a < m_nodxelem; ++a) {
      //~ int nid = m_elnod[e * m_nodxelem + a];
      //~ J_avg += (voln[nid] / voln_0[nid]); // J nodal promediado
    //~ }
    //~ J_avg /= m_nodxelem;
    //~ double J_bar = J_avg; // Usar J promedio para el elemento (F-bar volumétrico)
    

    //~ // 4. Presión física CORREGIDA (sin ANP!)
    //~ double p_physical = -mat[e]->Elastic().BulkMod() * log(J_bar);
    
    //~ //Artif Visc
    //~ double q = 0.0;
    //~ if (div_v < 0.0) {
      //~ double c = sqrt(K / rho_e); // Velocidad del sonido
      //~ //double h = elem_char_length[e]; // Longitud característica del elemento
      //~ //double h = getMinLength();
      //~ double h = pow(vol[e], 1.0/3.0); 
      //~ double alpha = 0.12;
      //~ double beta = 0.01;
      

      //~ q = rho_e * (-alpha * c * h * div_v + beta * h * h * div_v * div_v);
    //~ }
    
    //~ // Hourglass volumétrico (opcional, ajustable)


    //~ double hg_coeff = 0.0;
    //~ double p_hg = hg_coeff * K * (J - J_avg);
    
    //~ // Presión final
    //~ //p[e] = p_vol + q;

    //~ // 5. Combinar con PSPG
    //~ p[e] = p_physical + p_pspg + q + p_hg; // PSPG suprime oscilaciones
    
  //~ }
  //~ delete[] voln_0;
  //~ delete[] voln;
  
//~ }

dev_t void Domain_d::calcElemPressure() {
  
  // 1. Calcular volúmenes nodales acumulados
  double *voln_0 = new double[m_node_count];
  double *voln   = new double[m_node_count];

 // 1. Compute nodal volumes - parallel optimized
  par_loop(n, m_node_count) {
    voln_0[n] = voln[n] = 0.0;
    for (int i = 0; i < m_nodel_count[n]; ++i) {
      int e = m_nodel[m_nodel_offset[n] + i];
      voln_0[n] += vol_0[e];
      voln[n]   += vol[e];
    }
  }

    // 2. Material/Stabilization parameters
    //~ const double alpha_contact = 0.5;  
    //~ const double alpha_free = 0.2;      
    //~ double hg_coeff_free = 0.2;    // 
    //~ double hg_coeff_contact = 0.1;  // Slightly lower in contact
    //~ const double artvisc_coeff = 0.15;    // Artificial viscosity
    //~ const double log_factor = 0.8;
    //~ double pspg_scale = 0.5;  // Escalar a 10% del valor original
    //~ const double p_pspg_bulkfac = 0.1;
    
        
    // 3. Element loop - main computation
  par_loop(e, m_elem_count) {
    double K  = mat[e]->Elastic().BulkMod();
    double mu = mat[e]->Elastic().G();
    double rho_e = rho[e];
    double vol0 = vol_0[e];
    double vol1 = vol[e];
    double J_local = vol1 / vol0;
    double h = pow(vol1, 1.0/3.0);
    
    // Contact detection
    double contact_weight = 0.0;
    bool is_contact = false;
    for(int a = 0; a < m_nodxelem; ++a) {
      int nid = m_elnod[e*m_nodxelem + a];
      double3 cf = make_double3(contforce[m_dim*nid],contforce[m_dim*nid+1],contforce[m_dim*nid+2]);
        if(dot(cf,cf) > 0) {
            is_contact = true;
            break;
        }
        //contact_weight = std::max(contact_weight, std::min(1.0, length(cf) / (K * h * h)));
    }

    // F-bar adaptive blending
    double J_avg = 0.0;
    for(int a = 0; a < m_nodxelem; ++a) {
        int nid = m_elnod[e*m_nodxelem + a];
        J_avg += voln[nid] / voln_0[nid];
    }
    J_avg /= m_nodxelem;

    // Critical: Contact-adaptive blending
    //double alpha = is_contact ? 0.85 : 0.4;  // More local in contact
    double alpha = is_contact ? m_stab.alpha_contact : m_stab.alpha_free;  // More local in contact
    //double alpha = alpha_free + (alpha_contact - alpha_free) * contact_weight;
    double J_bar = alpha*J_local + (1-alpha)*J_avg;
     // IMPROVED PHYSICAL PRESSURE (Hybrid model)
    //double p_physical = -K * (m_stab.log_factor*log(J_bar) + (1.0-m_stab.log_factor)*(J_bar - 1.0));

    //double p_physical = -K * log(J_bar);

    double p_physical = -K * (log(J_bar) + (mu / K) * (J_bar * J_bar - 1.0));

    // Enhanced PSPG - dynamic tau calculation
    double c = sqrt(K / rho_e);  // Sound speed
    double tau = h / (2.0 * c);  // Dynamic stabilization
    
    // Velocity divergence
    double div_v = 0.0;
    for(int a = 0; a < m_nodxelem; ++a) {
        int nid = m_elnod[e*m_nodxelem + a];
        double3 va = getVelVec(nid);
        double3 gradNa =make_double3(getDerivative(e,0,0,a),getDerivative(e,0,1,a),getDerivative(e,0,2,a));
        div_v += dot(gradNa, va);
    }
    
    //div_v /= vol1;  // Normalización
    
    double p_pspg = 0.0;  // LIMITED TO COMPRESSION 
    //double p_pspg = std::min(pspg_scale * tau * div_v, p_pspg_bulkfac * K);  // Límite del 5% de 
            
    // Non-negative hourglass (stabilization only)
    double p_hg = (is_contact ? m_stab.hg_coeff_contact : m_stab.hg_coeff_free) * K * fabs(J_local - J_avg);
    
    // Artificial viscosity - compression only
    double p_q = 0.0;
    if(div_v < 0.0) {
        p_pspg = std::min(m_stab.pspg_scale * tau * div_v *K, m_stab.p_pspg_bulkfac * K);  // Límite del 5% de 
        //cout << "pspg "<<pspg_scale * tau * div_v *K<<endl;
        double q1 = m_stab.av_coeff_div * rho_e * h * c * (-div_v);
        double delta_J = 1.0 - J_local;
        double q2 = m_stab.av_coeff_bulk * K * delta_J;  // Volumetric term
        if (is_contact) {
            p_q = 0.5 * (q1 + q2);  // Mezcla en contacto
        } else {
            p_q = std::max(q1, q2);
        }
    }

    // Contact pressure boost (additional 10-15% in contact zones)
    
    // FINAL PRESSURE (contact boosted)
    p[e] = p_physical + p_pspg + p_hg + p_q;
    //p[e] = std::max(-10.0 * K, std::min(2.0 * K, p[e]));  // Límites conservadores
    
  }

  delete[] voln_0;
  delete[] voln;
}


void Domain_d::smoothPressureField(double gamma) {
  std::vector<double> p_new(m_elem_count, 0.0);
  
  for (int e = 0; e < m_elem_count; ++e) {
    double sum_p = 0.0;
    int count = 0;
    //printf("element nb%d\n",m_elem_neigh_count[e]);
    for (int i = 0; i < m_elem_neigh_count[e]; ++i) {
      int e_neigh = m_elem_neigh[4*e+i];
      sum_p += p[e_neigh] - p[e];
      count++;
    }
    if (count > 0) {
      p_new[e] = p[e] + gamma * sum_p / count;
    } else {
      p_new[e] = p[e];
    }
  }

  for (int e = 0; e < m_elem_count; ++e) p[e] = p_new[e];
}

dev_t void Domain_d::calcElemPressure_Hybrid() {
  // Compute nodal volumes (reused for hourglass control)
  double *voln_0 = new double[m_node_count];
  double *voln = new double[m_node_count];
  
  par_loop(n, m_node_count) {
    voln_0[n] = voln[n] = 0.0;
    for (int i = 0; i < m_nodel_count[n]; ++i) {
      int e = m_nodel[m_nodel_offset[n] + i];
      voln_0[n] += vol_0[e];  // Sum elemental ref volumes
      voln[n]   += vol[e];    // Sum elemental current volumes
    }
  }

  // Blended elemental-nodal pressure
  par_loop(e, m_elem_count) {
    double K = mat[e]->Elastic().BulkMod();
    double J_elem = vol[e] / vol_0[e];  // Elemental J
    double p_elem = K * (1.0 - J_elem); // Elemental pressure

    // Nodal J average (for stabilization)
    double J_nod = 0.0;
    for (int a = 0; a < m_nodxelem; ++a) {
      int nid = m_elnod[e * m_nodxelem + a];
      J_nod += voln[nid] / voln_0[nid]; 
    }
    J_nod /= m_nodxelem;
    double p_nod = K * (1.0 - J_nod);  // Nodal pressure

    // Blend 80% elemental + 20% nodal (adjust weights as needed)
    p[e] = 0.7 * p_elem + 0.3 * p_nod; 
  }

  delete[] voln_0;
  delete[] voln;
}

dev_t void Domain_d::calcElemPressure_Hybrid_VolHG() {
  // ---- Step 1: Compute Nodal Volumes (for HG control) ----
  double *voln_0 = new double[m_node_count];
  double *voln = new double[m_node_count];
  
  par_loop(n, m_node_count) {
    voln_0[n] = voln[n] = 0.0;
    for (int i = 0; i < m_nodel_count[n]; ++i) {
      int e = m_nodel[m_nodel_offset[n] + i];
      voln_0[n] += vol_0[e];  // Sum elemental ref volumes
      voln[n]   += vol[e];    // Sum elemental current volumes
    }
  }

  // ---- Step 2: Hybrid Pressure + Volumetric HG ----
  par_loop(e, m_elem_count) {
    double K = mat[e]->Elastic().BulkMod();
    double J_elem = vol[e] / vol_0[e];  // Elemental Jacobian
    double p_elem = K * (1.0 - J_elem); // Elemental pressure (hyperelastic)

    // Compute nodal J average for HG stabilization
    double J_nod = 0.0;
    for (int a = 0; a < m_nodxelem; ++a) {
      int nid = m_elnod[e * m_nodxelem + a];
      J_nod += voln[nid] / voln_0[nid]; 
    }
    J_nod /= m_nodxelem;

    // ---- Volumetric Hourglass Term (Critical for Stability) ----
    double hg_coeff = 0.0;  // Adjusted for metal forming (typical range: 0.02–0.05)
    double p_hg = hg_coeff * K * (J_elem - J_nod);  // HG pressure correction

    // ---- Final Blended Pressure ----
    p[e] = 0.7 * p_elem + 0.3 * (K * (1.0 - J_nod)) + p_hg;  // 80% elemental, 20% nodal + HG
  }

  delete[] voln_0;
  delete[] voln;
}

//Computational Methods in lagfangian & Eulerian Hydrocodes
//From Benson 1992
// Equation 1.3.12
//~ dev_t void Domain_d::calcElemPressureFromJ(){

  //~ par_loop(e,m_elem_count){
    //~ p[e] = mat[e]->Elastic().BulkMod() * ( 1.0 - vol[e]/vol_0[e] );
  //~ } // e< elem_count
//~ }


//~ dev_t void Domain_d::calcNodalPressureFromElemental() {
  //~ par_loop(n, m_node_count) {
    //~ p_node[n] = 0.0;
    //~ int count = 0;
    //~ for (int i = 0; i < m_nodel_count[n]; ++i) {
      //~ int e = m_nodel[m_nodel_offset[n] + i];
      //~ p_node[n] += p[e];
      //~ count++;
    //~ }
    //~ if (count > 0)
      //~ p_node[n] /= count;
  //~ }
//~ }

dev_t void Domain_d::calcNodalPressureFromElemental() {
  
  double* vol_acc = new double[m_node_count];

  par_loop(n, m_node_count) {
      p_node[n] = 0.0;
      vol_acc[n] = 0.0;
  }

//  par_loop(e, m_elem_count) {
  for (int e=0;e<m_elem_count;e++){
      for (int a = 0; a < m_nodxelem; ++a) {
          int nid = m_elnod[e*m_nodxelem + a];
          p_node[nid] += p[e] * vol[e];
          vol_acc[nid] += vol[e];
          //atomicAdd(&p_node[nid], p[e] * vol[e]);
          //atomicAdd(&vol_acc[nid], vol[e]);
      }
  }

  par_loop(n, m_node_count) {
      if (vol_acc[n] > 0.0)
          p_node[n] /= vol_acc[n];
  }
  
  delete[] vol_acc;
}

//FROM BENSON 1998
// SHOULDBE CALCULATE INITIAL AND CURRENT VOLUMES FIRST 
//ONLY FOR CONSTANT TETRA
//JOLDES, WITTEK, MILLER
//Non-locking tetrahedral finite element for surgical simulation
dev_t void Domain_d::calcElemPressureANP(){
  double *pn = new double [m_node_count];
  double *voln_0 = new double [m_node_count];
  double *voln = new double [m_node_count];

   //assume same bulk modulus
  double k = mat[0]->Elastic().BulkMod();

   par_loop(n, m_node_count){
     voln_0[n]=0.0;
     voln  [n]=0.0;
     pn[n]    = 0.0;
     for (int e=0; e<m_nodel_count[n];e++) {    
       int eglob   = m_nodel     [m_nodel_offset[n]+e]; //Element
       voln_0[n]+= /*0.25**/vol_0[eglob];
       voln[n] +=/*0.25**/vol[eglob]; 
     }
     //printf("Node %d vol %f\n", n,voln[n]);
     pn[n] = k*(1.0 - voln[n]/voln_0[n]); //0.25 is not necesary since is dividing
     p_node[n] = pn[n];
   } //NODE LOOP

   par_loop(e,m_elem_count){
     for (int ne=0;ne<m_nodxelem;ne++) //I.E. 4 
       p[e] += pn[m_elnod[e * m_nodxelem+ne]];    
     p[e] *= 0.25;
   }
   delete []pn;
   delete []voln_0;
   delete []voln;
 }

////CORRECTED
dev_t void Domain_d::calcElemPressureANP_Nodal() {
  double *pn = new double[m_node_count];
  double *voln_0 = new double[m_node_count];
  double *voln = new double[m_node_count];

  // Asumimos mismo bulk modulus para todos los elementos
  double k = mat[0]->Elastic().BulkMod();

  // Inicializar volúmenes nodales
  par_loop(n, m_node_count) {
    voln_0[n] = 0.0;
    voln[n]   = 0.0;
    pn[n]     = 0.0;

    for (int i = 0; i < m_nodel_count[n]; ++i) {
      int e = m_nodel[m_nodel_offset[n] + i];

      // Distribuir el volumen entre los nodos (tetra = 4 nodos)
      voln_0[n] += vol_0[e] / 4.0;
      voln[n]   += vol[e]   / 4.0;
    }

    // Cálculo del J nodal y presión nodal
    if (voln_0[n] > 1e-12) {
      double Jn = voln[n] / voln_0[n];
      pn[n] = k * (1.0 - Jn);
    } else {
      pn[n] = 0.0; // O alguna penalización por nodo inválido
    }

    p_node[n] = pn[n]; // Guardar para visualización
  }

  // Promediar presiones nodales para obtener presión elemental
  par_loop(e, m_elem_count) {
    p[e] = 0.0;
    for (int a = 0; a < m_nodxelem; ++a) {
      int nid = m_elnod[e * m_nodxelem + a];
      p[e] += pn[nid];
    }
    p[e] /= m_nodxelem;
  }

  delete[] pn;
  delete[] voln_0;
  delete[] voln;
}

dev_t void Domain_d::calcElemPressureANP_Nodal_HG() {
  double *pn = new double[m_node_count];
  double *voln_0 = new double[m_node_count];
  double *voln = new double[m_node_count];
  double *Jelem = new double[m_elem_count];  // Store elemental J for HG term

  // Assume uniform bulk modulus (or loop per element if needed)
  double k = mat[0]->Elastic().BulkMod();
  double hg_coeff = 0.02;  // Hourglass coefficient (tune as needed)

  // --- Step 1: Compute nodal volumes and pressures ---
  par_loop(n, m_node_count) {
    voln_0[n] = 0.0;
    voln[n] = 0.0;
    pn[n] = 0.0;

    for (int i = 0; i < m_nodel_count[n]; ++i) {
      int e = m_nodel[m_nodel_offset[n] + i];
      voln_0[n] += vol_0[e] / 4.0;  // Tet4 volume distribution
      voln[n] += vol[e] / 4.0;
    }

    // Nodal pressure (hyperelastic)
    if (voln_0[n] > 1e-12) {
      pn[n] = k * (1.0 - voln[n] / voln_0[n]);
    }
    p_node[n] = pn[n];  // For visualization
  }

  // --- Step 2: Compute elemental J (for HG term) ---
  par_loop(e, m_elem_count) {
    Jelem[e] = vol[e] / vol_0[e];  // Elemental Jacobian
  }

  // --- Step 3: Final elemental pressure (nodal avg + HG) ---
  par_loop(e, m_elem_count) {
    // Nodal-averaged J
    double Jnod = 0.0;
    for (int a = 0; a < m_nodxelem; ++a) {
      int nid = m_elnod[e * m_nodxelem + a];
      Jnod += voln[nid] / voln_0[nid];
    }
    Jnod /= m_nodxelem;

    // Hourglass stabilization term
    double p_hg = hg_coeff * k * (Jelem[e] - Jnod);

    // Final pressure: nodal avg + HG
    p[e] = 0.0;
    for (int a = 0; a < m_nodxelem; ++a) {
      int nid = m_elnod[e * m_nodxelem + a];
      p[e] += pn[nid];
    }
    p[e] = p[e] / m_nodxelem + p_hg;  // Add HG correction
  }

  delete[] pn;
  delete[] voln_0;
  delete[] voln;
  delete[] Jelem;
}

///// EXPERIMENTAL NEW
//~ dev_t void Domain_d::calcElemPressureANP_Element() {
  //~ double *voln_0 = new double[m_node_count];
  //~ double *voln = new double[m_node_count];

  //~ // Acumulación de volumen nodal (para calcular Jnod más adelante)
  //~ par_loop(n, m_node_count) {
    //~ voln_0[n] = 0.0;
    //~ voln[n] = 0.0;
    //~ for (int i = 0; i < m_nodel_count[n]; ++i) {
      //~ int e = m_nodel[m_nodel_offset[n] + i];
      //~ voln_0[n] += vol_0[e] / 4.0;
      //~ voln[n] += vol[e] / 4.0;
    //~ }
  //~ }

  //~ // Presión por elemento con estabilización tipo hourglass
  //~ par_loop(e, m_elem_count) {
    //~ double k = mat[e]->Elastic().BulkMod();

    //~ double J = vol[e] / vol_0[e];
    //~ double p_vol = k * (1.0 - J);

    //~ // Promedio nodal de J
    //~ double Jnod = 0.0;
    //~ for (int a = 0; a < m_nodxelem; ++a) {
      //~ int nid = m_elnod[e * m_nodxelem + a];
      //~ if (voln_0[nid] > 1e-12) {
        //~ Jnod += voln[nid] / voln_0[nid];
      //~ } else {
        //~ Jnod += 1.0; // Neutral (sin deformación)
      //~ }
    //~ }
    //~ Jnod /= m_nodxelem;

    //~ // Coeficiente de hourglass volumétrico (ajustable)
    //~ double hg_coeff = 0.05; // probar entre 0.01–0.1
    //~ double p_hg = hg_coeff * k * (J - Jnod);

    //~ p[e] = p_vol + p_hg;
  //~ }

  //~ // Presión nodal para postproceso
  //~ par_loop(n, m_node_count) {
    //~ p_node[n] = 0.0;
    //~ int count = 0;
    //~ for (int i = 0; i < m_nodel_count[n]; ++i) {
      //~ int e = m_nodel[m_nodel_offset[n] + i];
      //~ p_node[n] += p[e];
      //~ count++;
    //~ }
    //~ if (count > 0)
      //~ p_node[n] /= count;
  //~ }

  //~ delete[] voln_0;
  //~ delete[] voln;
//~ }

// Alternative: Use element-based approach (often more stable)
dev_t void Domain_d::calcElemPressureElementBased(){
  double k = mat[0]->Elastic().BulkMod() * 0.1; // Reduced bulk modulus
  
  par_loop(e, m_elem_count){
    if (vol_0[e] > 1e-15) {
      double vol_ratio = vol[e] / vol_0[e];
      vol_ratio = fmax(0.2, fmin(5.0, vol_ratio)); // Clamp ratio
      
      p[e] = k * (1.0 - vol_ratio);
      
      // Pressure limiting
      double max_pressure = 1.0 * k;
      p[e] = fmax(-max_pressure, fmin(max_pressure, p[e]));
    } else {
      p[e] = 0.0;
    }
  }
  
  // Update nodal pressures from element pressures
  par_loop(n, m_node_count){
    p_node[n] = 0.0;
    int elem_count = 0;
    
    for (int e=0; e<m_nodel_count[n]; e++) {    
      int eglob = m_nodel[m_nodel_offset[n]+e];
      p_node[n] += p[eglob];
      elem_count++;
    }
    
    if (elem_count > 0) {
      p_node[n] /= double(elem_count);
    }
  }
}


///// ASSUMING CONSTANT element node count
dev_t void Domain_d::CalcNodalVol(){
  double tot_vol = 0.0; //Only for verif
  par_loop(n, m_node_count){
    m_voln[n]=0.0;
    for (int e=0; e<m_nodel_count[n];e++) {    
      int eglob   = m_nodel     [m_nodel_offset[n]+e]; //Element
      //printf ("eglob %d, vol %e\n",eglob,vol[eglob]);
      m_voln[n] += /*1.0/m_nodxelem */ vol[eglob]; 
    }
    m_voln[n]/=m_nodxelem;
    //printf("Node %d vol %f ne count %d\n",n,m_voln[n],m_nodel_count[n]);
    tot_vol+=m_voln[n];
  } //NODE LOOP
  //printf("Total vol %f\n",tot_vol);
}

// //Assuming constant material
// dev_t void Domain_d::CalcNodalMassFromVol(){
  // double *rhon = new double [m_node_count];  

  
  // par_loop(n, m_node_count){
    // rhon[n]=0.0;
    // for (int e=0; e<m_nodel_count[n];e++) {    
      // int eglob   = m_nodel     [m_nodel_offset[n]+e]; //Element
      // rhon[n] += rho[eglob]; 
    // }
    // m_mdiag[n] = rhon[n]/(double)m_nodel_count[n] * m_voln[n];
    // //printf("Node %d mass %f rho %f vol %.4e\n",n,m_mdiag[n],rhon[n]/(double)m_nodel_count[n] , m_voln[n]);
    
  // } //NODE LOOP
  // double tot_mass = 0.0;
  // for (int n=0;n<m_node_count;n++)
    // tot_mass +=m_mdiag[n];
    
  // printf("Total Nodal Mass: %f\n",tot_mass);
  // delete rhon;
// }


//Assuming constant material
dev_t void Domain_d::CalcNodalMassFromVol(){
  double *rhon = new double [m_node_count];  

  
  par_loop(n, m_node_count){
    double mass = 0.0;
    double f = 1.0;

    for (int e=0; e<m_nodel_count[n];e++) {    
      int eglob   = m_nodel     [m_nodel_offset[n]+e]; //Element
      if (m_domtype == _Axi_Symm_){
        f = 1.0/m_radius[e];
      }
      mass += f*rho[eglob] * m_voln[n] / m_nodel_count[n];  //BENSON 1992.
    }    
    m_mdiag[n] = mass;
  } //NODE LOOP
    //printf("Node %d mass %f rho %f vol %.4e\n",n,m_mdiag[n],rhon[n]/(double)m_nodel_count[n] , m_voln[n]);

  double tot_mass = 0.0;
  for (int n=0;n<m_node_count;n++)
    tot_mass +=m_mdiag[n];
    
  printf("Total Nodal Mass: %f\n",tot_mass);
  delete rhon;
}

// subroutine Calc_Elastic_Stress(dom, dt)
dev_t void Domain_d::Calc_Elastic_Stress(const double dt){
  // integer :: e,gp
  // real(fp_kind), intent (in) :: dt
  // real(fp_kind) :: c
  // type (dom_type), intent (in) :: dom
  
  // !!!! PLAIN STRESS
  // !c = dom%mat_E / (1.0-dom%mat_nu*dom%mat_nu)
  
  // !!!! PLAIN STRAIN
  // c = dom%mat_E / ((1.0+dom%mat_nu)*(1.0-2.0*dom%mat_nu))
  // ! print *, "MAT C", c
  // ! print *, "MAT G", mat_G
  // do e = 1, elem_count 
    // do gp=1,elem%gausspc(e)
      // ! elem%str_inc(e,gp,:,:) = elem%str_inc(e,gp,:,:) + elem%str_inc(e,gp,:,:) * dt
      // ! elem%sigma(e,gp,1,1) = elem%sigma(e,gp,1,1) + c * (elem%str_inc(e,gp,1,1)+dom%mat_nu*elem%str_inc(e,gp,2,2))
      // ! elem%sigma(e,gp,2,2) = elem%sigma(e,gp,2,2) + c * (elem%str_inc(e,gp,2,2)+dom%mat_nu*elem%str_inc(e,gp,1,1))
      // ! elem%sigma(e,gp,1,2) = elem%sigma(e,gp,1,2) + 2.0* mat_G * elem%str_inc(e,gp,1,2)
      // ! elem%sigma(e,gp,2,1) = elem%sigma(e,gp,1,2)
      
      // !!!! PLAIN STRAIN  
      // elem%str_inc(e,gp,:,:) = elem%str_inc(e,gp,:,:) + elem%str_inc(e,gp,:,:) * dt
      // elem%sigma(e,gp,1,1) = elem%sigma(e,gp,1,1) + c * ((1.0-dom%mat_nu)*elem%str_inc(e,gp,1,1)+dom%mat_nu*elem%str_inc(e,gp,2,2))
      // elem%sigma(e,gp,2,2) = elem%sigma(e,gp,2,2) + c * ((1.0-dom%mat_nu)*elem%str_inc(e,gp,2,2)+dom%mat_nu*elem%str_inc(e,gp,1,1))
      // elem%sigma(e,gp,1,2) = elem%sigma(e,gp,1,2) + (1.0-2.0*dom%mat_nu) * elem%str_inc(e,gp,1,2)
      // elem%sigma(e,gp,2,1) = elem%sigma(e,gp,1,2)      
    // end do
  // end do 
  par_loop(e,m_elem_count){ 
    tensor3 Sigma;
    for (int gp=0;gp<m_gp_count;gp++){
      int offset_s = e * m_gp_count + gp;   //SCALAR
      int offset_t = offset_s * 6 ; //SYM TENSOR

       Sigma = FromFlatSym(m_sigma, offset_t );
    
      double c = mat[e]->Elastic().E() / ((1.0+mat[e]->Elastic().Poisson())*(1.0-2.0*mat[e]->Elastic().Poisson())) ;  
      tensor3 str_inc     = FromFlatSym(m_str_rate,     offset_t +gp)*dt;
      Sigma.xx +=  c * str_inc.xx + mat[e]->Elastic().Poisson()*str_inc.yy;
      Sigma.yy +=  c * str_inc.yy + mat[e]->Elastic().Poisson()*str_inc.xx;
      Sigma.xy +=  2.0 * mat[e]->Elastic().G() * str_inc.xy;
      //Sigma.zz += c * ((1.0 - mat[e]->Elastic().Poisson()*(str_inc.xx+str_inc.yy) *str_inc.zz);
      // ! elem%sigma(e,gp,1,1) = elem%sigma(e,gp,1,1) + c * (elem%str_inc(e,gp,1,1)+dom%mat_nu*elem%str_inc(e,gp,2,2))
      // ! elem%sigma(e,gp,2,2) = elem%sigma(e,gp,2,2) + c * (elem%str_inc(e,gp,2,2)+dom%mat_nu*elem%str_inc(e,gp,1,1))
      // ! elem%sigma(e,gp,1,2) = elem%sigma(e,gp,1,2) + 2.0* mat_G * elem%str_inc(e,gp,1,2)
      // ! elem%sigma(e,gp,2,1) = elem%sigma(e,gp,1,2)

      ToFlatSymPtr(Sigma, m_sigma,offset_t);  //TODO: CHECK IF RETURN VALUE IS SLOWER THAN PASS AS PARAM		
      
    }//gp
  }//el < elcount
    // end do !gauss point
  // end do      
      
}

// 
// !!!!!! IT ASSUMES PRESSURE AND STRAIN RATES ARE ALREADY CALCULATED
// !!!!!! (AT t+1/2 to avoid stress at rigid rotations, see Benson 1992)
dev_t void Domain_d::CalcStressStrain(double dt){

  par_loop(e,m_elem_count){
      //printf("calculating sigma \n");
        // Jaumann rate terms
      tensor3 RotationRateT,SRT,RS;
      tensor3 RotRate;
      tensor3 StrRate;
      tensor3 ShearStress;
      tensor3 Sigma;
      tensor3 Strain_pl_incr;
      tensor3 Strain;

    //printf("calculating sigma %d\n", e);
    for (int gp=0;gp<m_gp_count;gp++){
      int offset_s = e * m_gp_count + gp;   //SCALAR
      int offset_t = offset_s * 6 ; //SYM TENSOR
      ShearStress = FromFlatSym(m_tau,          offset_t );
      StrRate     = FromFlatSym(m_str_rate,     offset_t );
      RotRate     = FromFlatAntiSym(m_rot_rate, offset_t );
      Strain      = FromFlatSym(m_eps,          offset_t );

      // printf("TEST PREV SS\n");
      // print(ShearStress);
      

      SRT = ShearStress * Trans(RotRate);
      RS  = RotRate * ShearStress;
      // SRT = MatMul(elem%shear_stress(e,gp,:,:),transpose(elem%rot_rate(e,gp,:,:)))
      // RS  = MatMul(elem%rot_rate(e,gp,:,:), elem%shear_stress(e,gp,:,:))
      
      //ShearStress = dt * (2.0 * ( StrRate -SRT + RS));
      //tensor3 test = StrRate - 1.0/3.0*Trace(StrRate) * Identity();
      
      ShearStress	= ShearStress  + dt*(2.0* mat[e]->Elastic().G()*(StrRate - 1.0/3.0*Trace(StrRate) * Identity() ) + SRT+RS);
      
      //TODO: SAVE IN DOMAIN?
      double eff_strain_rate = sqrt ( 	0.5*( (StrRate.xx-StrRate.yy)*(StrRate.xx-StrRate.yy) +
                                  (StrRate.yy-StrRate.zz)*(StrRate.yy-StrRate.zz) +
                                  (StrRate.zz-StrRate.xx)*(StrRate.zz-StrRate.xx)) + 
                            3.0 * (StrRate.xy*StrRate.xy + StrRate.yz*StrRate.yz + StrRate.zx*StrRate.zx)
                          );

      double Et;
      
      tensor3 Sigma_trial = -p[offset_s] * Identity() + ShearStress;

    tensor3 s_trial = Sigma_trial - (1.0/3.0)*Trace(Sigma_trial)*Identity();


      double J2 = 0.5*(s_trial.xx*s_trial.xx +  2.0*s_trial.xy*s_trial.xy + 
                                      2.0*s_trial.xz*s_trial.xz + 
                     s_trial.yy*s_trial.yy+  
                                      2.0*s_trial.yz*s_trial.yz +               
                     s_trial.zz*s_trial.zz                 
                                     
                    );
      double sig_trial = sqrt(3.0*J2);
      
      
      if(mat[e]->Material_model == HOLLOMON ){
        sigma_y[e] = CalcHollomonYieldStress(pl_strain[e],mat[e]); 
      } else if (mat[e]->Material_model == JOHNSON_COOK){
        sigma_y[e] = CalcJohnsonCookYieldStress(pl_strain[e], eff_strain_rate, T[e], mat[e]); 
        }
        
        //printf("sy %.3f\n",sigma_y[e]);
  // Inside your plasticity block where (sigma_y[e] < sig_trial)
  double dep = 0.0; //Incremental plastic strain
  if (sigma_y[e] < sig_trial) {
	  double Ep = 0.0; //Hardening
      //printf("YIELD: %.5e\n", sigma_y[e]);
      // Calculate scaling factor
      double scale_factor = sigma_y[e] / sig_trial;
      
      // Scale the entire deviatoric stress tensor
      ShearStress = ShearStress * scale_factor;
      
      // Reconstruct the full stress tensor
      //Sigma = -p[offset_s] * Identity() + ShearStress;
      

      
      // Calculate new equivalent stress for verification
      //~ tensor3 s_new = Sigma - (1.0/3.0)*Trace(Sigma)*Identity();
      //~ double J2_new = 0.5*(s_new.xx*s_new.xx + 2.0*s_new.xy*s_new.xy + 
                          //~ 2.0*s_new.xz*s_new.xz + s_new.yy*s_new.yy +  
                          //~ 2.0*s_new.yz*s_new.yz + s_new.zz*s_new.zz);
      //~ double sig_new = sqrt(3.0*J2_new);
      

      
      // Update tangent modulus if needed
      if(mat[e]->Material_model == HOLLOMON) {
          Et = CalcHollomonTangentModulus(pl_strain[e], mat[e]);
      } else if (mat[e]->Material_model == JOHNSON_COOK){
		  Et = CalcJohnsonCookTangentModulus(pl_strain[e],eff_strain_rate, T[e], mat[e]);
      }
	  
	  Ep = mat[e]->Elastic().E()*Et/(mat[e]->Elastic().E()-Et);
	  
	  if (Ep<0) Ep = 1.*mat[e]->Elastic().E();
    dep = (sig_trial - sigma_y[e]) / (3.0 * mat[e]->Elastic().G()+Ep);//Fraser, Eq 3-49 TODO: MODIFY FOR TANGENT MODULUS = 0
      // Update plastic strain
    pl_strain[e] += dep;
  
    

  }//IF PLASTIC



      
      //printf("Shear Stress\n");
      //print(ShearStress);
      // elem%shear_stress(e,gp,:,:)	= dt * (2.0 * mat_G *(elem%str_rate(e,gp,:,:) - 1.0/3.0 * &
                                   // (elem%str_rate(e,gp,1,1)+elem%str_rate(e,gp,2,2)+elem%str_rate(e,gp,3,3))*ident) &
                                   // +SRT+RS) + elem%shear_stress(e,gp,:,:)
                                   
      // elem%sigma(e,gp,:,:) = -elem%pressure(e,gp) * ident + elem%shear_stress(e,gp,:,:)	!Fraser, eq 3.32
      Sigma = -p[offset_s] * Identity() + ShearStress;

       J2 = 0.5*(Sigma.xx*Sigma.xx +  2.0*Sigma.xy*Sigma.xy + 
                                      2.0*Sigma.xz*Sigma.xz + 
                     Sigma.yy*Sigma.yy+  
                                      2.0*Sigma.yz*Sigma.yz +               
                     Sigma.zz*Sigma.zz                 
                                     
                    );

  if (e==0){
    //printf("Epszz %.4e\n",Strain.zz);
    
    //printf("Eps R %.4e\n",StrRate.zz);
    
    //printf("Sigma zz %.4e radial %.4e mu x szz %.4e\n", Sigma.zz,sqrt(Sigma.xx*Sigma.xx+Sigma.yy*Sigma.yy),0.3*Sigma.zz);
  }      

      if (m_thermal)
        m_q_plheat[e]  = 0.0;
      
      if (dep > 0.0){
        ///////To calculate inelastic fraction and plastic work
        double f = dep/sigma_y[e];
        Strain_pl_incr.xx = f*(Sigma.xx-0.5*(Sigma.yy + Sigma.zz ));
        Strain_pl_incr.yy = f*(Sigma.yy-0.5*(Sigma.xx + Sigma.zz ));
        Strain_pl_incr.zz = f*(Sigma.zz-0.5*(Sigma.xx + Sigma.yy ));
        Strain_pl_incr.xy = Strain_pl_incr.yx = 1.5*f*(Sigma.xy);
        Strain_pl_incr.xz = Strain_pl_incr.zx = 1.5*f*(Sigma.xz);
        Strain_pl_incr.yz = Strain_pl_incr.zy = 1.5*f*(Sigma.yz);
      
      
        if (m_thermal){
          
            tensor3 depdt = 1./dt * Strain_pl_incr;
          //~ // // Double inner product, Fraser 3-106
          //printf("depdt %.3e\n", depdt);
          //~ // //cout << depdt<<endl;
          m_q_plheat[e] 	= 
            Sigma.xx * depdt.xx + 
            2.0*Sigma.xy * depdt.yx + 2.0 * Sigma.xz*depdt.zx +              //~ // ;
            Sigma.yy*depdt.yy +      //~ // //cout << "plastic heat "<<q_plheat<<endl;
            2.0*Sigma.yz*depdt.yz +     //~ // ps_energy[e] += q_plheat * dt;
            Sigma.zz*depdt.zz; //Parallel
          //~ if (m_q_plheat[e]>1.0e-5)
            //~ printf("m_q_plheat %.3e\n",m_q_plheat[e]);
        
              
            //printf("plastic heat per element %.3e\n",m_q_plheat[e]*vol[e]);
          
          }
      
      } 
      double Ep = 0;
      Strain = Strain + dt*StrRate;
      //printf("Strain %.4e\n",Strain.zz);
      
			//cout << "dep: "<<dep<<endl;
			//pl_strain[e] += dep;
			//delta_pl_strain = dep; // For heating work calculation
			
      //if (Material_model < JOHNSON_COOK ) //In johnson cook there are several fluences per T,eps,strain rate
			//if (Material_model == BILINEAR )
			//  Sigmay += dep*Ep;

      ///// OUTPUT TO Flatten arrays
      ToFlatSymPtr(Sigma, m_sigma,offset_t);  //TODO: CHECK IF RETURN VALUE IS SLOWER THAN PASS AS PARAM		
      //ToFlatSymPtr(Strain, 	strain,6*i);		
      ToFlatSymPtr(ShearStress, m_tau, offset_t);
      ToFlatSymPtr(Strain,      m_eps, offset_t);
      ToFlatSymPtr(StrRate,     m_epsr, offset_t);
      
      ToFlatSymPtr(Strain_pl_incr, m_strain_pl_incr, offset_t);
      
    }//gp
  }//el < elcount
    // end do !gauss point
  // end do
  //printf("ELEMENT %d SIGMA\n");
 
}



dev_t void Domain_d:: calcElemHourglassForces()
{
  if (m_dim == 3 && m_nodxelem == 4){
        // Tetrahedron hourglass control parameters
        const double hourglass_coeff = 0.15; // Typical value for stabilization
        const double small_volume = 1e-15;  // Threshold for degenerate elements
        
        // Hourglass base vectors for tetrahedra (4 modes)
        const double Sig[4][4] = {
            { 1.0, -1.0,  1.0, -1.0},
            { 1.0, -1.0, -1.0,  1.0},
            { 1.0,  1.0, -1.0, -1.0},
            {-1.0,  1.0, -1.0,  1.0}
        };

        par_loop(e, m_elem_count) {
            if (m_gp_count == 1) { // Single integration point
                int offset = e * m_nodxelem * m_dim;
                
                // Skip degenerate elements
                if (vol[e] < small_volume) {
                    for (int n = 0; n < m_nodxelem; n++) {
                        for (int d = 0; d < m_dim; d++) {
                            m_f_elem_hg[offset + n * m_dim + d] = 0.0;
                        }
                    }
                    continue;
                }

                // Calculate characteristic element size
                double elem_size = pow(vol[e], 1.0/3.0);
                
                // Compute stabilization coefficient
                double c_h = hourglass_coeff * rho[e] * mat[e]->cs0 * elem_size * elem_size;
                
                // Compute hourglass modes (4 modes for tetrahedra)
                double hmod[3][4] = {0}; // [dim][mode]
                
                for (int d = 0; d < m_dim; d++) {
                    for (int j = 0; j < 4; j++) {
                        for (int n = 0; n < 4; n++) {
                            hmod[d][j] += getVElem(e, n, d) * Sig[j][n];
                        }
                    }
                }
                
                // Apply hourglass forces
                for (int n = 0; n < 4; n++) {
                    for (int d = 0; d < m_dim; d++) {
                        double force = 0.0;
                        for (int j = 0; j < 4; j++) {
                            force -= hmod[d][j] * Sig[j][n];
                        }
                        m_f_elem_hg[offset + n * m_dim + d] = c_h * force;
                    }
                }
            }
        }
    }
    
  else {
    
  int jmax;
  if (m_dim==2) jmax = 1;
  else          jmax = 4;
  // real(fp_kind), dimension(4, nodxelem) :: Sig !! 4 COLUMNVECTORS IN 2D ONLY first is used
  // real(fp_kind), dimension(nodxelem,dim):: vel!!!!DIFFERENT FROM vele which is an 8 x 1 vector
  // real(fp_kind), dimension(dim,4) :: hmod
  double hmod[3][4] ={{0.0,0.0,0.0,0.0},
  {0.0,0.0,0.0,0.0},
  {0.0,0.0,0.0,0.0}}; //m_dim,4
  //Matrix hmod(m_dim,4);
    Matrix Sig(4,m_nodxelem);
 
  
      
  // //Sig nodxelem
  if (m_dim==3){
    
    double sig_[4][8] = { { 0.125, 0.125,-0.125,-0.125,-0.125,-0.125, 0.125, 0.125},
                          { 0.125,-0.125,-0.125, 0.125,-0.125, 0.125, 0.125,-0.125},
                          { 0.125,-0.125, 0.125,-0.125, 0.125,-0.125, 0.125,-0.125},
                          {-0.125, 0.125,-0.125, 0.125, 0.125,-0.125, 0.125,-0.125}};  
    
    for (int i=0;i<4;i++)
      for (int n=0;n<m_nodxelem;n++)
        Sig.Set(i,n,sig_[i][n]*8.0);

  } else if (m_dim == 2){
    double sig_[1][4] ={{0.25, -0.25, 0.25,-0.25}};    
    
    for (int n=0;n<m_nodxelem;n++)
      Sig.Set(0,n,sig_[0][n]*4.0);
      //Sig.Set(0,n,0.0);
  }

  //double f = 1/8;
  par_loop(e,m_elem_count){   
        
  if (m_gp_count==1){
      int offset = e * m_nodxelem*m_dim;   //SCALAR  
      //double hmod[m_dim][4];

    for (int d=0;d<m_dim;d++)
        for (int n=0;n<jmax;n++)
          hmod[d][n] = 0.0;
      
    for (int d=0;d<m_dim;d++)
        for (int n=0;n<m_nodxelem;n++)
          m_f_elem_hg[offset + n*m_dim + d] = 0.0;



      //for (int n=0;n<m_nodxelem;n++)
      //  printf("Elem %d Vel %.6e %.6e %.6e\n",e, getVElem(e,n,0),getVElem(e,n,1),getVElem(e,n,2)); ////DIM  

      for (int d=0;d<m_dim;d++)
        for (int j=0;j<jmax;j++)
          for (int n=0;n<m_nodxelem;n++){
            hmod[d][j] += getVElem(e,n,d) * Sig.getVal(j,n); ////DIM
          }

      // !!!!!!!!! GOUDREAU 1982
      for (int d=0;d<m_dim;d++)
        for (int j=0;j<jmax;j++)
          for (int n=0;n<m_nodxelem;n++)
            //elem%hourg_nodf(e,n,:) = elem%hourg_nodf(e,n,:) - hmod(:,j)*Sig(j,n)
            m_f_elem_hg[offset + n*m_dim + d] -= hmod[d][j] * Sig.getVal(j,n);
          // end do
      // end do
      // c_h  = 0.06 * elem%vol(e)**(0.6666666) * elem%rho(e,1) * 0.25 * mat_cs0
      double c_h = 0.15 * pow(vol[e], 0.6666666) * rho[e] * 0.2500 * mat[e]->cs0;
      //printf("c_h %.6e\n", c_h);


      for (int n=0;n<m_nodxelem;n++){      
        for (int d=0;d<m_dim;d++){
          m_f_elem_hg[offset + n*m_dim + d] *= c_h;
        }
        //printf("hg forces el %d node %d: %f %f %f\n",e, n,m_f_elem_hg[offset + n*m_dim],m_f_elem_hg[offset + n*m_dim + 1],m_f_elem_hg[offset + n*m_dim + 2]  );
      }

  } //gp ==1
  }//ELEM
  
  } //non tetra
}


  

// // Alternative simplified version focusing only on volumetric stabilization
// dev_t void Domain_d::calcElemHourglassForces_Simple()
// {
  // if (m_dim != 3 || m_nodxelem != 4) {
    // return; // Only for 3D tetrahedral elements
  // }
  
  // par_loop(e, m_elem_count) {
    
    // int offset = e * m_nodxelem * m_dim;
    
    // // Initialize forces
    // for (int n = 0; n < m_nodxelem; n++) {
      // for (int d = 0; d < m_dim; d++) {
        // m_f_elem_hg[offset + n * m_dim + d] = 0.0;
      // }
    // }
    
    // // Simple pressure-based stabilization
    // // This is similar to what you have in calcElemPressureANP
    // double vol_ratio = vol[e] / vol_0[e];
    // double pressure_dev = mat[e]->Elastic().BulkMod() * (1.0 - vol_ratio);
    
    // // Apply pressure-based hourglass control
    // double hg_coeff = 0.05 * fabs(pressure_dev) * vol[e];
    
    // // Get shape function derivatives (simplified for tetrahedra)
    // double dNdxi[4] = {-1.0, 1.0, 0.0, 0.0};
    // double dNdeta[4] = {-1.0, 0.0, 1.0, 0.0};  
    // double dNdzeta[4] = {-1.0, 0.0, 0.0, 1.0};
    
    // for (int n = 0; n < 4; n++) {
      // // Simple hourglass forces based on shape function derivatives
      // m_f_elem_hg[offset + n * m_dim + 0] = hg_coeff * dNdxi[n] * pressure_dev;
      // m_f_elem_hg[offset + n * m_dim + 1] = hg_coeff * dNdeta[n] * pressure_dev;
      // m_f_elem_hg[offset + n * m_dim + 2] = hg_coeff * dNdzeta[n] * pressure_dev;
    // }
  // }
// }

////// CLASSIC WILKINS
dev_t void Domain_d::calcArtificialViscosity() {
  //double alpha = 2.0;  // Increased from 1.0 (more linear damping)
  //double beta = 0.2;   // Increased from 0.05 (more quadratic damping)
  
  double alpha = m_artifvisc[0];
  double beta  = m_artifvisc[1];
    
  double q_max = 1e9;

  //par_loop(e, m_elem_count) {
  for (int e=0;e<m_elem_count;e++){
    int offset_t = e * 6;
    double c = sqrt(mat[e]->Elastic().BulkMod() / rho[e]);
    tensor3 StrRate = FromFlatSym(m_str_rate, e * m_gp_count * 6);
    double eps_v = Trace(StrRate);

    // Apply viscosity for BOTH compression and tension
    if (abs(eps_v) > 1e-12) {  // Small threshold to avoid noise
      double l = pow(vol[e], 1.0/3.0);
      double q = alpha * c * abs(eps_v) * l + beta * pow(eps_v * l, 2);
      q = std::min(q, q_max);

      // Subtract q from diagonal stresses (sign depends on expansion/compression)
      double q_signed = (eps_v > 0) ? -q : q;
      m_sigma[offset_t + 0] += q_signed;
      m_sigma[offset_t + 1] += q_signed;
      m_sigma[offset_t + 2] += q_signed;
    }
  }
}

  /////DUMMY IN CASE OF CPU
  __global__ void calcElemPressureKernel(Domain_d *dom_d){		
    dom_d->calcElemPressure();
  }

  __global__ void calcElemDensityKernel  (Domain_d *dom_d){
    dom_d->calcElemDensity();  
  }

  __global__ void calcElemVolKernel(Domain_d *dom_d){
		
		dom_d->CalcElemVol();
  }
  __global__ void calcElemInitialVolKernel(Domain_d *dom_d){
		
		dom_d->CalcElemInitialVol();
  }
  
  __global__ void calcStressStrainKernel(Domain_d *dom_d, double dt){
    dom_d->CalcStressStrain(dt);
  }

  __global__ void calcElemForcesKernel(Domain_d *dom_d){
		
		dom_d->calcElemForces();
  }

  __global__ void calcElemHourglassForcesKernel(Domain_d *dom_d){
		
		dom_d->calcElemHourglassForces();
  }
  
  __global__ void calcAccelKernel(Domain_d *dom_d){
		
		dom_d->calcAccel();
  }

__global__ void CalcNodalVolKernel        (Domain_d *dom_d){
  dom_d->CalcNodalVol();
  }
__global__ void CalcNodalMassFromVolKernel(Domain_d *dom_d){
    dom_d->CalcNodalMassFromVol();
  }

 
}; //Namespace MetFEM
