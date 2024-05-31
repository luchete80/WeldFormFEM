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
    Matrix *str_rate = new Matrix(m_dim,m_dim); //TODO: MAKE SYMM MATRIX
    Matrix *rot_rate = new Matrix(m_dim,m_dim); //TODO: MAKE SYMM MATRIX
  //Matrix *dHxy_detJ_loc = new Matrix(m_dim, m_nodxelem);

//  if (e < m_elem_count) {

    for (int gp=0;gp<m_gp_count;gp++){
      int offset = e * m_gp_count * 6 + gp;
	  int offset_det = e * m_gp_count;
		str_rate->SetZero();
    rot_rate->SetZero();
	//printf ("offset, %f , det %f\n", offset, m_detJ[offset + gp]);
      double f = 1.0 / m_detJ[offset_det + gp];
      //printf("f factor: %f\n",f);
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
          //printf("d %d n %d deriv %f vele %f\n",d, n, getDerivative(e,gp,d,n),  getVElem(e,n,d));
          
          str_rate->Set(d,d, str_rate->getVal(d,d) + getDerivative(e,gp,d,n) * f * getVElem(e,n,d));
          rot_rate->Set(d,d, 0.0);

          // elem%str_rate(e,gp, d,d) = elem%str_rate(e,gp, d,d) + temp(d,n) * elem%vele (e,dim*(n-1)+d,1) 
          // elem%rot_rate(e,gp, d,d) = 0.0d0
          
        }//dim
        // !!!! TO AVOID ALL MATMULT
        str_rate->Set(0,1, str_rate->getVal(0,1) + f *(getDerivative(e,gp,1,n) * getVElem(e,n,0) +
                                                       getDerivative(e,gp,0,n) * getVElem(e,n,1)));
        rot_rate->Set(0,1, rot_rate->getVal(0,1) + f* (getDerivative(e,gp,1,n) * getVElem(e,n,0) - 
                                                       getDerivative(e,gp,0,n) * getVElem(e,n,1)));
        if (m_dim == 3) {
          str_rate->Set(1,2, str_rate->getVal(1,2) + f *(getDerivative(e,gp,2,n) * getVElem(e,n,1) +
                                                         getDerivative(e,gp,1,n) * getVElem(e,n,2)));
          str_rate->Set(0,2, str_rate->getVal(0,2) + f *(getDerivative(e,gp,2,n) * getVElem(e,n,0) +
                                                         getDerivative(e,gp,0,n) * getVElem(e,n,2)));

          rot_rate->Set(1,2, rot_rate->getVal(1,2) + f *(getDerivative(e,gp,2,n) * getVElem(e,n,1) -
                                                         getDerivative(e,gp,1,n) * getVElem(e,n,2)));
          rot_rate->Set(0,2, rot_rate->getVal(0,2) + f *(getDerivative(e,gp,2,n) * getVElem(e,n,0) -
                                                         getDerivative(e,gp,0,n) * getVElem(e,n,2)));
        

        }// end if     
      }// end do !Nod x elem
      //printf ("test %fn", test);
      str_rate->Set(0,1, str_rate->getVal(0,1) *0.5);
      str_rate->Set(0,2, str_rate->getVal(0,2) *0.5);  
      str_rate->Set(1,2, str_rate->getVal(1,2) *0.5);
      
      rot_rate->Set(0,1, rot_rate->getVal(0,1) *0.5);
      rot_rate->Set(0,2, rot_rate->getVal(0,2) *0.5);  
      rot_rate->Set(1,2, rot_rate->getVal(1,2) *0.5);
      
      str_rate->Set(1,0,  str_rate->getVal(0,1));   
      str_rate->Set(2,1,  str_rate->getVal(1,2));  
      str_rate->Set(2,0,  str_rate->getVal(0,2));

      rot_rate->Set(1,0, -rot_rate->getVal(0,1));      
      rot_rate->Set(2,1, -rot_rate->getVal(1,2));  
      rot_rate->Set(2,0, -rot_rate->getVal(0,2)); 

      str_rate->ToFlatSymPtr(m_str_rate, offset);
      rot_rate->ToFlatSymPtr(m_rot_rate, offset); //UPPER PART

      printf("Strain Rate\n");
      str_rate->Print();

      // printf("Rot Rate\n");
      // str_rate->Print();
      
      // !elem%str_rate(e,gp,:,:) = matmul(elem%bl(e,gp,:,:),elem%vele (e,:,:)) 
      // !print *, "simlpified strain rate "

      //Inverse test
      // Matrix *test = new Matrix(3,3);
      // Matrix *invtest = new Matrix(3,3);
      // //printf("A\n");
      // test->Set(0,0,1);test->Set(0,1,1);test->Set(0,2,1);
      // test->Set(1,0,1);test->Set(1,1,2);test->Set(1,2,2);      
      // test->Set(2,0,1);test->Set(2,1,2);test->Set(2,2,3);
      // InvMat(*test,invtest);
      // ////printf("inv A\n");
      // test->Print();
      // invtest->Print();
        // delete test, invtest;
      } // Gauss Point
      delete str_rate,rot_rate;
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
    if (m_gp_count == 1)  w = pow(2.0, m_dim);      
    else                  w = 1.0;
    int offset = m_gp_count * e;
    vol[e] = 0.0;
    for (int gp=0;gp<m_gp_count;gp++){
      vol[e] += m_detJ[offset] * w;
    }    
  }
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
      //printf ("Elem %d vol %f\n", e, vol_0[e]);
    }      
  }  
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
      for (int i=0;i<3;i++)
        for (int j=0;j<3;j++){
          //printf("SIGMA %d %d %.6e\n",i,j,getSigma(e,gp,i,j));
        }
      for (int n=0; n<m_nodxelem;n++) {
        for (int d=0;d<m_dim;d++){
          m_f_elem[offset + n*m_dim + d] += getDerivative(e,gp,d,n) * getSigma(e,gp,d,d);
        }
      
        if (m_dim == 2){
          m_f_elem[offset + n*m_dim    ] +=  getDerivative(e,gp,1,n) * getSigma(e,gp,0,1);
          m_f_elem[offset + n*m_dim + 1] +=  getDerivative(e,gp,0,n) * getSigma(e,gp,0,1);
        } else {
          //printf("offset %d\n", offset + n*m_dim    );
          //printf ("sigma 0 1 %f\n", getSigma(e,gp,0,1));
          m_f_elem[offset + n*m_dim    ] +=  getDerivative(e,gp,1,n) * getSigma(e,gp,0,1) +
                                             getDerivative(e,gp,2,n) * getSigma(e,gp,0,2);
          m_f_elem[offset + n*m_dim + 1] +=  getDerivative(e,gp,0,n) * getSigma(e,gp,0,1) + 
                                             getDerivative(e,gp,2,n) * getSigma(e,gp,1,2);        
          m_f_elem[offset + n*m_dim + 2] +=  getDerivative(e,gp,1,n) * getSigma(e,gp,1,2) + 
                                             getDerivative(e,gp,0,n) * getSigma(e,gp,0,2);     
        }

      
        // if (m_dim == 2){
          // m_f_elem[offset + n*m_dim    ] +=  getDerivative(e,gp,1,n) * sigma.xy;
          // m_f_elem[offset + n*m_dim + 1] +=  getDerivative(e,gp,0,n) * sigma.xy;
        // } else {
          // //printf("offset %d\n", offset + n*m_dim    );
          // //printf ("sigma 0 1 %f\n", getSigma(e,gp,0,1));
          // m_f_elem[offset + n*m_dim    ] +=  getDerivative(e,gp,1,n) * sigma.xy +
                                             // getDerivative(e,gp,2,n) * sigma.xz;
          // m_f_elem[offset + n*m_dim + 1] +=  getDerivative(e,gp,0,n) * sigma.xy + 
                                             // getDerivative(e,gp,2,n) * sigma.yz;        
          // m_f_elem[offset + n*m_dim + 2] +=  getDerivative(e,gp,1,n) * sigma.yz + 
                                             // getDerivative(e,gp,0,n) * sigma.xz;     
        // }

        
      }// nod x elem
      
      /// TO REPLACE EASIER!!!!!!
        // for (int i = 0; i < m_nodxelem; i++) {
            // for (int j = 0; j < m_dim; j++) {
              // m_f_elem[m_dim*i+j] = 0.0;
              // for (int k = 0; k < m_dim; k++) 
                // //forces[i][j] += dNdX[gp][k][i] * s2[k][j]*m_detJ[gp] * gauss_weights[gp];
                // m_f_elem[m_dim*i+j] += getDerivative(e,gp,k,i) * getSigma(e,gp,k,j)*m_detJ[gp] ;
              
              // //printf ("forces %e",forces[i][j]);
            // }
        // }

    } // Gauss Point
    for (int n=0; n<m_nodxelem;n++) {
      for (int d=0;d<m_dim;d++){
        m_f_elem[offset + n*m_dim + d] *= w;
      }
       //printf ("elem %d forces %f %f %f\n",e,m_f_elem[offset + n*m_dim + 0],
       //                                    m_f_elem[offset + n*m_dim + 1],
       //                                    m_f_elem[offset + n*m_dim + 2]);
    }  
  }//if e<elem_count
}

dev_t void Domain_d::calcElemPressure(){

  par_loop(e,m_elem_count){
    printf("calc pressure \n");
    int offset_t = e * m_gp_count *6;
    Matrix *sigma   = new Matrix(3,3);
    //Matrix str_inc(m_dim,m_dim);
    //sigma.FromFlatSymPtr();
    double trace;
    double press_inc = 0.0;
    for (int gp=0;gp<m_gp_count;gp++){
      trace = 0.0;
      tensor3 str_inc     = FromFlatSym(m_str_rate,     offset_t +gp)*dt;
      //printf("str inc, dt %f\n", dt);print(str_inc);
      press_inc += Trace(str_inc);
    }//gauss point
    press_inc = -press_inc/m_gp_count;
    //  //printf("trace %f\n",trace);
    int offset = e * m_gp_count ;
    for (int gp=0;gp<m_gp_count;gp++){
      //  //printf("bulk mod:%f, press inc%f\n", mat[e]->Elastic().BulkMod(),press_inc);
      trace = 0.0;
      for (int d = 0; d<3;d++) trace += getSigma(e,gp,d,d);
      
      p[offset + gp] = -1.0/3.0 * trace + mat[e]->Elastic().BulkMod() * press_inc;
      printf("pressure %f\n",p[offset + gp]);
    }
    delete sigma;
  } // e< elem_count
}
    
// subroutine calc_elem_pressure_from_strain (modK)
  // implicit none
  // real(fp_kind), intent(in) :: modK
  // integer :: e, gp  
  // real(fp_kind) :: press_inc
  
  // gp = 1
  // do e = 1, elem_count
    // press_inc = 0.0d0
    // do gp = 1, elem%gausspc(e)
      // press_inc = press_inc + trace (elem%str_inc(e,gp,:,:))
    // end do
    // press_inc = -press_inc/elem%gausspc(e)
    // ! print *, "press inc ", press_inc
    // do gp = 1, elem%gausspc(e)    
          // elem%pressure(e,gp) = -1.0/3.0 * trace (elem%sigma(e,gp,:,:)) + press_inc * modK
    // end do
    // ! print *, "mod K", modK
    // ! print *, "strain inc", elem%str_inc(e,gp,:,:)
    // ! print *, "press_inc ", press_inc
    // ! print *, "elem%pressure(e,gp) FROM STRAIN", elem%pressure(e,1)
  // end do
// end subroutine 

///////////////////// ONLY FOR TESTING, ELASTIC STRAIN

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

  tensor3 Sigma;
  par_loop(e,m_elem_count){ 
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


	
  // implicit none
  // real(fp_kind) :: SRT(3,3), RS(3,3), ident(3,3)
  // integer :: e,gp
  // real(fp_kind) ,intent(in):: dt
  
  // real(fp_kind) :: p00
  
  // p00 = 0.
  
  // ident = 0.0d0
  // ident (1,1) = 1.0d0; ident (2,2) = 1.0d0;  ident (3,3) = 1.0d0
  printf("calculating sigma \n");
  par_loop(e,m_elem_count){

        // Jaumann rate terms
      tensor3 RotationRateT,SRT,RS;
      tensor3 RotRate;
      tensor3 StrRate;
      tensor3 ShearStress;
      tensor3 Sigma;

    //printf("calculating sigma %d\n", e);
    for (int gp=0;gp<m_gp_count;gp++){
      int offset_s = e * m_gp_count + gp;   //SCALAR
      int offset_t = offset_s * 6 ; //SYM TENSOR
      ShearStress = FromFlatSym(m_tau,          offset_t );
      StrRate     = FromFlatSym(m_str_rate,     offset_t );
      RotRate     = FromFlatAntiSym(m_rot_rate, offset_t );

      // printf("TEST PREV SS\n");
      // print(ShearStress);
      

      SRT = ShearStress * Trans(RotRate);
      RS  = RotRate * ShearStress;
      // SRT = MatMul(elem%shear_stress(e,gp,:,:),transpose(elem%rot_rate(e,gp,:,:)))
      // RS  = MatMul(elem%rot_rate(e,gp,:,:), elem%shear_stress(e,gp,:,:))
      
      //ShearStress = dt * (2.0 * ( StrRate -SRT + RS));
      tensor3 test = StrRate-1.0/3.0*(StrRate.xx+StrRate.yy+StrRate.zz)*Identity();

      ShearStress	= ShearStress  + dt*(2.0* mat[e]->Elastic().G()*(StrRate - 1.0/3.0*Trace(StrRate) * Identity() ) + SRT+RS);
      //printf("Shear Stress\n");
      //print(ShearStress);
      // elem%shear_stress(e,gp,:,:)	= dt * (2.0 * mat_G *(elem%str_rate(e,gp,:,:) - 1.0/3.0 * &
                                   // (elem%str_rate(e,gp,1,1)+elem%str_rate(e,gp,2,2)+elem%str_rate(e,gp,3,3))*ident) &
                                   // +SRT+RS) + elem%shear_stress(e,gp,:,:)
                                   
      // elem%sigma(e,gp,:,:) = -elem%pressure(e,gp) * ident + elem%shear_stress(e,gp,:,:)	!Fraser, eq 3.32
      Sigma = -p[offset_s] * Identity() + ShearStress;
      // printf("SHEAR STRESS\n");
      // print(ShearStress);

      // printf("STR RATE\n");
      // print(StrRate);
      
      // printf("ELEMENT %d SIGMA\n");
      // print(Sigma);

      ///// OUTPUT TO Flatten arrays
      ToFlatSymPtr(Sigma, m_sigma,offset_t);  //TODO: CHECK IF RETURN VALUE IS SLOWER THAN PASS AS PARAM		
      //ToFlatSymPtr(Strain, 	strain,6*i);		
      ToFlatSymPtr(ShearStress, m_tau, offset_t);
      
    }//gp
  }//el < elcount
    // end do !gauss point
  // end do
  printf("ELEMENT %d SIGMA\n");
 
}


// dev_t void Domain_d::CalcStressStrain(const double dt){

    // double srt[m_gp_count][3][3];
    // double rs[m_gp_count][3][3];
    // double d[3][3];
    // double tau[3][3];
    // for (int gp = 0; gp < m_gp_count; gp++) {
        // for (int i = 0; i < 3; i++) {
            // for (int j = 0; j < 3; j++) {
                // srt[gp][i][j] = rs[gp][i][j] = 0.0;
                // for (int k=0;k<m_dim;k++){
                  // srt[gp][i][j] += tau[gp][i][k] * rot_rate[gp][j][k];
                  // rs[gp][i][j] += rot_rate[gp][i][k] * tau[gp][k][j];
                // }
                // dev(str_rate[gp],d);
                // tau[gp][i][j] += dt * ((2.0 * mat[0]->Elastic().G() *d[i][j]) + rs[gp][i][j] + srt[gp][i][j]);
                // stress[i][j] = tau[i][j] - p[gp] * (i == j);
                // //Sigma_tst->Set
                // printf ("stress %e",stress[gp][i][j]);
            // }
            // printf("\n");
        // }
    // }
    
    // m_sigma[0]=stress[0][0];     m_sigma[1]=stress[1][1];    m_sigma[2]=stress[2][2];  
    // m_sigma[3]=stress[0][1];     m_sigma[4]=stress[1][2];    m_sigma[5]=stress[0][2];  
    
    
// }

// //////////////////////////////////////////////////////////////////////
// __device__ void Domain_d::StressStrain(int i) {
	
	// //int i = threadIdx.x + blockDim.x*blockIdx.x;
	// double dep = 0.;
  // double Ep;  //STORAGE THIS OR THE TANGENT?
	
	// if ( i < solid_part_count ) {	
		// //Pressure = EOS(PresEq, Cs, P0,Density, RefDensity); //CALL BEFORE!

		// // Jaumann rate terms
		// tensor3 RotationRateT,SRT,RS;
		// tensor3 RotationRate;
		// tensor3 StrainRate;
		// tensor3 ShearStress;
		// tensor3 Sigma;
		// tensor3 Strain,Straina,Strainb;
		
		// double temprr[6],tempss[6],tempsr[6];
		// double tempssa[6],tempssb[6];
		// for (int k=0;k<6;k++){ //First the diagonal
			// temprr[k] = rotrate[6*i+k];
			// tempss[k] = shearstress[6*i+k];
			// tempsr[k] = strrate[6*i+k];
		// }
		// ShearStress   = FromFlatSym (tempss);	
		// StrainRate    = FromFlatSym(tempsr);
		// RotationRate  = FromFlatAntiSym(temprr);
   // // //printf(" Strain rate xx %f \n",StrainRate.zz);
		// RotationRateT = Trans(RotationRate);
		
		// SRT = ShearStress * RotationRateT;
		// RS = RotationRate * ShearStress;

		// ShearStress	= deltat*(2.0*G[i]*(StrainRate-1.0/3.0*(StrainRate.xx+StrainRate.yy+StrainRate.zz)*Identity())+SRT+RS) + ShearStress;	

    // eff_strain_rate[i] = sqrt ( 	0.5*( (StrainRate.xx-StrainRate.yy)*(StrainRate.xx-StrainRate.yy) +
                                // (StrainRate.yy-StrainRate.zz)*(StrainRate.yy-StrainRate.zz) +
                                // (StrainRate.zz-StrainRate.xx)*(StrainRate.zz-StrainRate.xx)) + 
                          // 3.0 * (StrainRate.xy*StrainRate.xy + StrainRate.yz*StrainRate.yz + StrainRate.zx*StrainRate.zx)
                        // );
                        
		// double J2	= 0.5*(ShearStress.xx*ShearStress.xx + 2.0*ShearStress.xy*ShearStress.yx +
					// 2.0*ShearStress.xz*ShearStress.zx + ShearStress.yy*ShearStress.yy +
					// 2.0*ShearStress.yz*ShearStress.zy + ShearStress.zz*ShearStress.zz);

    // //Scale back, Fraser Eqn 3-53
		// double sig_trial = sqrt(3.0*J2); 
    // if ( sigma_y[i] < sig_trial ) ShearStress = sigma_y[i]/sig_trial * ShearStress; //Yielding      
    // //std::min((Sigmay/sqrt(3.0*J2)),1.0)*ShearStressa;

   // if       (mat[i]->Material_model == HOLLOMON )       {
     
      // ////printf("Hollomon!");
      // ////printf("pl strain %f\n",pl_strain[i]);
      // //sigma_y[i] = mat [i]->CalcYieldStress(pl_strain[i]);
      // //ShowProps(mat[i]);
      // sigma_y[i] = CalcHollomonYieldStress(pl_strain[i],mat [i]);
      // ////printf ("sigma_y %f\n",sigma_y[i]);
      // //sigma_y[i] = mat [i]->testret();
      // //(*materials_ptr)->testret();
      // //sigma_y[i] = mat [i]->CalcYieldStress(pl_strain[i]); 
    // } 
		// else if  (mat[i]->Material_model == JOHNSON_COOK )   sigma_y[i] = CalcJohnsonCookYieldStress(pl_strain[i],eff_strain_rate[i],T[i], mat[i]);
		
		// sigma_eq[i] = sig_trial;	
		
		// if ( sig_trial > sigma_y[i]) {
      // if              (mat[i]->Material_model == HOLLOMON ) {
				// Et[i] = CalcHollomonTangentModulus(pl_strain[i],mat[i]); //Fraser 3.54
				// // Et_m = Et;        
      // } else if       (mat[i]->Material_model == JOHNSON_COOK ) {
        // Et[i] = CalcJohnsonCookTangentModulus(pl_strain[i], eff_strain_rate[i], T[i],mat[i]); //Fraser 3.54
      // } else if       (mat[i]->Material_model == BILINEAR ) {
        // //Ep = mat->Elastic().E()*Et/(mat->Elastic().E()-Et);
      // }
			// if (mat[i]->Material_model > BILINEAR ) {//Else Ep = 0
        // //cout << "Calculating Ep"<<endl;
				// Ep = mat[i]->Elastic().E()*Et[i]/(mat[i]->Elastic().E()-Et[i]);
				// // if (Ep < 0)
					// // cout << "ATTENTION Material Ep <0 "<<Ep<<", Et" << Et <<", platrain"<<pl_strain<<"effstrrate"<<eff_strain_rate<<endl;
			// }
      // if (Ep<0) Ep = 1.*mat[i]->Elastic().E();
      
			// dep=( sig_trial - sigma_y[i])/ (3.*G[i] + Ep);	//Fraser, Eq 3-49 TODO: MODIFY FOR TANGENT MODULUS = 0
			// pl_strain[i] += dep;	
      // ////printf("Particle %d, dep %.1e, sigtrial %.1e\n",i,dep,sig_trial);
			// sigma_eq[i] = sigma_y[i];
		// }

    // // if (mat[i]->Material_model == BILINEAR )
      // // sigma_y[i] += dep*mat[i]->Ep;
    
		// Sigma = -p[i] * Identity() + ShearStress;	//Fraser, eq 3.32

		// Strain	= deltat * StrainRate + Strain;
    
    // // if (mat[i]->Material_model==JOHNSON_COOK){
      // // //printf("JOHNSON_COOK!\n"); //test
    // // } elsr     
    // // if (mat[i]->Material_model==HOLLOMON){
      // // //printf("HOLLOMON!\n"); //test
    // // }    
    // // if (mat[i]->Material_model==BILINEAR){
      // // //printf("BILINEAR!\n"); //test
    // // }

		// ///// OUTPUT TO Flatten arrays
		// ToFlatSymPtr(Sigma, sigma,6*i);  //TODO: CHECK IF RETURN VALUE IS SLOWER THAN PASS AS PARAM		
		// ToFlatSymPtr(Strain, 	strain,6*i);		
		// ToFlatSymPtr(ShearStress, shearstress,6*i);

	// }//particle count
// }






// !!!! AFTER CALCULATING VELE 
// !!!!!THIS HOURGLASS CONTROL IS BASED ON 
// !!!!! GOUDREAU 1982 --> Used this simplified hourglass correction
// !!!!! FLANAGAN 1981
// !!!!! ATTENTION: IN FLANAGAN INTRINSIC COORDINATES ARE FROM -1/2 to 1/2
// !!!!! SO: h1=1/4(1-2r)(1-2s) (Flanagan Eqn 55). 
// !!!!! With our instrinsic from -1 to 1 , in  Table 2 Sigma should be 
// !!!! Sigma is quadratic (2D) or qubic(3D) coefficient of 
// subroutine calc_hourglass_forces
  // implicit none
  // integer :: e, n, j, d, gp, jmax
  // real(fp_kind) :: c_h
  // real(fp_kind), dimension(4, nodxelem) :: Sig !! 4 COLUMNVECTORS IN 2D ONLY first is used
  // real(fp_kind), dimension(nodxelem,dim):: vel!!!!DIFFERENT FROM vele which is an 8 x 1 vector
  // real(fp_kind), dimension(dim,4) :: hmod
  // real(fp_kind) :: test
// !real(fp_kind), dimension(1,4) :: test
  // if (dim .eq. 2) then
    // jmax = 1
  // else
    // jmax = 4
  // end if
  // !!! TODO: TO BE DEFINED IN DOMAIN ONCE hmod(:,:) = 0.0d0
  // elem%hourg_nodf(:,:,:) = 0.0d0
  // if (dim .eq. 3) then
    // !Also in Flanagan Appendix (data GB in program)
    // Sig(1,:) = [ 0.125, 0.125,-0.125,-0.125,-0.125,-0.125, 0.125, 0.125]
    // Sig(2,:) = [ 0.125,-0.125,-0.125, 0.125,-0.125, 0.125, 0.125,-0.125]
    // Sig(3,:) = [ 0.125,-0.125, 0.125,-0.125, 0.125,-0.125, 0.125,-0.125]
    // Sig(4,:) = [-0.125, 0.125,-0.125, 0.125, 0.125,-0.125, 0.125,-0.125]
    // Sig(:,:) = Sig(:,:) * 8
  // else 
    // Sig(1,:) = [ 0.25, -0.25, 0.25,-0.25] !!!  
    // Sig(:,:) = Sig(:,:) * 4
  // end if
  
  // gp = 1
  // do e=1, elem_count    
    // if (elem%gausspc(e) .eq. 1) then
      // hmod(:,:) = 0.0d0
      // !test = matmul (elem%dHxy(e,gp,:,:),transpose(Sig(:,:))) !!!!SHOULD BE ORTHOGONAL
      // !print *, "test ", test
      // !print *, "dHxy ", elem%dHxy(e,gp,:,:)
          
      // do n=1,nodxelem
        // do d=1,dim
          // vel (n,d) = nod%v(elem%elnod(e,n),d)    
        // end do
      // end do
      // do j =1,jmax !MODE
        // do n=1,nodxelem !1:4 or 8
          // !print *, "elem v ", vel (n,:)
          // hmod(:,j) = hmod(:,j) + vel (n,:)*Sig(j,n) !!!!! ":" is on dimension d, GOUDREAU EQN (30)
        // end do
      // end do
      
      
      
      // !!!!!!!!! GOUDREAU 1982
      // do n=1,nodxelem
        // do j =1,jmax 
            // elem%hourg_nodf(e,n,:) = elem%hourg_nodf(e,n,:) - hmod(:,j)*Sig(j,n)
          // end do
      // end do
      // c_h  = 0.06 * elem%vol(e)**(0.6666666) * elem%rho(e,1) * 0.25 * mat_cs0
      
      // !print *, "hourglass c ", c_h
      // elem%hourg_nodf(e,:,:) = elem%hourg_nodf(e,:,:) * c_h
      
      
      // !!!!!!!!! FLANAGAN 1981
      
      
      // ! do n=1,nodxelem
      // ! print *, "hourglass forces", elem%hourg_nodf(e,n,:) 
      // ! end do
      // ! ! elem%hourg_nodf(e,:,:) = - matmul(matmul(transpose(Sig(:,:)),Sig(:,:)),vel (:,:)) * c_h 
      
      // ! print *, "alt hourglass forces", elem%hourg_nodf(e,:,:) 
      // ! do n=1,nodxelem
        // ! test = 0
        // ! do d=1,3
          // ! test = test + elem%hourg_nodf(e,n,d) *  vel (n,d)
        // ! end do
        // ! print *, "dot vel and hg forces ", test
      // ! end do

      // ! !dimension and modes
      // ! do d=1,3
        // ! do j =1,jmax !MODE
          // ! test = 0
          // ! do n=1,nodxelem !1:4 or 8
            // ! !print *, "elem v ", vel (n,:)
            // ! !hmod(:,j) = hmod(:,j) + vel (n,:)*Sig(j,n) !!!!! ":" is on dimension d, GOUDREAU EQN (30)
            // ! test = test + vel (n,d)*Sig(j,n) !!!!! ":" is on dimension d, GOUDREAU EQN (30)
          // ! end do
          // ! ! print *, "Mode ", j, "Dim", d
          // ! print *, "mode", j, "dim ", d, "sum v x Sig", test, ", hmod ", hmod(d,j)
          // ! ! print *, "force "
        // ! end do
      // ! end do

      // ! do n=1,nodxelem
          // ! print *, "v ", vel (n,:)
      // ! end do
  // else
    // !print *, "NO HOURGLASS"
    // end if
  // end do !elemen
// end subroutine

dev_t void Domain_d:: calcElemHourglassForces()
{
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
      
    //VA_LIST NOT WORKING PROPERLY
    // SetMatVals(&Sig,32, 0.125, 0.125,-0.125,-0.125,-0.125,-0.125, 0.125, 0.125,
    // 0.125,-0.125,-0.125, 0.125,-0.125, 0.125, 0.125,-0.125,
    // 0.125,-0.125, 0.125,-0.125, 0.125,-0.125, 0.125,-0.125,
    // -0.125, 0.125,-0.125, 0.125, 0.125,-0.125, 0.125,-0.125);
                       
    ////printf("Sigma mat HG\n");
    //Sig.Print();
    // double Sig[4][8] = { { 0.125, 0.125,-0.125,-0.125,-0.125,-0.125, 0.125, 0.125},
                                  // { 0.125,-0.125,-0.125, 0.125,-0.125, 0.125, 0.125,-0.125},
                                  // { 0.125,-0.125, 0.125,-0.125, 0.125,-0.125, 0.125,-0.125},
                                  // {-0.125, 0.125,-0.125, 0.125, 0.125,-0.125, 0.125,-0.125}};   
  } else if (m_dim == 2){
    // double Sig[4][8] ={{0.25, -0.25, 0.25,-0.25},
    // {0.25, -0.25, 0.25,-0.25},
    // {0.25, -0.25, 0.25,-0.25},
    // {0.25, -0.25, 0.25,-0.25}};    
  }

  //double f = 1/8;
  par_loop(e,m_elem_count){   
        
  if (m_gp_count==1){
      int offset = e * m_nodxelem*m_dim;   //SCALAR  
      //double hmod[m_dim][4];
      
    for (int d=0;d<m_dim;d++)
        for (int n=0;n<m_nodxelem;n++)
          m_f_elem_hg[offset + n*m_dim + d] = 0.0;
  

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
      double c_h = 0.1* pow(vol[e], 0.6666666) * rho[e] * 0.2500 * mat[e]->cs0;
      ////printf("c_h %.6e\n", c_h);

      for (int n=0;n<m_nodxelem;n++){      
        for (int d=0;d<m_dim;d++){
          m_f_elem_hg[offset + n*m_dim + d] *= c_h;
        ////printf("hg forces: %f %f %f\n",m_f_elem_hg[offset + n*m_dim],m_f_elem_hg[offset + n*m_dim + 1],m_f_elem_hg[offset + n*m_dim + 2]  );
        }
      }
      
  } //gp ==1
  }//ELEM
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

 
}; //Namespace MetFEM
