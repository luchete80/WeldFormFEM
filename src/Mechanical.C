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

        //!!! er hoop = vr/r
        //if (dim .eq. 2 .and. bind_dom_type .eq. 3) then 
        //  ! if (elem%gausspc(e) .eq. 1) then
        //    elem%str_rate(e,gp, 3,3) = elem%str_rate(e,gp, 3,3) + 0.25d0*elem%vele (e,dim*(n-1)+1,1) / elem%radius(e,gp) !! 0.25 is shapemat
        //  ! print *, "hoop er", elem%str_rate(e,gp, 3,3) 
        //  elem%rot_rate(e,gp, 3,3) = 0.0d0
        //  ! end if
        //end if 

                                                       getDerivative(e,gp,0,n) * getVElem(e,n,1)));
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
      //printf("Element %d Vol %f\n",e,vol[e]);
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
      for (int i=0;i<3;i++)
        for (int j=0;j<3;j++){
          //printf("SIGMA %d %d %.6e\n",i,j,getSigma(e,gp,i,j));
        }
      for (int n=0; n<m_nodxelem;n++) {
        for (int d=0;d<m_dim;d++){
          m_f_elem[offset + n*m_dim + d] += getDerivative(e,gp,d,n) * getSigma(e,gp,d,d);
        }
      
        if (m_dim == 2){
          if (m_domtype != _Axi_Symm_){
          m_f_elem[offset + n*m_dim    ] +=  getDerivative(e,gp,1,n) * getSigma(e,gp,0,1);
          m_f_elem[offset + n*m_dim + 1] +=  getDerivative(e,gp,0,n) * getSigma(e,gp,0,1);
          } else {
//              fa = 0.25d0/elem%radius(e,gp) * elem%detJ(e,gp) !!! THEN IS WEIGHTED BY 4 in case of gauss point =1
//              !!! AREA WEIGHTED, BENSON EQN 2.4.3.2
 //             !!! 2.4.3.2 remains sig * Area/(4 r0), which is (4detJ)/(4r0) = detJ /r0
//              !!! LATER IS MULTIPLIED BY WEIGHT WICH GIVES THE AREA
                double fa = 0.25 / m_radius[e] * m_detJ[offset_det + gp];; //TODO: CHANGE According to element data
//
                m_f_elem[offset + n*m_dim    ] += getDerivative(e,gp,1,n) * getSigma(e,gp,0,1) - 
                                                  getSigma(e,gp,0,0) - getSigma(e,gp,2,2);

                m_f_elem[offset + n*m_dim    ] += getDerivative(e,gp,0,n) * getSigma(e,gp,0,1) - 
                                                  getSigma(e,gp,0,1) *fa;
                                                  
//              elem%f_int(e,n,1) = elem%f_int(e,n,1) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 1,2) - &
//                                                     (elem%sigma (e,gp, 1,1) - elem%sigma (e,gp, 3,3) ) * fa
//                                                     
//              elem%f_int(e,n,2) = elem%f_int(e,n,2) + elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,2) - &
//                                                     elem%sigma (e,gp, 1,2) * fa       
            
          }
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

      
 
        
      }// nod x elem


    } // Gauss Point
    
    for (int n=0; n<m_nodxelem;n++) {
      for (int d=0;d<m_dim;d++){
        m_f_elem[offset + n*m_dim + d] *= w;
      }
    }  
    
    //TO CHECK
    /*
    for (int n=0; n<m_nodxelem;n++) {
      printf("Element %d forces\n",e);
      printf("%.3e %.3e %.3e\n",m_f_elem[offset + n*m_dim ],m_f_elem[offset + n*m_dim + 1] ,m_f_elem[offset + n*m_dim + 2] );
    } 
    */ 
    
  }//if e<elem_count
}

dev_t void Domain_d::calcElemPressure(){

  par_loop(e,m_elem_count){
    //printf("calc pressure \n");
    int offset_t = e * m_gp_count *6;

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
      //printf("pressure %f\n",p[offset + gp]);
    }

  } // e< elem_count
}

//Computational Methods in lagfangian & Eulerian Hydrocodes
//From Benson 1992
// Equation 1.3.12
dev_t void Domain_d::calcElemPressureFromJ(){

  par_loop(e,m_elem_count){
    p[e] = mat[e]->Elastic().BulkMod() * ( 1.0 - vol[e]/vol_0[e] );
  } // e< elem_count
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
  } //NODE LOOP

  par_loop(e,m_elem_count){
    for (int ne=0;ne<m_nodxelem;ne++) //I.E. 4 
      p[e] += pn[m_elnod[e * m_nodxelem+ne]];    
    p[e] *= 0.25;
  }
  delete pn,voln_0,voln;
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
  printf("Total vol %f\n",tot_vol);
}

//Assuming constant material
dev_t void Domain_d::CalcNodalMassFromVol(){
  double *rhon = new double [m_node_count];  
  double tot_mass = 0.0;
  
  par_loop(n, m_node_count){
    rhon[n]=0.0;
    for (int e=0; e<m_nodel_count[n];e++) {    
      int eglob   = m_nodel     [m_nodel_offset[n]+e]; //Element
      rhon[n] += rho[eglob]; 
    }
    m_mdiag[n] = rhon[n]/(double)m_nodel_count[n] * m_voln[n];
    tot_mass +=m_mdiag[n];
    //printf("Node %d mass %f rho %f vol %f\n",n,m_mdiag[n],rhon[n]/(double)m_nodel_count[n] , m_voln[n]);
    
  } //NODE LOOP
  printf("Tot mass: %f\n",tot_mass);
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


	
  // implicit none
  // real(fp_kind) :: SRT(3,3), RS(3,3), ident(3,3)
  // integer :: e,gp
  // real(fp_kind) ,intent(in):: dt
  
  // real(fp_kind) :: p00
  
  // p00 = 0.
  
  // ident = 0.0d0
  // ident (1,1) = 1.0d0; ident (2,2) = 1.0d0;  ident (3,3) = 1.0d0

  par_loop(e,m_elem_count){
      //printf("calculating sigma \n");
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
      //printf("mat G %f\n",mat[e]->Elastic().G());
      ShearStress	= ShearStress  + dt*(2.0* mat[e]->Elastic().G()*(StrRate - 1.0/3.0*Trace(StrRate) * Identity() ) + SRT+RS);
      
      double J2 = 0.5*(ShearStress.xx*ShearStress.xx +  2.0*ShearStress.xy*ShearStress.xy + 
                                      2.0*ShearStress.xz*ShearStress.xz + 
                     ShearStress.yy*ShearStress.yy+  
                                      2.0*ShearStress.yz*ShearStress.yz +               
                     ShearStress.zz*ShearStress.zz                 
                                     
                    );
      double sig_trial = sqrt(3.0*J2);
      //printf("sy %f\n", sigma_y[e]);
      //printf("Sigmay %.f\n",pl_strain[e]);
      if (sigma_y[e]<sig_trial){
        //printf("Yield elem %d, %f, sig_trial %f, yield stress %f\n",e,pl_strain[e],sig_trial, sigma_y[e]);
        //elem%shear_stress(e,gp,:,:) = elem%shear_stress(e,gp,:,:) * elem%sigma_y(e,gp) / sig_trial
        //elem%pl_strain(e,gp) = elem%pl_strain(e,gp) + (sig_trial - elem%sigma_y(e,gp)) / (3.0d0 * mat_G) !
        ShearStress = ShearStress * (sigma_y[e] / sig_trial);
        pl_strain[e] += (sig_trial - sigma_y[e]) / (3.0 *  mat[e]->Elastic().G());

/*
To calculate inelastic fraction and plastic work
double f = dep/Sigmay;
		Strain_pl_incr (0,0) = f*(Sigma(0,0)-0.5*(Sigma(1,1) + Sigma(2,2) ));
		Strain_pl_incr (1,1) = f*(Sigma(1,1)-0.5*(Sigma(0,0) + Sigma(2,2) ));
		Strain_pl_incr (2,2) = f*(Sigma(2,2)-0.5*(Sigma(0,0) + Sigma(1,1) ));
		Strain_pl_incr (0,1) = Strain_pl_incr (1,0) = 1.5*f*(Sigma(0,1));
		Strain_pl_incr (0,2) = Strain_pl_incr (2,0) = 1.5*f*(Sigma(0,2));
		Strain_pl_incr (1,2) = Strain_pl_incr (2,1) = 1.5*f*(Sigma(1,2));
		
		Strain_pl = Strain_pl + Strain_pl_incr;
*/

      }
      //printf("Shear Stress\n");
      //print(ShearStress);
      // elem%shear_stress(e,gp,:,:)	= dt * (2.0 * mat_G *(elem%str_rate(e,gp,:,:) - 1.0/3.0 * &
                                   // (elem%str_rate(e,gp,1,1)+elem%str_rate(e,gp,2,2)+elem%str_rate(e,gp,3,3))*ident) &
                                   // +SRT+RS) + elem%shear_stress(e,gp,:,:)
                                   
      // elem%sigma(e,gp,:,:) = -elem%pressure(e,gp) * ident + elem%shear_stress(e,gp,:,:)	!Fraser, eq 3.32
      Sigma = -p[offset_s] * Identity() + ShearStress;
      //printf("SHEAR STRESS\n");
      //print(ShearStress);

      // printf("STR RATE\n");
      // print(StrRate);
      
      //printf("ELEMENT SIGMA\n");
      //print(Sigma);
      double Ep = 0;
			double dep=( sig_trial - sigma_y[e])/ (3.*mat[e]->Elastic().G() + Ep);	//Fraser, Eq 3-49 TODO: MODIFY FOR TANGENT MODULUS = 0
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
      
    }//gp
  }//el < elcount
    // end do !gauss point
  // end do
  //printf("ELEMENT %d SIGMA\n");
 
}



dev_t void Domain_d:: calcElemHourglassForces()
{
  if (m_dim == 3 && m_nodxelem == 4)
    return;
    
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
      double c_h = 0.0 * pow(vol[e], 0.6666666) * rho[e] * 0.2500 * mat[e]->cs0;
      //printf("c_h %.6e\n", c_h);


      for (int n=0;n<m_nodxelem;n++){      
        for (int d=0;d<m_dim;d++){
          m_f_elem_hg[offset + n*m_dim + d] *= c_h;
        }
        //printf("hg forces el %d node %d: %f %f %f\n",e, n,m_f_elem_hg[offset + n*m_dim],m_f_elem_hg[offset + n*m_dim + 1],m_f_elem_hg[offset + n*m_dim + 2]  );
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
