#ifndef _MATERIAL_CUH_
#define _MATERIAL_CUH_

#include "defs.h"
#include "Matrix.h"


//TODO: CHANGE TO ENUM
#define BILINEAR				0
#define HOLLOMON				1 //POWER LAW
#define JOHNSON_COOK		2
#define _GMT_         	3

class Elastic_{
	private:
	double E_m, nu_m;	//Poisson and young
	double K_m, G_m;
  double cs_m;
  double rho_m;
	
	public:
	 spec_ Elastic_(){}
	 spec_ Elastic_(const double &e, const double &nu)
                       :E_m(e),nu_m(nu){
                         K_m= e / ( 3.0*(1.0 -2.0*nu) );
                         G_m= e / (2.0* (1.0 + nu));
                        }
	 spec_ const double& Poisson()const{return nu_m;}
	 spec_ const double& E()const{return E_m;}
	 spec_ const double& G()const{return G_m;}
   spec_ const double& BulkMod()const{return K_m;}
	
};

class Material_{
	
	protected:
	Elastic_ elastic_m;

	double E_m, nu;	//TODO, move to elastic class

  
	public:
  double cs0;
  double Ep;  //If Bilinear this is constant, 
  int			Material_model;	//TODO: Change to enum

  ////// NO //virtual FUNCTIONS, HOLLOMON MATERIAL ///////
  double K, m;
	double eps0, eps1;
  double sy0;
  void InitHollomon(){}

  ////// NO //virtual FUNCTIONS, JOHNSON COOK MATERIAL ///////
	double T_t,T_m;	//transition and melting temps
  
  double T_min, T_max;
  double e_min, e_max;
  double er_min, er_max;
  
  //JOHNSON COOK
	double A, B, C;
	double n/*, m*/;
	double eps_0;
  
  //GMT
  double C1, C2, m1, m2, n1,n2, I1,I2;
  
  //THERMAL
  double k_T, cp_T; ///MAYBE MOVE TO element or nodes


  spec_ void InitHollomon(const Elastic_ &el, const double sy0_, const double &k_, const double &m_)
  {
    elastic_m = el;
    K = k_;
    m = m_;
  Material_model = HOLLOMON;
  eps0 = sy0_/el.E(); 
  sy0  = sy0_;
  eps1 = pow(sy0_/k_, 1./m);
  printf( "eps_0, %.2e eps_1 %.2e\n",eps0, eps1);
  if (eps0 > eps1){
    printf ("ERROR, Hollomon material bad definition, please correct Yield Stress, Elastic Modulus or Material hardening constants.");
  }
  
  }
  void Init_JohnsonCook(const Elastic_ &el,const double &a, const double &b, const double &n_, 
              const double &c, const double &eps_0_,
              const double &m_, const double &T_m_, const double &T_t_)
  //:
	// Material_(el),A(a),B(b),C(c),
  // m(m_),n(n_),eps_0(eps_0_),T_m(T_m_),T_t(T_t_)
  { 
    elastic_m = el;
    Material_model = JOHNSON_COOK;
    A=a;B=b;C=c;
    m=m_;n=n_;eps_0=eps_0_;
    T_m=T_m_;T_t=T_t_;
  }

	Material_(){
    sy0=1.0e10;
  }
  spec_ void test(){printf("test\n");}
//  //virtual spec_ double testret(){return 2.0;}
	Material_(const Elastic_ el):elastic_m(el){
     sy0=1.0e10;
  }
	// //virtual  __device__ inline double CalcTangentModulus(){};
	// //virtual  __device__ inline double CalcTangentModulus(const double &strain, const double &strain_rate, const double &temp){};
	// //virtual  __device__ inline double CalcTangentModulus(const double &strain){};
	// //virtual  __device__ inline double CalcYieldStress(){};
	// //virtual  __device__ inline double CalcYieldStress(const double &strain){return 0.0;};
	// //virtual  __device__ inline double CalcYieldStress(const double &strain, const double &strain_rate, const double &temp){};
	 spec_ const Elastic_& Elastic()const{return elastic_m;}

  Matrix getElasticMatrix(){
      Matrix D(6,6);
      double E = Elastic().E();
      double nu = Elastic().Poisson();
      double G = Elastic().G();
      
      double lambda  = (E * nu) /((1.0+nu)*(1.0-2.0*nu)); 
      D.Set(0,1, lambda);               D.Set(0,2, lambda);
      D.Set(1,0, lambda);               D.Set(1,2, lambda);
      D.Set(2,0, lambda);               D.Set(2,1, lambda);
      
      printf("Cmat\n");
      
      for (int d=0;d<3;d++) D.Set(d,d,lambda+2.0*G);
      for (int d=3;d<6;d++) D.Set(d,d,G);    
      return D;
  }

}; //MATERIAL 

class _Plastic{
	
	public:
	//virtual inline double CalcYieldStress();	
	////virtual inline double CalcYieldStress();
};

class Bilinear:
public Material_{

 	public:
	Bilinear(const double &ep){ //THIS IS DIFFERENT FROM WELDFORM CPU, IN WHICH BILINEAR IS NOT A MATERIAL
    Ep = ep;
    Material_model = BILINEAR;
  }
};


//TODO: derive johnson cook as plastic material flow
class JohnsonCook:
public Material_{
	double T_t,T_m;	//transition and melting temps
	double A, B, C;
	double n, m;
	double eps_0;
	
	public:
	JohnsonCook(){
    Material_model = JOHNSON_COOK;
  }
	//You provide the values of A, B, n, m, 
	//θmelt, and  θ_transition
	//as part of the metal plasticity material definition.
	JohnsonCook(const Elastic_ &el,const double &a, const double &b, const double &n_, 
              const double &c, const double &eps_0_,
              const double &m_, const double &T_m_, const double &T_t_):
	Material_(el),A(a),B(b),C(c),
  m(m_),n(n_),eps_0(eps_0_),T_m(T_m_),T_t(T_t_)
  { Material_model = JOHNSON_COOK;
  }
	inline double dev_t CalcYieldStress(){return 0.0;}	
	inline double dev_t CalcYieldStress(const double &plstrain){
     double Et =0.;

    if (plstrain > 0.)
      Et = n * B * pow(plstrain,n-1.);
    else 
      Et = Elastic().E()*0.1; //ARBITRARY! TODO: CHECK MATHEMATICALLY
    return Et;
  } //TODO: SEE IF INCLUDE	
	inline double  dev_t CalcYieldStress(const double &strain, const double &strain_rate, const double &temp);	
	inline double  dev_t CalcTangentModulus(const double &strain, const double &strain_rate, const double &temp);
  //~JohnsonCook(){}
};

class Hollomon:
public Material_{
	double K, m;
	double eps0, eps1;
  double sy0;
	
	public:
	spec_ Hollomon(){Material_model = HOLLOMON;}
	//You provide the values of A, B, n, m, 
	//θmelt, and  θ_transition
	//as part of the metal plasticity material definition.
	//ASSUMING AT FIRST COEFFICIENTS ARE GIVEN TO TOTAL STRAIN-STRESS
	spec_ Hollomon(const double eps0_, const double &k_, const double &m_):
    K(k_), m(m_){ 
    eps0 = eps0_;
    Material_model = HOLLOMON;}
  spec_  Hollomon(const Elastic_ &el, const double sy0_, const double &k_, const double &m_):
  Material_(el),K(k_), m(m_) {
  Material_model = HOLLOMON;
  eps0 = sy0_/el.E(); 
  sy0  = sy0_;
  eps1 = pow(sy0_/k_, 1./m);
  printf( "eps_0, %.2e eps_1 %.2e\n",eps0, eps1);
  if (eps0 > eps1){
    printf ("ERROR, Hollomon material bad definition, please correct Yield Stress, Elastic Modulus or Material hardening constants.");
  }
  
  }
  spec_ double testret(){printf("hollomon testret\n"); return 2.0;}
	
  inline double  dev_t CalcTangentModulus(const double &strain);
	inline double  dev_t CalcYieldStress(){}	
	inline double  dev_t CalcYieldStress(const double &strain);	
};


//TODO: derive johnson cook as plastic material flow
class GMT:
public Material_{
	double C1, C2;
  double n1, n2; //Strain hardening exponent
	double m1, m2;
  double I1, I2; //EXPONENTIAL TERMS
	double eps_0; //ONLY FOR JC DAMAGE , CORRECT THIS
  
	
	public:
	GMT(){
    Material_model = _GMT_;
  }
	//You provide the values of A, B, n, m, 
	//θmelt, and  θ_transition
	//as part of the metal plasticity material definition.
	GMT(const Elastic_ &el, 
              const double &n1_,  const double &n2_,  
              const double &C1_,  const double &C2_, 
              const double &m1_,  const double &m2_,
              const double &I1_,  const double &I2_,
              const double &e_min_  = 0.0,const double &e_max_  =1.0e10, /*const double &e_0 = 1.0,*/
              const double &er_min_ = 0.0,const double &er_max_ =1.0e10,
              const double &T_min_  = 0.0,const double &T_max_  =1.0e10):
	Material_(el),
  C1(C1_),C2(C2_),
  n1(n1_),n2(n2_),
  m1(m1_),m2(m2_),
  I1(I1_),I2(I2_)
  {
    e_min =e_min_; e_max =e_max_;
    er_min=er_min_;er_max=er_max_;
    T_min =T_min_; T_max =T_max_;
		// T_m=T_m_;
		// T_t=T_t_;
	}
	inline double CalcYieldStress(){return 0.0;}	
	inline double CalcYieldStress(const double &plstrain){
     // double Et =0.;

    // if (plstrain > 0.)
      // Et = n * B * pow(plstrain,n-1.);
    // else 
      // Et = Elastic().E()*0.1; //ARBITRARY! TODO: CHECK MATHEMATICALLY
    // return Et;
  } //TODO: SEE IF INCLUDE	
	inline double CalcYieldStress(const double &strain, const double &strain_rate, const double &temp);	
	inline double CalcTangentModulus(const double &strain, const double &strain_rate, const double &temp);
  double &getRefStrainRate(){return eps_0;}//only for JC
  //~JohnsonCook(){}
};


/////// VIRTUAL FUNCTIONS NOT SUPPORTED ON CUDA ///////

dev_t inline double CalcHollomonYieldStress(const double &strain, Material_ *mat) //IN CASE OF NOT USING //virtual FUNCTIONS
{
  double sy = 0.0;
  //printf("K %f eps0 %f , eps1 %f sy0 %f m %f\n",mat->K,mat->eps0,mat->eps1,mat->sy0,mat->m);
   if (strain + mat->eps0 > mat->eps1) sy = mat->K*pow(strain + mat->eps0,mat->m); //plateau surpassed. If no plateau, eps1=eps0 so 
   //if (strain + mat->eps0 > mat->eps1) sy = mat->K*pow(strain /*+ mat->eps0*/,mat->m); //plateau surpassed. If no plateau, eps1=eps0 so 
   //if (strain  > mat->eps0) sy = mat->K*pow(strain /*+ mat->eps0*/,mat->m); //plateau surpassed. If no plateau, eps1=eps0 so 
   else                      sy = mat->sy0; 
  //if (sy>mat->sy0)printf("SY %.3f\n", sy);
	return sy; 
  
}

dev_t inline void ShowProps(Material_ *mat) //IN CASE OF NOT USING //virtual FUNCTIONS
{
  printf("K %f \n eps0 %f \n eps1 \n%f, sy0 \n%f, m %f \n", mat->K,mat->eps0,mat->eps1,mat->sy0,mat->m);
  
}

dev_t inline double CalcJohnsonCookYieldStress(const double &strain, const double &strain_rate, const double &temp, Material_ *mat) //IN CASE OF NOT USING //virtual FUNCTIONS
{
	double T_h = (temp - mat->T_t) / (mat->T_m - mat->T_t);
	double sr = strain_rate;
	if (strain_rate == 0.0)
		sr = 1.e-5;
	
	double sy = (mat->A+mat->B*pow(strain, mat->n))*(1.0 + mat->C * log (sr/ mat->eps_0) ) * (1.0 - pow(T_h,mat->m));
	
	return sy;
}

inline double dev_t CalcHollomonTangentModulus(const double &strain, Material_ *mat) {
	double Et;
  if (strain + mat->eps0 > mat->eps1) Et = mat->K * mat->m * pow(strain + mat->eps0, (mat->m-1.0));
  else                      Et = 0.;
	//cout << "ET: "<<Et<<endl;
	return Et;
}

inline double dev_t CalcJohnsonCookTangentModulus(const double &plstrain, const double &strain_rate, const double &temp, Material_ *mat)	{
	double sy, T_h;
  //cout << "n, B, C, eps_0, T_t, m"<<n<<", "<<B<<", "<<C<<"eps0, "<<eps_0<<", "<<", "<<T_t<<", "<<m<<endl;
	T_h = (temp - mat->T_t) / (mat->T_m - mat->T_t);
	
  //double sy = (A+B*pow(strain, n))*(1.0 + C * log (strain_rate/ eps_0) ) * (1.0 - pow(T_h,m));
  double Et =0.;

  if (plstrain > 0.)
    Et = mat->n * mat->B * pow(plstrain,mat->n-1.)*(1.0 + mat->C*log(strain_rate/ mat->eps_0)) * (1.0-pow (T_h,mat->m));
  else 
    Et = mat->Elastic().E()*0.1; //ARBITRARY! TODO: CHECK MATHEMATICALLY
  return Et;
}	


///////////////////////////////////////////////////
//sy = [A + B(epl^n)] [1 + C ln(e_dot pl/e_dot 0) (1 - pow)]
////////////// TEMPERATURE SHOULD BE IN CELSIUS

inline double dev_t CalcGMTYieldStress(const double &strain, const double &strain_rate, const double &temp, Material_ *mat)	{
  // OLD////////////////////////
	// double T_h = (temp - T_t) / (T_m - T_t);
	// double sr = strain_rate;
	// if (strain_rate == 0.0)
		// sr = 1.e-5;
    
	// double sy = (A+B*pow(strain, n))*(1.0 + C * log (sr/ eps_0) ) * (1.0 - pow(T_h,m));  
  
  // NEW /////////////////////
	double T_h = (temp - mat->T_t) / (mat->T_m - mat->T_t);
	double sr = strain_rate;
  double f = 1.0;
	// if (strain_rate > eps_0)
		// f = (1.0 + C * log(strain_rate/eps_0));
	
	double e,er,T, sy;
  e = strain; er = strain_rate; T = temp;
  if      (e < mat->e_min) e = mat->e_min;
  else if (e > mat->e_max) e = mat->e_max;

  if      (er < mat->er_min) er = mat->er_min;
  else if (er > mat->er_max) er = mat->er_max;
  
  if      (T < mat->T_min) T = mat->T_min;
  else if (T > mat->T_max) T = mat->T_max;
  
  sy = mat->C1 * exp(mat->C2*T)*pow(e,mat->n1*T + mat->n2) * exp((mat->I1*T + mat->I2)/e) * pow(er,mat->m1 * T + mat->m2);
  // cout << "n1 n2 C1 c2 m1 m2 I1 i2" << n1 << ", "<<n2 << ", "
                                    // << C1 << ", "<<C2 << ", "
                                    // << m1 << ", "<<m2 << ", "
                                    // << I1 << ", "<<I2 << ", "<<endl;
  // // cout << "e, er, T " << e <<", "<<er<<", "<<T<<endl; 
	// cout << "sy GMT"<< sy<<", e, er, T" << e << ", " <<er << ", " <<T<<endl;
	return sy;
}	

inline double dev_t CalcGMTTangentModulus(const double &plstrain, const double &strain_rate, const double &temp, Material_ *mat)	{
	double T_h;
  
	double e,er,T, sy;
  e = plstrain; er = strain_rate; T = temp;
  if      (plstrain < mat->e_min) e = mat->e_min;
  else if (plstrain > mat->e_max) e = mat->e_max;

  if      (strain_rate < mat->er_min) er = mat->er_min;
  else if (strain_rate > mat->er_max) er = mat->er_max;

  if      (temp < mat->T_min) T = mat->T_min;
  else if (temp > mat->T_max) T = mat->T_max;
  
  //double sy = (A+B*pow(strain, n))*(1.0 + C * log (strain_rate/ eps_0) ) * (1.0 - pow(T_h,m));
  double Et =0.;

	Et = mat->C1*exp(mat->C2*T) * pow(er,mat->m1*T + mat->m2)* //constant part
       pow(e,T*mat->n1 + mat->n2-2.0)*(-mat->I1*T - mat->I2 + e*(mat->n1*T + mat->n2))*exp((mat->I1*T + mat->I2)/e);

  // if (strain_rate > eps_0)
		// f = (1.0 + C * log(strain_rate/eps_0));
   // if (plstrain > 0.)
    // Et = n * B * pow(plstrain,n-1.) * f * (1.0-pow (T_h,m));
   // else 
     // //Et = Elastic().E()*0.1; //ARBITRARY! TODO: CHECK MATHEMATICALLY
		// Et = Elastic().E();
	 
  return Et;
}	



#endif
