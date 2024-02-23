///////////////////////////////////////////////////
//sy = [A + B(epl^n)] [1 + C ln(e_dot pl/e_dot 0) (1 - pow)]

inline double __device__ JohnsonCook::CalcYieldStress(const double &strain, const double &strain_rate, const double &temp)	{
	double T_h = (temp - T_t) / (T_m - T_t);
	double sr = strain_rate;
	if (strain_rate == 0.0)
		sr = 1.e-5;
	
	double sy = (A+B*pow(strain, n))*(1.0 + C * log (sr/ eps_0) ) * (1.0 - pow(T_h,m));
	
	return sy;
}	

#include<iostream>
using namespace std;
inline double __device__ JohnsonCook::CalcTangentModulus(const double &plstrain, const double &strain_rate, const double &temp)	{
	double sy, T_h;
  //cout << "n, B, C, eps_0, T_t, m"<<n<<", "<<B<<", "<<C<<"eps0, "<<eps_0<<", "<<", "<<T_t<<", "<<m<<endl;
	T_h = (temp - T_t) / (T_m - T_t);
	
  //double sy = (A+B*pow(strain, n))*(1.0 + C * log (strain_rate/ eps_0) ) * (1.0 - pow(T_h,m));
  double Et =0.;

  if (plstrain > 0.)
    Et = n * B * pow(plstrain,n-1.)*(1.0 + C*log(strain_rate/ eps_0)) * (1.0-pow (T_h,m));
  else 
    Et = Elastic().E()*0.1; //ARBITRARY! TODO: CHECK MATHEMATICALLY
  return Et;
}	

inline double __device__ Hollomon::CalcYieldStress(const double &strain){
  double sy = 0.0;
  printf("K %f eps0 %f , eps1 %f sy0 %f m %f\n",K,eps0,eps1,sy0,m);
   if (strain + eps0 > eps1) sy = K*pow(strain + eps0, m); //plateau surpassed. If no plateau, eps1=eps0 so 
   else                      sy = sy0; 
	return sy;
}	

inline double __device__ Hollomon::CalcTangentModulus(const double &strain) {
	double Et;
  if (strain + eps0 > eps1) Et = K*m*pow(strain + eps0, (m-1.0));
  else                      Et = 0.;
	//cout << "ET: "<<Et<<endl;
	return Et;
}

// IF NO VIRTUAL FUNCTIONS ARE USED 

// __device__ double CalcHollomonYieldStress(const double &strain, Material_ *mat){
  // double sy = 0.0;
  // //printf("K %f eps0 %f , eps1 %f sy0 %f m %f\n",K,eps0,eps1,sy0,m);
   // if (strain + mat->eps0 > mat->eps1) sy = mat->K*pow(strain + mat->eps0,mat->m); //plateau surpassed. If no plateau, eps1=eps0 so 
   // else                      sy = mat->sy0; 
	// return sy; 
  
// }