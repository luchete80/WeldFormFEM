/*  Copyright (c) 2013-2019 INGV, EDF, UniCT, JHU

    Istituto Nazionale di Geofisica e Vulcanologia, Sezione di Catania, Italy
    Électricité de France, Paris, France
    Università di Catania, Catania, Italy
    Johns Hopkins University, Baltimore (MD), USA

    This file is part of GPUSPH. Project founders:
        Alexis Hérault, Giuseppe Bilotta, Robert A. Dalrymple,
        Eugenio Rustico, Ciro Del Negro
    For a full list of authors and project partners, consult the logs
    and the project website <https://www.gpusph.org>

    GPUSPH is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GPUSPH is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GPUSPH.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef TENSOR_IMPL
#define TENSOR_IMPL

#include "tensor.cuh"
#include "vector_math.h"
//#define __spec __device__ __forceinline__
#define __spec __device__ __inline__


#include <stdio.h> //Debug


//Converts tensor to flat symm
// __device__ inline symtensor3 FromFlatSym(double flat[]){
	// symtensor3 ret;
	// // for (int i=0;i<3;i++)
		// // ret.m_data [i][i] = flat[i];
	// // ret.m_data .xy = ret.m_data [1][0] = flat[3]; 
	// // ret.m_data [1][2] = ret.m_data [2][1] = flat[4]; 
	// // ret.m_data .xz = ret.m_data [2][0] = flat[5]; 
	
	// return ret;
// }

// __device__ tensor3::tensor3(double flat[]){	//Six components
	// // for (int i=0;i<3;i++)
		// // m_data [i][i] = flat[i];
	// //Check
	// // if (sizeof(flat/sizeof(double)==6)
		// // FromFlatSym(flat);
// }


//Converts tensor to flat symm
//TODO: CHECK IF THE RETURN VALUE IS IN PARAMETERS THIS IS FASTER
__spec tensor3 FromFlatSym(double flat[]){
	tensor3 ret;
  ret.xx = flat[0];		ret.yy = flat[1];		ret.zz = flat[2];
  ret.xy = ret.yx = flat[3];
  ret.yz = ret.zy = flat[4];
  ret.xz = ret.zx = flat[5];  
	
	return ret;
}

__spec tensor3 FromFlatAntiSym(double flat[]){
	tensor3 ret;
  ret.xx = flat[0];		ret.yy = flat[1];		ret.zz = flat[2];  

	ret.xy = flat[3]; 
	ret.yz = flat[4]; 
	ret.xz = flat[5]; 
	
	ret.yx = -flat[3]; 
	ret.zy = -flat[4]; 
	ret.zx = -flat[5]; 
  
  return ret;
}

// __device__ void tensor3::FromFlatSymPtr(double *flat){
	// for (int i=0;i<3;i++)
		// m_data [i][i] = flat[i];
	// m_data .xy = m_data [1][0] = flat[3]; 
	// m_data [1][2] = m_data [2][1] = flat[4]; 
	// m_data .xz = m_data [2][0] = flat[5]; 
// }

__spec void ToFlatSymPtr(const tensor3 &m_data, double *flat, int initial){
	flat [initial + 0] = m_data.xx; flat [initial + 1] = m_data .yy; flat [initial + 2] = m_data.zz;
  flat [initial + 3] = m_data.xy;
  flat [initial + 4] = m_data.yz;
  flat [initial + 5] = m_data.xz;
}

__spec tensor3 Identity(){
	tensor3 ret;
	ret.xx = ret.yy = ret.zz = 1.;
	ret.xy = ret.xz = ret.yx = ret.yz = ret.zx = ret.zy = 0.;
	//ret[1][1]=ret[2][2]=1.;
	
	return ret;
}

__spec
void
clear(tensor3& T)
{
	T.xx = T.xy = T.xz =
	T.yx = T.yy = T.yz = 
	T.zx = T.zy = T.zz = 0.0f;
}

/**** Methods for loading/storing tensors from textures and array ****/

// // // #include "textures.cuh"

// // // //! Fetch tau tensor from texture
// // // /*!
 // // // an auxiliary function that fetches the tau tensor
 // // // for particle i from the textures where it's stored
// // // */
// // // __device__
// // // symtensor3 fetchTau(uint i)
// // // {
	// // // symtensor3 tau;
	// // // double2 temp = tex1Dfetch(tau0Tex, i);
	// // // tau.xx = temp.x;
	// // // tau.xy = temp.y;
	// // // temp = tex1Dfetch(tau1Tex, i);
	// // // tau.xz = temp.x;
	// // // tau.yy = temp.y;
	// // // temp = tex1Dfetch(tau2Tex, i);
	// // // tau.yz = temp.x;
	// // // tau.zz = temp.y;
	// // // return tau;
// // // }

//! Fetch tau tensor from split arrays
/*!
 an auxiliary function that fetches the tau tensor
 for particle i from the arrays where it's stored
*/

// __device__ double& tensor3::operator()(int row, int col)
// {
    // // assert(col >= 0 && col < 4);
    // // assert(row >= 0 && row < 4);

    // return m_data[row][col];
// }

 // __device__ double  tensor3::operator()(int row, int col) const
// {
    // // assert(col >= 0 && col < 4);
    // // assert(row >= 0 && row < 4);

    // return m_data[row][col];
// }

// __device__ void tensor3::operator()() {
    // // reset all elements of the tensor3 to 0.0
    // for (int row{ 0 }; row < 4; ++row) {
        // for (int col{ 0 }; col < 4; ++col) {
            // m_data[row][col] = 0.0;
        // }
    // }
// }

// // __device__ tensor3::tensor3(){
	// // data = new double[3];
	// // for (int i=0;i<9;i++) data[i] = 0.;
// // }



// __device__ double& tensor3::operator[](const int &i){

	// //return m_data[3*i];
// }

// __device__ tensor3 tensor3:: operator* (const tensor3 &b){
	// tensor3 ret;
	// for (int i=0;i<3;i++)
		// for (int j=0;j<3;j++)
			// for (int k=0;k<3;k++)
				// ret.m_data[i][j]+=m_data[i][k]*b.m_data[k][j];
	
	// return ret;
// }

// __device__ tensor3 tensor3:: operator+ (const tensor3 &b){
	// tensor3 ret;
	// for (int i=0;i<3;i++)
		// for (int j=0;j<3;j++)
				// ret.m_data[i][j]=m_data[i][j]+b(i,j);
	
	// return ret;
// }

// __device__ tensor3&  tensor3::operator+= (const tensor3 &b){
	// //tensor3 ret;
	// for (int i=0;i<3;i++)
		// for (int j=0;j<3;j++)
				// m_data[i][j]=m_data[i][j]+b(i,j);
	// // ret = this;
	// return *this;
// }

// __device__ tensor3 tensor3:: operator- (const tensor3 &b){
	// tensor3 ret;
	// for (int i=0;i<3;i++)
		// for (int j=0;j<3;j++)
				// ret.m_data[i][j]=m_data[i][j]-b.m_data[i][j];
	
	// return ret;
// }

__spec double3 operator* (const tensor3 &m_data, const double3 &v){
	double3 ret;
	ret.x = m_data.xx*v.x + m_data.xy * v.y+m_data.xz*v.z;	
	ret.y = m_data.yx*v.x + m_data.yy * v.y+m_data.yz*v.z;	
	ret.z = m_data.zx*v.x + m_data.zy * v.y+m_data.zz*v.z;	
	return ret;
}

__spec double3 operator* (const double3 &v, const tensor3 &m_data){
	double3 ret;
	ret.x = m_data.xx*v.x+m_data.xy*v.y+m_data.xz*v.z;	
	ret.y = m_data.yx*v.x+m_data.yy*v.y+m_data.yz*v.z;	
	ret.z = m_data.zx*v.x+m_data.zy*v.y+m_data.zz*v.z;	
	return ret;
}

__spec tensor3 operator* (const tensor3 &a, const tensor3 &b){
	tensor3 ret;
	ret.xx = a.xx*b.xx + a.xy * b.yx + a.xz * b.zx;	
	ret.xy = a.xx*b.yx + a.xy * b.yy + a.xz * b.yz;
	ret.xz = a.xx*b.zx + a.xy * b.zy + a.xz * b.zz;
	
	ret.yx = a.yx*b.xx + a.yy * b.yx + a.yz * b.zx;	
	ret.yy = a.yx*b.yx + a.yy * b.yy + a.yz * b.yz;
	ret.yz = a.yx*b.zx + a.yy * b.zy + a.yz * b.zz;
	
	ret.zx = a.zx*b.xx + a.zy * b.yx + a.zz * b.zx;	
	ret.zy = a.zx*b.yx + a.zy * b.yy + a.zz * b.yz;
	ret.zz = a.zx*b.zx + a.zy * b.zy + a.zz * b.zz;
	return ret;
}


// __device__ tensor3 tensor3::operator*= (const double &v){
	// tensor3 ret;
	// for (int i=0;i<3;i++)
		// for (int j=0;j<3;j++){
				// ret.m_data[i][j]=m_data[i][j]*v;
				// m_data[i][j] = ret.m_data[i][j];
		// }
	// return ret;
// }

// __device__ tensor3& tensor3::operator* (const double &f){
	// for (int i=0;i<3;i++)
		// for (int j=0;j<3;j++)
			// m_data[i][j] *= f;	
	// return *this;
// }


// __device__ tensor3 tensor3:: operator- (const double &f){
	// tensor3 ret;
	// for (int i=0;i<3;i++)
		// for (int j=0;j<3;j++)
				// ret(i,j)=m_data[i][j]-f;
	
	// return ret;
// }

__device__ tensor3 Trans (const tensor3 &m_data){
	tensor3 ret;
	ret.xx = m_data.xx; 	ret.yy = m_data.yy; 	ret.zz = m_data.zz;
	ret.xy = m_data.yx; 	ret.xz = m_data.zx; 	
  ret.yx = m_data.xy;   ret.yz = m_data.zy; 	
  ret.zx = m_data.xz; 	ret.zy = m_data.yz;	
	return ret;
}

__spec tensor3 operator* (const double &f, const tensor3 &m_data){
	tensor3 ret;
	ret.xx = m_data.xx * f;	ret.xy = m_data.xy * f;	 ret.xz = m_data.xz * f;	
	ret.yx = m_data.yx * f;	ret.yy = m_data.yy * f;	 ret.yz = m_data.yz * f;	
	ret.zx = m_data.zx * f;	ret.zy = m_data.zy * f;	 ret.zz = m_data.zz * f;		
	return ret;
}

__spec tensor3 operator* (const tensor3 &m_data, const double &f){
	tensor3 ret;
	ret.xx = m_data.xx * f;	ret.xy = m_data.xy * f;	 ret.xz = m_data.xz * f;	
	ret.yx = m_data.yx * f;	ret.yy = m_data.yy * f;	 ret.yz = m_data.yz * f;	
	ret.zx = m_data.zx * f;	ret.zy = m_data.zy * f;	 ret.zz = m_data.zz * f;		
	return ret;
}

__spec
tensor3
operator +(tensor3 const& T1, tensor3 const& T2)
{
	tensor3 R;
	R.xx = T1.xx + T2.xx;
	R.xy = T1.xy + T2.xy;
	R.xz = T1.xz + T2.xz;
  
	R.yx = T1.yx + T2.yx;
	R.yy = T1.yy + T2.yy;
	R.yz = T1.yz + T2.yz;
	
  R.zx = T1.zx + T2.zx;
  R.zy = T1.zy + T2.zy;
  R.zz = T1.zz + T2.zz;
	return R;
}

__spec
tensor3
operator -(tensor3 const& T1, tensor3 const& T2)
{
	tensor3 R;
	R.xx = T1.xx - T2.xx;
	R.xy = T1.xy - T2.xy;
	R.xz = T1.xz - T2.xz;
  
	R.yx = T1.yx - T2.yx;
	R.yy = T1.yy - T2.yy;
	R.yz = T1.yz - T2.yz;
	
  R.zx = T1.zx - T2.zx;
  R.zy = T1.zy - T2.zy;
  R.zz = T1.zz - T2.zz;
	return R;
}

// __device__ tensor3 tensor3::operator= (const double &f){
	// tensor3 ret;
	// for (int i=0;i<3;i++)
		// for (int j=0;j<3;j++)
			// m_data[i][j] = f;
	// return ret;
// }

// // __device__ void tensor3::operator= (tensor3 &t){
		// // for (int i=0;i<3;i++)
		// // for (int j=0;j<3;j++)
			// // m_data[i][j] = t.m_data[i][j];
// // }


// __device__ /*__forceinline__*/ tensor3 Identity(){
	// tensor3 ret;
	// ret(0,0) = ret(1,1) = ret(2,2) = 1.;
	// //ret[1][1]=ret[2][2]=1.;
	
	// return ret;
// }

// /*__spec*/
// __device__
// double3
// dot(tensor3 const& T, double3 const& v)
// {
	// return make_double3(
			// T(0,0)*v.x + T(0,1)*v.y +T(0,2)*v.z,
			// T(1,0)*v.x + T(1,1)*v.y +T(1,2)*v.z,			
			// T(2,0)*v.x + T(2,1)*v.y +T(2,2)*v.z
			// );

// }

// __device__ tensor3 operator/ (const tensor3 &b, const double &f){

	// tensor3 ret;
	// for (int i=0;i<3;i++)
		// for (int j=0;j<3;j++)
			// ret(i,j) = b(i,j)/f;	
			
// }

__spec void print(const tensor3 &m_data){
	printf("[%f %f %f],\n[%f %f %f],\n[%f %f %f].\n",
	m_data.xx,m_data.xy,m_data.xz,
	m_data.yx,m_data.yy,m_data.yz,
	m_data.zx,m_data.zy,m_data.zz);
}

#undef __spec

#endif
