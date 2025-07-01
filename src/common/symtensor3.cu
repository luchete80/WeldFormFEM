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
#define __spec __device__ __forceinline__

#include <stdio.h> //Debug

__spec
void
clear(symtensor3& T)
{
	T.xx = T.xy = T.xz =
		T.yy = T.yz = T.zz = 0.0f;
}

__spec
void
clear(symtensor4& T)
{
	T.xx = T.xy = T.xz = T.xw =
		T.yy = T.yz = T.yw =
		T.zz = T.zw = T.ww = 0.0f;
}


// determinant of a 3x3 symmetric tensor
__spec
double
det(symtensor3 const& T)
{
	double ret = 0;
	ret += T.xx*(T.yy*T.zz - T.yz*T.yz);
	ret -= T.xy*(T.xy*T.zz - T.xz*T.yz);
	ret += T.xz*(T.xy*T.yz - T.xz*T.yy);
	return ret;
}

// determinant of a 4x4 symmetric tensor
__spec
double
det(symtensor4 const& T)
{
	double ret = 0;

	// first minor: ww * (xyz × xyz)
	double M = 0;
	M += T.xx*(T.yy*T.zz - T.yz*T.yz);
	M -= T.xy*(T.xy*T.zz - T.xz*T.yz);
	M += T.xz*(T.xy*T.yz - T.xz*T.yy);
	ret += M*T.ww;

	// second minor: -zw * (xyz × xyw)
	M = 0;
	M += T.xx*(T.yy*T.zw - T.yz*T.yw);
	M -= T.xy*(T.xy*T.zw - T.xz*T.yw);
	M += T.xw*(T.xy*T.yz - T.xz*T.yy);
	ret -= M*T.zw;

	// third minor: yw * (xyz × xzw)
	M = 0;
	M += T.xx*(T.yz*T.zw - T.zz*T.yw);
	M -= T.xz*(T.xy*T.zw - T.xz*T.yw);
	M += T.xw*(T.xy*T.zz - T.xz*T.yz);
	ret += M*T.yw;

	// last minor: xw * (xyz × yzw)
	M = 0;
	M += T.xy*(T.yz*T.zw - T.zz*T.yw);
	M -= T.xz*(T.yy*T.zw - T.yz*T.yw);
	M += T.xw*(T.yy*T.zz - T.yz*T.yz);
	ret -= M*T.xw;

	return ret;
}

// L-infinity norm of a symmetric 4x4 tensor
__spec
double
norm_inf(symtensor4 const& T)
{
	double m = fmaxf(T.xx, T.xy);
	m = fmaxf(m, T.xz);
	m = fmaxf(m, T.xw);
	m = fmaxf(m, T.yy);
	m = fmaxf(m, T.yz);
	m = fmaxf(m, T.yw);
	m = fmaxf(m, T.zz);
	m = fmaxf(m, T.zw);
	m = fmaxf(m, T.ww);
	return m;
}

__spec
symtensor3
inverse(symtensor3 const& T)
{
	symtensor3 R;
	double D(det(T));
	R.xx = (T.yy*T.zz - T.yz*T.yz)/D;
	R.xy = (T.xz*T.yz - T.xy*T.zz)/D;
	R.xz = (T.xy*T.yz - T.xz*T.yy)/D;
	R.yy = (T.xx*T.zz - T.xz*T.xz)/D;
	R.yz = (T.xz*T.xy - T.xx*T.yz)/D;
	R.zz = (T.xx*T.yy - T.xy*T.xy)/D;

	return R;
}

__spec
symtensor3
operator -(symtensor3 const& T1, symtensor3 const& T2)
{
	symtensor3 R;
	R.xx = T1.xx - T2.xx;
	R.xy = T1.xy - T2.xy;
	R.xz = T1.xz - T2.xz;
	R.yy = T1.yy - T2.yy;
	R.yz = T1.yz - T2.yz;
	R.zz = T1.zz - T2.zz;
	return R;
}

__spec
symtensor3
operator +(symtensor3 const& T1, symtensor3 const& T2)
{
	symtensor3 R;
	R.xx = T1.xx + T2.xx;
	R.xy = T1.xy + T2.xy;
	R.xz = T1.xz + T2.xz;
	R.yy = T1.yy + T2.yy;
	R.yz = T1.yz + T2.yz;
	R.zz = T1.zz + T2.zz;
	return R;
}

__spec
symtensor3 &
operator -=(symtensor3 &T1, symtensor3 const& T2)
{
	T1.xx -= T2.xx;
	T1.xy -= T2.xy;
	T1.xz -= T2.xz;
	T1.yy -= T2.yy;
	T1.yz -= T2.yz;
	T1.zz -= T2.zz;
	return T1;
}

__spec
symtensor3 &
operator +=(symtensor3 &T1, symtensor3 const& T2)
{
	T1.xx += T2.xx;
	T1.xy += T2.xy;
	T1.xz += T2.xz;
	T1.yy += T2.yy;
	T1.yz += T2.yz;
	T1.zz += T2.zz;
	return T1;
}

__spec
symtensor3
operator /(symtensor3 const& T1, double f)
{
	symtensor3 R;
	R.xx = T1.xx/f;
	R.xy = T1.xy/f;
	R.xz = T1.xz/f;
	R.yy = T1.yy/f;
	R.yz = T1.yz/f;
	R.zz = T1.zz/f;
	return R;
}


__spec
symtensor3 &
operator /=(symtensor3 &T1, double f)
{
	T1.xx /= f;
	T1.xy /= f;
	T1.xz /= f;
	T1.yy /= f;
	T1.yz /= f;
	T1.zz /= f;
	return T1;
}

__spec
double3
dot(symtensor3 const& T, double3 const& v)
{
	return make_double3(
			T.xx*v.x + T.xy*v.y + T.xz*v.z,
			T.xy*v.y + T.yy*v.y + T.yz*v.z,
			T.xz*v.x + T.yz*v.y + T.zz*v.z);

}

__spec
double3
dot(symtensor3 const& T, double4 const& v)
{
	return make_double3(
			T.xx*v.x + T.xy*v.y + T.xz*v.z,
			T.xy*v.y + T.yy*v.y + T.yz*v.z,
			T.xz*v.x + T.yz*v.y + T.zz*v.z);

}

// T.v
__spec
double4
dot(symtensor4 const& T, double4 const& v)
{
	return make_double4(
			T.xx*v.x + T.xy*v.y + T.xz*v.z + T.xw*v.w,
			T.xy*v.x + T.yy*v.y + T.yz*v.z + T.yw*v.w,
			T.xz*v.x + T.yz*v.y + T.zz*v.z + T.zw*v.w,
			T.xw*v.x + T.yw*v.y + T.zw*v.z + T.ww*v.w);

}

// v.T.w
__spec
double
dot(double4 const& v, symtensor4 const& T, double4 const& w)
{
	return dot(v, dot(T,w));
}

// v.T.v
__spec
double
ddot(symtensor4 const& T, double4 const& v)
{
	return T.xx*v.x*v.x + T.yy*v.y*v.y + T.zz*v.z*v.z + T.ww*v.w*v.w +
		2*(
			(T.xy*v.y + T.xw*v.w)*v.x +
			(T.yz*v.z + T.yw*v.w)*v.y +
			(T.xz*v.x + T.zw*v.w)*v.z);
}

// First row of the adjugate of a given tensor3
__spec
double4
adjugate_row1(symtensor4 const& T)
{
	return make_double4(
		T.yy*T.zz*T.ww + T.yz*T.zw*T.yw + T.yw*T.yz*T.zw - T.yy*T.zw*T.zw - T.yz*T.yz*T.ww - T.yw*T.zz*T.yw,
		T.xy*T.zw*T.zw + T.yz*T.xz*T.ww + T.yw*T.zz*T.xw - T.xy*T.zz*T.ww - T.yz*T.zw*T.xw - T.yw*T.xz*T.zw,
		T.xy*T.yz*T.ww + T.yy*T.zw*T.xw + T.yw*T.xz*T.yw - T.xy*T.zw*T.yw - T.yy*T.xz*T.ww - T.yw*T.yz*T.xw,
		T.xy*T.zz*T.yw + T.yy*T.xz*T.zw + T.yz*T.yz*T.xw - T.xy*T.yz*T.zw - T.yy*T.zz*T.xw - T.yz*T.xz*T.yw);
}

__spec tensor3 tensor3::Identity(){
	tensor3 ret;
	ret(0,0) = ret(1,1) = ret(2,2) = 1.;
	//ret[1][1]=ret[2][2]=1.;
	
	return ret;
}
#undef __spec

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
__device__
symtensor3 fetchTau(uint i,
	const double2 *__restrict__ tau0,
	const double2 *__restrict__ tau1,
	const double2 *__restrict__ tau2)
{
	symtensor3 tau;
	double2 temp = tau0[i];
	tau.xx = temp.x;
	tau.xy = temp.y;
	temp = tau1[i];
	tau.xz = temp.x;
	tau.yy = temp.y;
	temp = tau2[i];
	tau.yz = temp.x;
	tau.zz = temp.y;
	return tau;
}

//! Store tau tensor to split arrays
__device__
void storeTau(symtensor3 const& tau, uint i,
	double2 *__restrict__ tau0,
	double2 *__restrict__ tau1,
	double2 *__restrict__ tau2)
{
	tau0[i] = make_double2(tau.xx, tau.xy);
	tau1[i] = make_double2(tau.xz, tau.yy);
	tau2[i] = make_double2(tau.yz, tau.zz);
}


__device__ double& tensor3::operator()(int row, int col)
{
    // assert(col >= 0 && col < 4);
    // assert(row >= 0 && row < 4);

    return m_data[row][col];
}

 __device__ double  tensor3::operator()(int row, int col) const
{
    // assert(col >= 0 && col < 4);
    // assert(row >= 0 && row < 4);

    return m_data[row][col];
}

__device__ void tensor3::operator()() {
    // reset all elements of the tensor3 to 0.0
    for (int row{ 0 }; row < 4; ++row) {
        for (int col{ 0 }; col < 4; ++col) {
            m_data[row][col] = 0.0;
        }
    }
}

// __device__ tensor3::tensor3(){
	// data = new double[3];
	// for (int i=0;i<9;i++) data[i] = 0.;
// }



__device__ double& tensor3::operator[](const int &i){

	//return m_data[3*i];
}

__device__ tensor3 tensor3:: operator* (const tensor3 &b){
	tensor3 ret;
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
			for (int k=0;k<3;k++)
				ret.m_data[i][j]+=m_data[i][k]*b.m_data[k][j];
	
	return ret;
}

__device__ tensor3 tensor3:: operator+ (const tensor3 &b){
	tensor3 ret;
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
				ret.m_data[i][j]=m_data[i][j]+b(i,j);
	
	return ret;
}

__device__ tensor3&  tensor3::operator+= (const tensor3 &b){
	//tensor3 ret;
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
				m_data[i][j]=m_data[i][j]+b(i,j);
	// ret = this;
	return *this;
}

__device__ tensor3 tensor3:: operator- (const tensor3 &b){
	tensor3 ret;
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
				ret.m_data[i][j]=m_data[i][j]-b.m_data[i][j];
	
	return ret;
}

__device__ double3 tensor3::operator* (const double3 &v){
	double3 ret;
	ret.x = m_data[0][0]*v.x+m_data[0][1]*v.y+m_data[0][2]*v.z;	
	ret.y = m_data[1][0]*v.x+m_data[1][1]*v.y+m_data[1][2]*v.z;	
	ret.z = m_data[2][0]*v.x+m_data[2][1]*v.y+m_data[2][2]*v.z;	
	return ret;
}

__device__ tensor3 tensor3::operator*= (const double &v){
	tensor3 ret;
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++){
				ret.m_data[i][j]=m_data[i][j]*v;
				m_data[i][j] = ret.m_data[i][j];
		}
	return ret;
}

__device__ tensor3& tensor3::operator* (const double &f){
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
			m_data[i][j] *= f;	
	return *this;
}


__device__ tensor3 tensor3:: operator- (const double &f){
	tensor3 ret;
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
				ret(i,j)=m_data[i][j]-f;
	
	return ret;
}

__device__ tensor3 tensor3:: Trans (){
	tensor3 ret;
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
			ret.m_data[i][j] = m_data[j][i];
	
	return ret;
}

__device__ tensor3 operator* (const double &f, const tensor3 &b){
	tensor3 ret;
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
			ret.m_data[i][j] = b(i,j)*f;	
	return ret;
}

__device__ tensor3 tensor3::operator= (const double &f){
	tensor3 ret;
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
			m_data[i][j] = f;
	return ret;
}

// __device__ void tensor3::operator= (tensor3 &t){
		// for (int i=0;i<3;i++)
		// for (int j=0;j<3;j++)
			// m_data[i][j] = t.m_data[i][j];
// }


__device__ /*__forceinline__*/ tensor3 Identity(){
	tensor3 ret;
	ret(0,0) = ret(1,1) = ret(2,2) = 1.;
	//ret[1][1]=ret[2][2]=1.;
	
	return ret;
}

/*__spec*/
__device__
double3
dot(tensor3 const& T, double3 const& v)
{
	return make_double3(
			T(0,0)*v.x + T(0,1)*v.y +T(0,2)*v.z,
			T(1,0)*v.x + T(1,1)*v.y +T(1,2)*v.z,			
			T(2,0)*v.x + T(2,1)*v.y +T(2,2)*v.z
			);

}

__device__ tensor3 operator/ (const tensor3 &b, const double &f){

	tensor3 ret;
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
			ret(i,j) = b(i,j)/f;	
			
}

__device__ void tensor3::print(){
	printf("[%f %f %f],\n[%f %f %f],\n[%f %f %f].\n",
	m_data[0][0],m_data[0][1],m_data[0][2],
	m_data[1][0],m_data[1][1],m_data[1][2],
	m_data[2][0],m_data[2][1],m_data[2][2]);
}

#endif
