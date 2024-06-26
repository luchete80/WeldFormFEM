#ifndef _DOUBLE3_C_H
#define _DOUBLE3_C_H


struct double3 {
  double x, y, z;
  double3(){x=y=z=0.0;}
  //double3(double x_, double y_, double z_){x=x_;y=y_;z=z_;}
};

inline double3 make_double3(double x_,double y_,double z_){double3 ret; ret.x=x_;
                                                                        ret.y=y_;
                                                                        ret.z=z_;
                                                          return ret;}

inline double3 operator+(double3 v1, double3 v2){return make_double3(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z);}
inline double3 operator-(double3 v1, double3 v2){return make_double3(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z);}
inline double3 operator*(double  s, double3 v  ){return make_double3(v.x*s,v.y*s,v.z*s);}

#endif