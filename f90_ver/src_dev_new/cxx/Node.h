#ifndef _NODE_H_
#define _NODE_H_

#include <list>
#include "Math/Vec3D.h"

class Element;

class Node{
	
public:  
	Node(){}
	Node(const Vec3D &v){m_pos=v;}
  //const double & operator(int) const;
    
  Vec3D m_pos;	//Coordinates, from DynELA format
  
  ~Node(){}
};


#endif