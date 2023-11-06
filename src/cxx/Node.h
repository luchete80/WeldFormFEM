#ifndef _NODE_H_
#define _NODE_H_

#incldue <list>
#include "Math/Vec3D.h"

class Element;

class Node{

public:  
  const double & operator(int) const;
    
  Vec3D m_pos;	//Coordinates, from DynELA format
  
  
};


#endif