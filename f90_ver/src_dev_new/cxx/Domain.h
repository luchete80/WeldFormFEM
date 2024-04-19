#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "Math/Vec3D.h"

class Element;
class Node;

#include <vector>;

class Domain {
public:
  Domain();
  Domain(const int &dim):m_dim(dim){}
	void AddBoxLength(Vec3D const & V, Vec3D const & L, const double &r);
	void AllocateNodes(const int &nc);
	Element* 	getElement(const int &e) {return m_element[e];}
	Node* 		getNode		(const int &n) {return m_node[n];}
	
protected:
  std::vector <Element*>  m_element;
  std::vector <Node*>     m_node;  
	int 										m_dim;
  
	
																//IN CUDA BEING double3
	// std::vector <Vec3D>		m_u;		//VELOCITY, IN  ORDER TO BE SIMILAR TO GPU VERSION
	// std::vector <Vec3D>		m_v;		//VELOCITY, IN  ORDER TO BE SIMILAR TO GPU VERSION	
	// std::vector <Vec3D>		m_a;		//VELOCITY, IN  ORDER TO BE SIMILAR TO GPU VERSION	
	
	
  unsigned long elem_count;
  unsigned long node_count;

};


#endif