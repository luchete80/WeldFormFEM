#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include "Node.h"

class Element {
	
protected:
  std::vector <Node *> m_node; //CHECK LIST PERFORMANCE
	
public:
  /*inline*/ virtual void computeShapeFunctions(){}; //INLINE

  
};

#endif