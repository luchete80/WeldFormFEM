#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include "Node.h"

#include "IntegrationPoint.h"

class Element {
	
protected:
  std::vector <Node *> 							m_node; //CHECK LIST PERFORMANCE
  std::vector <IntegrationPoint *> 	m_int_point; //CHECK LIST PERFORMANCE
	
public:
  /*inline*/ virtual void computeShapeFunctions(){}; //INLINE

  
};

#endif