#ifndef _ELEMENT_H_
#define _ELEMENT_H_

class Node;

class Element {

  std::vector <Node *> m_node; //CHECK LIST PERFORMANCE
	
public:
  /*inline*/ virtual void computeShapeFunctions(){}; //INLINE

  
};