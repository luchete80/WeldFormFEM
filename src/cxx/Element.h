#ifndef _ELEMENT_H_
#define _ELEMENT_H_

class Node;

class Element {
public:
  /*inline*/ virtual void computeShapeFunctions(){}; //INLINE
  
  std::vector <Node *> m_node; //CHECK LIST PERFORMANCE
  
};