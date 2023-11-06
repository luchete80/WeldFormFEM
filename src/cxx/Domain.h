#ifndef _DOMAIN_H_
#define _DOMAIN_H_

class Element;
class Node;
#include <vector>;

class Domain {
public:
  Domain();
	void AddBoxLength();

protected:
  std::vector <Element*>  m_element;
  std::vector <Node*>     m_node;  
	int 										m_dim;
	
  unsigned long elem_count;
  unsigned long node_count;

public:
	Element* 	getElement(const int &e) {return m_element[e];}
	Node* 		getNode		(const int &n) {return m_node[n];}
};


#endif