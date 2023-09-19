#ifndef _MESH_H_
#define _MESH_H_

class Element;
class Node;

class Mesh {
public:
  Mesh();

protected:
  std::vector <Element*>  m_element;
  std::vector <Node*>     m_node;  
  
  unsigned long elem_count;
  unsigned long node_count;
};


#endif