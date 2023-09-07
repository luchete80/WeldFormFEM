#ifndef _NODE_DATA_H_
#define _NODE_DATA_H_

struct NodeData {

  //// NODE DATA, SEPARATE?
  double3* x; //Vector is double
	double3* v;
	double3* a;
	double3* u;
  
};

void AllocateNodeData(NodeData *node, const int &nc);

#endif