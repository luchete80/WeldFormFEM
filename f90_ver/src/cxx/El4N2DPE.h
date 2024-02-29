#ifndef _EL4N2DPE_H_
#define _EL4N2DPE_H_

#include "Element.h"

class El4N2DPE:
public Element{

public:   
  
	El4N2DPE(std::vector <Node *>n) {
			for (int i=0;i<4;i++){
				m_node.push_back(n[i]);
			}
	}

	~El4N2DPE	(){}
 };


#endif