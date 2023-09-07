#ifndef _ELEMENTDATA_H_
#define _ELEMENTDATA_H_

struct ElementData {
  public:
  //// ELEMENT DATA, WITH GAUSS POINT
  double *pressure;
  double *str_rate;
  double *str_inc;
  
  double *dHxy_detJ; //(e,gp,dim,n)
  
};


////// NOT A FUNCTION OF STRUCT!!!!
////// APPENDING
void AllocateElements(ElementData *elem, const int &dim, const int &el_count);

#endif