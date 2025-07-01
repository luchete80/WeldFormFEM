#ifndef _FUNCTION_H_
#define _FUNCTION_H_


struct function {
	int id;
	std::vector <double> time; ///set??
	//std::vector <T> value;
  //T getValAtTime(const double &t){
  std::vector <double> value;
  double getValAtTime(const double &t){
    double ret;
    //assumed ordered
    int i=time.size()-1;
    bool end = false;
    while (!end){
      i--;
      if (i==0 || t > time[i]) end =true;
    }
    ret =  value[i]+ (value[i+1] - value[i])/(time[i+1] - time[i])*(t-time[i]);
    return ret;
    
  }
	//std::map;
};

#endif
