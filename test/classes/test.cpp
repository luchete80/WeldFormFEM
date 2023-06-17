class Obj{

};

extern "C" {
  void hello_();


  void readstep_(char* name, Obj** ptr){
    *ptr = new Obj(); //use name in the actual process
  }
  void pass_it_(Obj** ptr){
    hello_();
    delete *ptr; //some usage of the object
  }

}