#ifndef SYM_TEST_MUTEX_H
#define SYM_TEST_MUTEX_H

#include <map>
#include <vector>
#include <string>

#include "cuddObj.hh"


class SymManager;

class GSTMutex{

  BDD notMutexBDD;
 public:
 GSTMutex(){}
  
  void check_mutexes(SymManager & manager);
};


extern GSTMutex gst_mutex;

#endif
