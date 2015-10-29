#ifndef __NUCLEUS_HH
#define __NUCLEUS_HH

#include <iostream>
#include <math.h>
#include "TObject.h"
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <string>

using namespace std;

class Nucleus : public TObject{
 public:
  Nucleus(char*);
  Nucleus(int, int, double, const char*);
  Nucleus(int Z, int N, char*);
  void SetZ(int);
  void SetN(int);
  void SetMassExcess(double);  
  void SetMass(double);  
  void SetMass();  
  void SetSymbol(const char*);  
  int GetZfromSymbol(char*);  
  int GetZ();
  int GetN();
  int GetA();
  double GetMassExcess();
  double GetMass();
  double GetRadius();
  const char* GetSymbol();
 private:
  int fZ;
  int fN;
  double fMass;
  double fMassExcess;
  string fSymbol;
  ClassDef(Nucleus,1);
};
#endif
