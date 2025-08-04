#include "flash_atten.h"

#include "mc_scverify.h"

CCS_MAIN(int arc, char** argv) {
  fptype query[SeqLen][_d];
  fptype key[SeqLen][_d];
  fptype value[SeqLen][_d];
  fptype output[SeqLen][_d];

  read_txt_2d<fptype, SeqLen, _d>(query, "./input/q.txt");
  read_txt_2d<fptype, SeqLen, _d>(key  , "./input/k.txt");
  read_txt_2d<fptype, SeqLen, _d>(value, "./input/v.txt");

  std::cout << "-----------------FLASH-----------------\n";

  for(int i = 0; i < _Tc*_Bc; i++)
  flashAtten2<_Tc, _Bc, _d>(query[i], key, value, output[i]);
  
  for(int i = 0; i < _Tc*_Bc; i++) {
    for(int j = 0; j < _d; j++) {
      printf("%.4f ", output[i][j].to_float());
    }
    printf("\n");
  }
  
  std::cout << "-----------------FLASH_WEIGHT-----------------\n";
  
  for(int i = 0; i < _Tc*_Bc; i++)
  weightAtten2<_Tc, _Bc, _d>(query[i], key, value, output[i]);
  
  for(int i = 0; i < _Tc*_Bc; i++) {
    for(int j = 0; j < _d; j++) {
      printf("%.4f ", output[i][j].to_float());
    }
    printf("\n");
  }

  CCS_RETURN(0);
}