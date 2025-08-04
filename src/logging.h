#ifndef __LOGGING_H__
#define __LOGGING_H__

template<class T, int N=1, int M=1>
void
print(
  const char *name,
  T          *buff
){
  std::cout << "\n|||||||||| " << name << " ||||||||||\n";
  for(int i = 0; i < N; i++) {
    std::cout << "------------------------\n";
    for(int j = 0; j < M; j++) {
      if (j < M-1) {
        printf("%.4f ", buff[i*M+j].to_float());
      } else {
        printf("%.4f\n", buff[i*M+j].to_float());
      }
    }
  }
  std::cout << "____________________________________\n";
}

#endif