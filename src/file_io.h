#ifndef __FILE_IO_HH__
#define __FILE_IO_HH__

#include <sstream>
#include <fstream>
#include <iostream>

/// @brief Read 2-d txt array
/// and return a pointer array
template<class _ty, int N, int M>
void
read_txt_2d(_ty arr[N][M], const char *filename) {
  float temp_arr[N][M];

  std::ifstream file(filename);

  for(int row = 0; row < N; ++row)
  {
    std::string line;
    getline(file, line);
    if ( !file.good() ) {
      std::cout << "[ERROR]: File not found or open!!\n";
      exit(-1);
    }
    
    std::stringstream iss(line);

    for (int col = 0; col < M; ++col)
    {
      std::string val;
      getline(iss, val, ',');
      if (!iss)
          break;
      std::stringstream convertor(val);
      convertor >> temp_arr[row][col];
      arr[row][col] = (_ty)temp_arr[row][col];
    }
  }

  std::cout << filename << ": " << N << ", " << M << "\n";
}

/// @brief Write 2-d array to .txt file
template<int N, int M>
void
write_txt_2d(const char *filename, float *arr) {
  std::ofstream file(filename);

  if (file.is_open()) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
        file << arr[i*M + j] << ",\n"[j==(N-1)];
      }
    }
  } else {
    std::cout << "[FILE NOT OPEN]: " << filename << std::endl;
  }
}

/// @brief Read 3-d txt array
/// and return a pointer array
template<class _ty, int N, int M, int K>
void
read_txt_3d(_ty arr[N][M*K], const char *filename) {
  float temp_arr[N][M*K];

  std::ifstream file(filename);

  for(int row = 0; row < N; ++row)
  {
    for (int m = 0; m < M; ++m) {
      std::string line;
      getline(file, line);
      if ( !file.good() ) {
        std::cout << "[ERROR]: File not found or open!!\n";
        exit(-1);
      }
      
      std::stringstream iss(line);

      for (int k = 0; k < K; ++k)
      {
        std::string val;
        getline(iss, val, ',');
        if (!iss)
        break;
        std::stringstream convertor(val);
        convertor >> temp_arr[row][m*K+k];
        arr[row][m*K+k] = (_ty)temp_arr[row][m*K+k];
      }
    }
  }

  std::cout << filename << ": " << N << ", " << M << ", " << K << "\n";
}
#endif