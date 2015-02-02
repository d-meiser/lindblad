int main(int argn, char **argv) {
  int a[10];
#pragma omp simd
  for (int i = 0; i < 10; ++i) {
    a[i] += 100;
  }
}

