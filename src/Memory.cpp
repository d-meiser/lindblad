#include <Memory.hpp>

#if defined(_WIN32)
#include <malloc.h>
#define LINDBLAD_ALLOCA _alloca
#else
#include <alloca.h>
#define LINDBLAD_ALLOCA alloca
#endif

void *alignedAlloca(size_t size, size_t alignment) {
  void *result = LINDBLAD_ALLOCA(size + alignment);
  int misalignment = ((size_t)result) % alignment;
  if (misalignment) {
    result = (void*)((char*)result + (alignment - misalignment));
  }
  return result;
}


