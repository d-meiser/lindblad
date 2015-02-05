#ifndef MEMORY_HPP
#define MEMORY_HPP

#include <cstddef>

#if defined(_WIN32)
#include <malloc.h>
#define LINDBLAD_ALLOCA _alloca
#else
#include <alloca.h>
#define LINDBLAD_ALLOCA alloca
#endif

#define LINDBLAD_MAX_STACK_ALLOCATION_SIZE (1<<16)

#define LINDBLAD_ALIGNED_ALLOCA(size, alignment)                               \
  (void *)(                                                                    \
      (alignment) *                                                            \
      (((long long)(LINDBLAD_ALLOCA((size) + (alignment))) + (alignment)) /    \
       (alignment)));

#endif

