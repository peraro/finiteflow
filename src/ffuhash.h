#ifndef FF_UHASH_H
#define FF_UHASH_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    unsigned int val[5];
  } FFUHashVal;

  FFUHashVal ffUHash(const void * data, size_t n_bytes);

#ifdef __cplusplus
}
#endif

#endif // FF_UHASH_H
