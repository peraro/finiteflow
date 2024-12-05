#include "ffuhash.h"
#include "WjCryptLib_Sha1.h"

_Static_assert(sizeof(SHA1_HASH) == sizeof(FFUHashVal),
               "FFUHashVal and SHA1_HASH have incompatible sizes");

FFUHashVal ffUHash(const void * data, size_t n_bytes)
{
  FFUHashVal ret = {0};
  Sha1Calculate(data, n_bytes, (SHA1_HASH*)(&ret));
  return ret;
}
