#include "spktyp.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"  
#include <cassert>
#include <iostream>

using namespace LGMC_NS;

Memory::Memory() {}

/* ----------------------------------------------------------------------
   safe malloc 
------------------------------------------------------------------------- */

void *Memory::smalloc(bigint nbytes, const char *name)
{
    if (nbytes == 0) return NULL;

    void *ptr = malloc(size_t(nbytes));
    if(ptr == NULL) {
        std::cout << name << '\n';
        assert(false);
    }               
    return ptr;
}

/* ----------------------------------------------------------------------
   safe realloc 
------------------------------------------------------------------------- */

void *Memory::srealloc(void *ptr, bigint nbytes, const char *name)
{
    if (nbytes == 0) {
        destroy(ptr);
        return NULL;
    }

    ptr = realloc(ptr, size_t(nbytes));
    if(ptr == NULL) {
        std::cout << name << '\n';
        assert(false);
    }
    
    return ptr;
}

/* ----------------------------------------------------------------------
   safe free 
------------------------------------------------------------------------- */

void Memory::sfree(void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}

/* ----------------------------------------------------------------------
   erroneous usage of templated create/grow functions
------------------------------------------------------------------------- */

void Memory::fail(const char *name)
{
    std::cout << name << '\n';
    assert(false);
}
