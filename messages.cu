#include <stdlib.h>
#include <stdio.h>
#include "gpen.h"

void _error(const char *file, const int line,
            const char *func, const char *message)
{
  fprintf(stderr, "ERROR: %s (%d): %s: %s\n", file, line, func, message);
  abort();
}
