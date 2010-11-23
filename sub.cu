#include <stdlib.h>
#include "gpen.h"

Z atox(const char *arg)
{
  Z n = atoi(arg);
  return ((n + TILE_X - 1) / TILE_X) * TILE_X;
}

Z atoy(const char *arg)
{
  Z n = atoi(arg);
  return ((n + TILE_Y - 1) / TILE_Y) * TILE_Y;
}
