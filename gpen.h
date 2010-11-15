#ifndef GPEN_H
#define GPEN_H

typedef int Z;
#if defined(DOUBLE) || defined(OUBLE) /* so -DOUBLE works */
typedef double R;
#else
typedef float R;
#endif
typedef struct {R lnrho, ux, uy, uz;} Q;

#define RADIUS ( 3) /* Half of the derivative order == width of ghost zone */
#define TILE_X (16) /* x-thickness of the thick pencil */
#define TILE_Y (16) /* y-thickness of the thick pencil */

#endif
