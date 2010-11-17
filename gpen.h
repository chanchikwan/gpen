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

void _error(const char *, const int, const char *, const char *);
#define error(message) _error(__FILE__, __LINE__, __func__, message)

Q *initialize_modules(const Z, const Z, const Z);

Z output(Z, const Q *);

void initial_condition(Q *, Q (*)(R, R, R));

void rk_2n(Q *, const R);

#endif
