#ifndef GPEN_H
#define GPEN_H

typedef int Z;
#if defined(DOUBLE) || defined(OUBLE) /* so -DOUBLE works */
typedef double R;
#else
typedef float R;
#endif
typedef struct {R lnrho, ux, uy, uz;} Q;

#define N_VAR  ( 4) /* Number of variables */
#define RADIUS ( 3) /* Half of the derivative order == width of ghost zone */
#define TILE_X (16) /* x-thickness of the pencil bundles */
#define TILE_Y ( 8) /* y-thickness of the pencil bundles */

void _error(const char *, const int, const char *, const char *);
#define error(message) _error(__FILE__, __LINE__, __func__, message)

R *initialize_modules(const Z, const Z, const Z);

Z output(Z, const R *);

void initial_condition(R *, Q (*)(R, R, R));
void update_ghosts(R *);

void rk_2n(R *, const R);
void pde(const R *, R *);

#endif
