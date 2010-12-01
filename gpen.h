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

#if defined(DOUBLE) || defined(OUBLE)
#  define TILE_X ( 8) /* x-thickness of the pencil bundles */
#else
#  define TILE_X (16) /* x-thickness of the pencil bundles */
#endif
#define TILE_Y ( 8) /* y-thickness of the pencil bundles */

void _error(const char *, const int, const char *, const char *);
#define error(message) _error(__FILE__, __LINE__, __func__, message)

R *initialize_modules(const R, const Z, const Z, const Z,
                               const R, const R, const R);

Z output_raw(Z, const R *);
Z output_vtk(Z, const R *);

void initial_condition(R *, Q (*)(double, double, double));
void update_ghosts(R *);

void rk_2n(R *, const R);
void pde(const R *, R *);

#endif
