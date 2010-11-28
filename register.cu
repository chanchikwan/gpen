#include "gpen.h"

extern void initialize_io               (void *, const Z, const Z, const Z);
extern void initialize_initial_condition(void *, const Z, const Z, const Z,
                                                 const R, const R, const R);

extern void initialize_rk_2n            (void *, const Z, const Z, const Z);
extern void initialize_pde              (        const Z, const Z, const Z);

static void *Func, *Res, *Host;

static void done(void)
{
  cudaFree(Func);
  cudaFree(Res);
  free(Host);
}

R *initialize_modules(const Z nx, const Z ny, const Z nz,
                      const R lx, const R ly, const R lz)
{
  cudaError_t err;

  const Z ndata  = nx * ny * nz;
  const Z hghost = (nx * ny + ny * nz + nz * nx) * RADIUS;
  const Z ntotal = ndata + 2 * hghost;

  err = cudaMalloc(&Func, sizeof(R) * ntotal * N_VAR);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  err = cudaMalloc(&Res,  sizeof(R) * ndata  * N_VAR);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  Host = malloc(sizeof(R) * ndata * N_VAR);
  if(!Host) error("fail to allocate host memory");

  atexit(done);

  initialize_io               (Host, nx, ny, nz);
  initialize_initial_condition(Host, nx, ny, nz,
                                     lx, ly, lz);
  initialize_rk_2n            (Res,  nx, ny, nz);
  initialize_pde              (      nx, ny, nz);

  return (R *)Func;
}
