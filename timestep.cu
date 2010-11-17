#include "gpen.h"

static __constant__ Z nx, ny;

static uint3 gsz, bsz;

static Q *res;

__global__ void zero(Q *f)
{
  const Q q0 = {0, 0, 0, 0};

  const Z n = (blockIdx.y * ny +  blockIdx.x) * nx + threadIdx.x;

  f[n] = q0;
}

__global__ void kernel(Q *f, Q *g, const R dt_beta, const R alpha)
{
  const Z n = (blockIdx.y * ny +  blockIdx.x) * nx + threadIdx.x;

  Q F = f[n]; /* sizeof(Q) read */
  Q G = g[n]; /* sizeof(Q) read */

  F.lnrho += dt_beta * G.lnrho;
  F.ux    += dt_beta * G.ux   ;
  F.uy    += dt_beta * G.uy   ;
  F.uz    += dt_beta * G.uz   ;
  /* 8 floating-point operations */

  G.lnrho *= alpha;
  G.ux    *= alpha;
  G.uy    *= alpha;
  G.uz    *= alpha;
  /* 4 floating-point operations */

  f[n] = F; /* sizeof(Q) write */
  g[n] = G; /* sizeof(Q) write */
}

void initialize_rk_2n(void *g, const Z nx, const Z ny, const Z nz)
{
  cudaError_t err;

  err = cudaMemcpyToSymbol("ny", &ny, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  err = cudaMemcpyToSymbol("nx", &nx, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  bsz = make_uint3(nx, 1,  1);
  gsz = make_uint3(ny, nz, 1);

  res = (Q *)g;
}

void rk_2n(Q *f, const R dt)
{
  const R alpha[] = {-5./9., -153./128., 0.     };
  const R beta [] = { 1./3.,   15./16.,  8./ 15.};

  static Z i = 0;

  if(i == 0) /* once rk_2n() is called, i == 3 */
    zero<<<gsz, bsz>>>(res);

  for(i = 0; i < 3; ++i) {
    /* TODO: boundary condition */
    /* TODO: get res */
    kernel<<<gsz, bsz>>>(f, res, dt * beta[i], alpha[i]);
    cudaThreadSynchronize();
  }
}
