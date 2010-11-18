#include "gpen.h"

static __constant__ Z n;

static Z gsz, bsz;

static R *res;

__global__ void zero(R *f)
{
  const Z l = blockIdx.x * blockDim.x + threadIdx.x;

  if(l < n * N_VAR) f[l] = 0.0;
}

__global__ void kernel(R *f, R *g, const R dt_beta, const R alpha)
{
  const Z l = blockIdx.x * blockDim.x + threadIdx.x;

  if(l < n * N_VAR) {
    R F = f[l];
    R G = g[l];
    F += dt_beta * G;
    G *= alpha;
    f[l] = F;
    g[l] = G;
  }
}

void initialize_rk_2n(void *g, const Z nx, const Z ny, const Z nz)
{
  cudaError_t err;

  const Z n = nx * ny * nz;

  cudaDeviceProp dev;
  err = cudaGetDeviceProperties(&dev, 0);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  err = cudaMemcpyToSymbol("n", &n, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  bsz = dev.maxThreadsPerBlock;
  gsz = (n * N_VAR + bsz - 1) / bsz;

  res = (R *)g;
}

void rk_2n(R *f, const R dt)
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
