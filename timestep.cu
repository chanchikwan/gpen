#include "gpen.h"

static __constant__ Z n, m;

static uint3 gsz, bsz;

static Z offset;

static R *res;

__global__ void zero(R *f)
{
  const Z l = blockIdx.x * blockDim.x + threadIdx.x;

  if(l < n) f[blockIdx.y * n + l] = 0.0;
}

__global__ void kernel(R *f, R *g, const R dt_beta, const R alpha)
{
  const Z l = blockIdx.x * blockDim.x + threadIdx.x;

  if(l < n) {
    const Z lf = blockIdx.y * m + l;
    const Z lg = blockIdx.y * n + l;

    R F = f[lf];
    R G = g[lg];
    F += dt_beta * G;
    G *= alpha;
    f[lf] = F;
    g[lg] = G;
  }
}

void initialize_rk_2n(void *g, const Z nx, const Z ny, const Z nz)
{
  cudaError_t err;

  const Z n = nx * ny * nz;
  const Z m = n + (nx * ny + ny * nz + nz * nx) * (2 * RADIUS);

  cudaDeviceProp dev;
  err = cudaGetDeviceProperties(&dev, 0);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  err = cudaMemcpyToSymbol("n", &n, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  err = cudaMemcpyToSymbol("m", &m, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  bsz.x = dev.maxThreadsPerBlock / 4; /* This increases performance...? */
  bsz.y = 1;
  bsz.z = 1;

  gsz.x = (n + bsz.x - 1) / bsz.x;
  gsz.y = N_VAR;
  gsz.z = 1;

  offset = nx * ny * RADIUS;

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
    update_ghosts(f);
    pde(f, res);
    kernel<<<gsz, bsz>>>(f + offset, res, dt * beta[i], alpha[i]);
    cudaThreadSynchronize();
  }
}
