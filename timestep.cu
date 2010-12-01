#include "gpen.h"

static Z Ndata, Hghost, Ntotal;

static uint3 Gsz, Bsz;

static R *res;

__global__ void zero(R *f, const Z ndata)
{
  const Z l = blockIdx.x * blockDim.x + threadIdx.x;

  if(l < ndata) f[blockIdx.y * ndata + l] = 0.0;
}

__global__ void kernel(R *f, R *g, const R dt_beta, const R alpha,
                       const Z ndata, const Z hghost, const Z ntotal)
{
  const Z l = blockIdx.x * blockDim.x + threadIdx.x;

  if(l < ndata) {
    const Z lf = blockIdx.y * ntotal + l + hghost;
    const Z lg = blockIdx.y * ndata  + l;

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

  cudaDeviceProp dev;
  err = cudaGetDeviceProperties(&dev, 0);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  Ndata  = nx * ny * nz;
  Hghost = (nx * ny + ny * nz + nz * nx) * RADIUS;
  Ntotal = Ndata + 2 * Hghost;

  Bsz.x = dev.maxThreadsPerBlock / 4; /* This increases performance...? */
  Bsz.y = 1;
  Bsz.z = 1;

  Gsz.x = (Ndata + Bsz.x - 1) / Bsz.x;
  Gsz.y = N_VAR;
  Gsz.z = 1;

  res = (R *)g;
}

void rk_2n(R *f, const R dt)
{
  const R alpha[] = {-5./9., -153./128., 0.     };
  const R beta [] = { 1./3.,   15./16.,  8./ 15.};

  static Z i = 0;

  if(i == 0) /* once rk_2n() is called, i == 3 */
    zero<<<Gsz, Bsz>>>(res, Ndata);

  for(i = 0; i < 3; ++i) {
    update_ghosts(f);
    pde(f, res);
    kernel<<<Gsz, Bsz>>>(f, res, dt * beta[i], alpha[i], Ndata, Hghost, Ntotal);
    cudaThreadSynchronize();
  }
}
