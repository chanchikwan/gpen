#include "gpen.h"

#define HACK(x, y) if(i < 0) (x) = (y)

static __constant__ Z nx, ny, nz, stride, ntotal;

static uint3 gsz, bsz;

__global__ void rolling_cache(R *res, const R *f)
{
  /* For tile[][], shifted by RADIUS to take care the ghost cells */
  const Z i = threadIdx.x + RADIUS;
  const Z j = threadIdx.y + RADIUS;

  /* For global array f[], just like the mn_loop in the pencil code */
  const Z I = blockIdx.x * blockDim.x + threadIdx.x;
  const Z J = blockIdx.y * blockDim.y + threadIdx.y;

  /* Cached data */
  Q gpen[1 + 2 * RADIUS];

  /* Tile with x- and y-ghost so we can compute derivatives */
  __shared__ R tile[TILE_Y + 2 * RADIUS][TILE_X + 2 * RADIUS];

  if(I < nx && J < ny) {
    Z k, in, out;

    /* Pick the input  variable */
    out = in = J * nx + I;

    /* Start up: pre-load G-Pens */
    #pragma unroll
    for(k = 0; k < 2 * RADIUS; ++k) {
      gpen[k+1].lnrho = f[in + 0 * ntotal];
      gpen[k+1].ux    = f[in + 1 * ntotal];
      gpen[k+1].uy    = f[in + 2 * ntotal];
      gpen[k+1].uz    = f[in + 3 * ntotal];
      in += stride;
    }

    for(k = 0; k < nz; ++k) {
      Z l;

      /* Data shifting */
      #pragma unroll
      for(l = 0; l < 2 * RADIUS; ++l) gpen[l] = gpen[l+1];

      /* Load the next z-slide */
      gpen[2 * RADIUS].lnrho = f[in + 0 * ntotal];
      gpen[2 * RADIUS].ux    = f[in + 1 * ntotal];
      gpen[2 * RADIUS].uy    = f[in + 2 * ntotal];
      gpen[2 * RADIUS].uz    = f[in + 3 * ntotal];

      /* TODO: compute res */

      /* HACK: never reach, trick the compiler so we can measure performance */
      HACK(res[out], gpen[j].lnrho + gpen[j].ux + gpen[j].uy + gpen[j].uz);

      in  += stride;
      out += stride;
    }
  }
}

void initialize_pde(const Z nx, const Z ny, const Z nz)
{
  cudaError_t err;

  const Z s = nx * ny;
  const Z n = s * nz;
  const Z m = n + (nx * ny + ny * nz + nz * nx) * 2 * RADIUS;

  err = cudaMemcpyToSymbol("nx", &nx, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  err = cudaMemcpyToSymbol("ny", &ny, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  err = cudaMemcpyToSymbol("nz", &nz, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  err = cudaMemcpyToSymbol("stride", &s, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  err = cudaMemcpyToSymbol("ntotal", &m, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  bsz.x = TILE_X;
  bsz.y = TILE_Y;
  bsz.z = 1;

  gsz.x = (nx + TILE_X - 1) / TILE_X;
  gsz.y = (ny + TILE_Y - 1) / TILE_Y;
  gsz.z = 1;
}

void pde(const R *f, R *res)
{
  /* Use the memcpy() convention --- the output is the first argument */
  rolling_cache<<<gsz, bsz>>>(res, f);
}
