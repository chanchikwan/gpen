#include "gpen.h"

#define HACK(x, y) if(i < 0) (x) = (y)
#define GPEN(l, k) gpen[l][k][threadIdx.y][threadIdx.x]

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

  /* Cached data, tile[][] includes ghosts so we can compute derivatives */
  __shared__ R gpen[N_VAR][1 + 2 * RADIUS][TILE_Y][TILE_X];
  __shared__ R tile[TILE_Y + 2 * RADIUS][TILE_X + 2 * RADIUS];

  if(I < nx && J < ny) {
    Z k, in, out;

    /* Pick the input  variable */
    out = in = J * nx + I;

    /* Start up: pre-load G-Pens */
    #pragma unroll
    for(k = 0; k < 2 * RADIUS; ++k) {
      Z l;
      #pragma unroll
      for(l = 0; l < N_VAR; ++l)
        GPEN(l, k + 1) = f[in + l * ntotal];
      in += stride;
    }

    for(k = 0; k < nz; ++k) {
      Z l;

      /* Data shifting and load the next z-slide */
      #pragma unroll
      for(l = 0; l < N_VAR; ++l) {
        Z m;
        #pragma unroll
        for(m = 0; m < 2 * RADIUS; ++m)
          GPEN(l, m) = GPEN(l, m + 1);
        GPEN(l, 2 * RADIUS) = f[in + l * ntotal];
      }

      /* TODO: compute res */

      /* HACK: never reach, trick the compiler so we can measure performance */
      HACK(res[out], GPEN(0, j) + GPEN(1, j) + GPEN(2, j) + GPEN(3, j));

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
