#include "gpen.h"

#define HACK(x, y) if(i < 0) (x) = (y)
#define GPEN(l, k) gpen[l][k][threadIdx.y][threadIdx.x]

static __constant__ Z nx, ny, nz, stride, ndata, ntotal;

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

      #pragma unroll
      for(l = 0; l < N_VAR; ++l) {
        const R *data  = f + l * ntotal;
        const R *ghost = data + ndata;

        /* Load y-ghost zones */
        ghost += nx * ny * (2 * RADIUS);
        if(threadIdx.y < RADIUS) {
          tile[j - RADIUS][i] = (blockIdx.y == 0) ?
            ghost[(k * 2 * RADIUS + j - RADIUS) * nx + I] :
            data [in - RADIUS * stride - RADIUS * nx];
          tile[j + TILE_Y][i] = (blockIdx.y == gridDim.y - 1) ?
            ghost[(k * 2 * RADIUS + j         ) * nx + I] :
            data [in - RADIUS * stride + TILE_Y * nx];
        }

        /* Load x-ghost zones */
        ghost += nx * (2 * RADIUS) * nz;
        if(threadIdx.x < RADIUS) {
          tile[j][i - RADIUS] = (blockIdx.x == 0) ?
            ghost[(k * ny + J) * 2 * RADIUS + i - RADIUS] :
            data [in - RADIUS * stride - RADIUS];
          tile[j][i + TILE_X] = (blockIdx.x == gridDim.x - 1) ?
            ghost[(k * ny + J) * 2 * RADIUS + i         ] :
            data [in - RADIUS * stride + TILE_X];
        }

        /* Copy saved data to tile[][] in order to compute derivatives */
        tile[j][i] = GPEN(l, RADIUS);
        __syncthreads();

        /* HACK: never reach, trick the compiler to measure performance */
        HACK(res[out + l * ndata], tile[i + RADIUS][j + RADIUS]);
      }

      /* TODO: compute res */

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

  err = cudaMemcpyToSymbol("ndata",  &n, sizeof(Z));
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
