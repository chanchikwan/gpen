#include "gpen.h"
#include "deriv.h"
#include "hydro.h"

#define XGHOST(l, k) xghost[(l) * Ntotal + (k) * xstride]
#define YGHOST(l, k) yghost[(l) * Ntotal + (k) * ystride]
#define  CACHE(l, k)  cache[l][k][threadIdx.y][threadIdx.x]
#define   DATA(l, k)   data[(l) * Ntotal + (k) * Stride]
#define    RES(l, k)    res[(l) * Ndata  + (k) * Stride]

static __constant__ Z Nx, Ny, Nz, Stride, Ndata, Ntotal; /* sizes */
static __constant__ Z Xghost, Yghost, Hghost; /* offsets */
static __constant__ R Dx[RADIUS], Dy[RADIUS], Dz[RADIUS];

static uint3 Gsz, Bsz;

__global__ void rolling_cache(R *res, const R *f)
{
  /* Cache and tile[][] includes ghosts so we can compute derivatives */
  __shared__ R cache[N_VAR][1 + 2 * RADIUS][TILE_Y][TILE_X];
  __shared__ R tile[TILE_Y + 2 * RADIUS][TILE_X + 2 * RADIUS];

  /* For tile[][], shifted by RADIUS to take care the ghost cells */
  const Z i = threadIdx.x + RADIUS;
  const Z j = threadIdx.y + RADIUS;

  /* Store derivatives */
  R fx[N_VAR], fy[N_VAR], fz[N_VAR];

  /* Free running index along z direction */
  Z k;

  /* Thread-dependent variables --- very complicated!!! */
  const R *data = NULL, *xghost = NULL, *yghost = NULL;
  Z xstride, ystride, ghosti, ghostj;
  {
    /* Global indices only needed to compute thread-dependent variables */
    const Z I = blockIdx.x * TILE_X + threadIdx.x;
    const Z J = blockIdx.y * TILE_Y + threadIdx.y;

    /* This thread is in the domain, offset data and res to correct
       locations for input and output */
    if(I < Nx && J < Ny) {
      Z in = J * Nx + I;
      data = f + Xghost + Yghost + in;
      res += in;
    }

    /* This thread needs to load -x ghost */
    if(threadIdx.x < RADIUS)
    {
      ghosti = threadIdx.x;
      if(blockIdx.x) {
        xghost  = f + Hghost + J * Nx + I - RADIUS;
        xstride = Stride;
      } else {
        xghost  = f + J * RADIUS + threadIdx.x;
        xstride = Ny * RADIUS;
      }
    }
    /* This thread needs to load +x ghost */
    else if(i >= TILE_X)
    {
      ghosti = i + RADIUS;
      if(I + RADIUS < Nx) {
        xghost  = f + Hghost + J * Nx + I + RADIUS;
        xstride = Stride;
      } else {
        xghost  = f + Ntotal - Xghost + J * RADIUS + (I + RADIUS - Nx);
        xstride = Ny * RADIUS;
      }
      if(blockIdx.x == gridDim.x - 1) {
        xghost += (i - TILE_X) - (I + RADIUS - Nx);
        ghosti -= gridDim.x * TILE_X - Nx;
      }
    }

    /* This thread needs to load -y ghost */
    if(threadIdx.y < RADIUS)
    {
      ghostj = threadIdx.y;
      if(blockIdx.y) {
        yghost  = f + Hghost + (J - RADIUS) * Nx + I;
        ystride = Stride;
      } else {
        yghost  = f + Xghost + threadIdx.y * Nx + I;
        ystride = RADIUS * Nx;
      }
    }
    /* This thread needs to load +y ghost */
    else if(j >= TILE_Y)
    {
      ghostj = j + RADIUS;
      if(J + RADIUS < Ny) {
        yghost  = f + Hghost + (J + RADIUS) * Nx + I;
        ystride = Stride;
      } else {
        yghost  = f + Ntotal - Xghost - Yghost + (J + RADIUS - Ny) * Nx + I;
        ystride = RADIUS * Nx;
      }
      if(blockIdx.y == gridDim.y - 1) {
        yghost += ((j - TILE_Y) - (J + RADIUS - Ny)) * Nx;
        ghostj -= gridDim.y * TILE_Y - Ny;
      }
    }

    /* Unwind the stack for I and J since they are not needed anymore */
  }

  /*========================================================================*/
  /*                             ROLLING  CACHE                             */
  /*------------------------------------------------------------------------*/
  /* Start up: pre-load G-Pens */
  if(data) {
    Z l;
    #pragma unroll
    for(l = 0; l < N_VAR; ++l) {
      #pragma unroll
      for(k = 0; k < 2 * RADIUS; ++k)
        CACHE(l, k + 1) = DATA(l, k);
    }
    data += (2 * RADIUS) * Stride; /* shift by z-ghost */
  }

  /*------------------------------------------------------------------------*/
  /* Main z-loop: scan over z-slides; note that we place the k-loop
     outside the l-loop in order to save memory */
  for(k = 0; k < Nz; ++k) {
    Z l;

    #pragma unroll
    for(l = 0; l < N_VAR; ++l) {
      Z m;
      #pragma unroll
      /* Shift data */
      for(m = 0; m < 2 * RADIUS; ++m)
        CACHE(l, m) = CACHE(l, m + 1);
      /* Load the next slide */
      if(data)
        CACHE(l, m) = DATA(l, k);
    }

    for(l = 0; l < N_VAR; ++l) {
      /* Copy cache into tile */
      tile[j][i] = CACHE(l, RADIUS);

      /* Load x-ghost */
      if(xghost) tile[j][ghosti] = XGHOST(l, k);

      /* Load y-ghost */
      if(yghost) tile[ghostj][i] = YGHOST(l, k);

      __syncthreads();

      /* Compute derivatives */
      fx[l] = D_X;
      fy[l] = D_Y;
      fz[l] = D_Z;

      __syncthreads();
    }

    /*--------------------------------------------------------------------*/
    /* Compute res */
    if(data) {
      /* Continuity equation */
      RES(0, k) -= ADVECTION(0) + fx[1] + fy[2] + fz[3];
      /* x-momentum equation */
      RES(1, k) -= ADVECTION(1);
      /* y-momentum equation */
      RES(2, k) -= ADVECTION(2);
      /* z-momentum equation */
      RES(3, k) -= ADVECTION(3);
    }
    __syncthreads();
  }
  /*------------------------------------------------------------------------*/
  /*                             ROLLING  CACHE                             */
  /*========================================================================*/

}

void initialize_pde(const Z nx, const Z ny, const Z nz,
                    const R lx, const R ly, const R lz)
{
  cudaError_t err;

  const Z xghost = ny * nz * RADIUS;
  const Z yghost = nz * nx * RADIUS;
  const Z hghost = nx * ny * RADIUS + xghost + yghost;

  const Z stride = nx * ny;
  const Z ndata  = nx * ny * nz;
  const Z ntotal = ndata + 2 * hghost;

  const R coef[] = {45.0/60.0, -9.0/60.0, 1.0/60.0};
  R d[RADIUS];
  Z i;

  err = cudaMemcpyToSymbol("Nx", &nx, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));
  err = cudaMemcpyToSymbol("Ny", &ny, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));
  err = cudaMemcpyToSymbol("Nz", &nz, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  err = cudaMemcpyToSymbol("Stride", &stride, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));
  err = cudaMemcpyToSymbol("Ndata",  &ndata,  sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));
  err = cudaMemcpyToSymbol("Ntotal", &ntotal, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  err = cudaMemcpyToSymbol("Xghost", &xghost, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));
  err = cudaMemcpyToSymbol("Yghost", &yghost, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));
  err = cudaMemcpyToSymbol("Hghost", &hghost, sizeof(Z));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  for(i = 0; i < RADIUS; ++i) d[i] = coef[i] / (lx / nx);
  err = cudaMemcpyToSymbol("Dx", d, RADIUS * sizeof(R));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  for(i = 0; i < RADIUS; ++i) d[i] = coef[i] / (ly / ny);
  err = cudaMemcpyToSymbol("Dy", d, RADIUS * sizeof(R));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  for(i = 0; i < RADIUS; ++i) d[i] = coef[i] / (lz / nz);
  err = cudaMemcpyToSymbol("Dz", d, RADIUS * sizeof(R));
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  Bsz.x = TILE_X;
  Bsz.y = TILE_Y;
  Bsz.z = 1;

  Gsz.x = (nx + TILE_X - 1) / TILE_X;
  Gsz.y = (ny + TILE_Y - 1) / TILE_Y;
  Gsz.z = 1;
}

void pde(const R *f, R *res)
{
  /* Use the memcpy() convention --- the output is the first argument */
  rolling_cache<<<Gsz, Bsz>>>(res, f);
}
