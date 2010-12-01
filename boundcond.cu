#include "gpen.h"

static Z Nx, Ny, Nz;

static uint3 Gsz_x, Bsz_x;
static uint3 Gsz_y, Bsz_y;
static uint3 Gsz_z, Bsz_z;

__global__ void periodic_x(R *f, const Z nx, const Z ny, const Z nz)
{
  const Z xghost = ny * nz * RADIUS;
  const Z hghost = (ny + nz) * nx * RADIUS + xghost;
  const Z ntotal = nx * ny * nz + 2 * hghost;

  const Z l = blockIdx.x * blockDim.y + threadIdx.y;

  R *mghost = f + blockIdx.y * ntotal;
  R *data   = mghost + hghost;
  R *pghost = mghost + ntotal - xghost;

  const Z i = threadIdx.x;

  if(l < ny * nz) {
    mghost += l * RADIUS;
    data   += l * nx;
    pghost += l * RADIUS;

    mghost[i] = data[i + nx - RADIUS];
    pghost[i] = data[i];
  }
}

__global__ void periodic_y(R *f, const Z nx, const Z ny, const Z nz)
{
  const Z ndata  = nx * ny * nz;
  const Z xghost = ny * nz * RADIUS;
  const Z yghost = nz * nx * RADIUS;
  const Z zghost = nx * ny * RADIUS;
  const Z ntotal = ndata + 2 * (xghost + yghost + zghost);

  const Z l = blockIdx.x * blockDim.x + threadIdx.x;

  R *mghost = f + blockIdx.y * ntotal + xghost;
  R *data   = mghost + yghost + zghost;
  R *pghost = data + ndata + zghost;

  const Z count  = nx * RADIUS;
  const Z stride = nx * ny;
  Z k;

  if(l < count) for(k = 0; k < nz; ++k) {
    mghost[l] = data[l + nx * (ny - RADIUS)];
    pghost[l] = data[l];

    mghost += count;
    data   += stride;
    pghost += count;
  }
}

__global__ void periodic_z(R *f, const Z nx, const Z ny, const Z nz)
{
  const Z ndata  = nx * ny * nz;
  const Z zghost = nx * ny * RADIUS;
  const Z hghost = (nx + ny) * nz * RADIUS + zghost;
  const Z ntotal = ndata + 2 * hghost;

  const Z l = blockIdx.x * blockDim.x + threadIdx.x;

  R *data   = f + blockIdx.y * ntotal + hghost;
  R *mghost = data - zghost;
  R *pghost = data + ndata;

  if(l < zghost) {
    mghost[l] = data[l + ndata - zghost];
    pghost[l] = data[l];
  }
}

void initialize_boundcond(const Z nx, const Z ny, const Z nz)
{
  cudaError_t err;

  cudaDeviceProp dev;
  err = cudaGetDeviceProperties(&dev, 0);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  Nx = nx;
  Ny = ny;
  Nz = nz;

  Bsz_x.x = RADIUS;
  Bsz_x.y = dev.maxThreadsPerBlock / RADIUS;
  Bsz_x.z = 1;
  Gsz_x.x = (ny * nz + Bsz_x.y - 1) / Bsz_x.y;
  Gsz_x.y = N_VAR;
  Gsz_x.z = 1;

  Bsz_y.x = dev.maxThreadsPerBlock / 4;
  Bsz_y.y = 1;
  Bsz_y.z = 1;
  Gsz_y.x = (nx * RADIUS + Bsz_y.x - 1) / Bsz_y.x;
  Gsz_y.y = N_VAR;
  Gsz_y.z = 1;

  Bsz_z.x = dev.maxThreadsPerBlock / 4;
  Bsz_z.y = 1;
  Bsz_z.z = 1;
  Gsz_z.x = (nx * ny * RADIUS + Bsz_z.x - 1) / Bsz_z.x;
  Gsz_z.y = N_VAR;
  Gsz_z.z = 1;
}

void update_ghosts(R *f)
{
  periodic_x<<<Gsz_x, Bsz_x>>>(f, Nx, Ny, Nz);
  periodic_y<<<Gsz_y, Bsz_y>>>(f, Nx, Ny, Nz);
  periodic_z<<<Gsz_z, Bsz_z>>>(f, Nx, Ny, Nz);
}
