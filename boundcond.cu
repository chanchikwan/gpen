#include "gpen.h"

static Z nx, ny, nz;

static uint3 gsz_x, bsz_x;
static uint3 gsz_y, bsz_y;
static uint3 gsz_z, bsz_z;

__global__ void periodic_x(R *f, const Z nx, const Z ny, const Z nz)
{
  const Z offset = nx * ny * RADIUS;
  const Z ndata  = nx * ny * nz;
  const Z ntotal = ndata + (nx * ny + ny * nz + nz * nx) * (2 * RADIUS);

  const Z l = blockIdx.x * blockDim.y + threadIdx.y;
  const Z i = threadIdx.x;

  R *data  = f + blockIdx.y * ntotal + offset;
  R *ghost = data + ndata + offset + nx * (2 * RADIUS) * nz;

  if(l < ny * nz) {
    data  += l * nx;
    ghost += l * (2 * RADIUS);

    ghost[i] = data[i - RADIUS + nx];
    ghost[i + RADIUS] = data[i];
  }
}

__global__ void periodic_y(R *f, const Z nx, const Z ny, const Z nz)
{
  const Z offset = nx * ny * RADIUS;
  const Z ndata  = nx * ny * nz;
  const Z ntotal = ndata + (nx * ny + ny * nz + nz * nx) * (2 * RADIUS);

  const Z l = blockIdx.x * blockDim.x + threadIdx.x;
  Z k;

  R *data  = f + blockIdx.y * ntotal + offset;
  R *ghost = data + ndata + offset;

  if(l < nx * RADIUS) for(k = 0; k < nz; ++k) {
    ghost[l] = data[l - nx * RADIUS + nx * ny];
    ghost[l + nx * RADIUS] = data[l];

    data  += nx * ny;
    ghost += nx * (2 * RADIUS);
  }
}

__global__ void periodic_z(R *f, const Z nx, const Z ny, const Z nz)
{
  const Z offset = nx * ny * RADIUS;
  const Z ndata  = nx * ny * nz;
  const Z ntotal = ndata + (nx * ny + ny * nz + nz * nx) * (2 * RADIUS);

  const Z l = blockIdx.x * blockDim.x + threadIdx.x;
  const Z L = l + offset;

  f += blockIdx.y * ntotal;

  if(l < offset) {
    f[l] = f[l + ndata];
    f[L + ndata] = f[L];
  }
}

void initialize_boundcond(const Z n, const Z m, const Z l)
{
  cudaError_t err;

  cudaDeviceProp dev;
  err = cudaGetDeviceProperties(&dev, 0);
  if(cudaSuccess != err) error(cudaGetErrorString(err));

  nx = n;
  ny = m;
  nz = l;

  /* No assumption */
  bsz_x.x = RADIUS;
  bsz_x.y = dev.maxThreadsPerBlock / RADIUS;
  bsz_x.z = 1;
  gsz_x.x = (ny * nz * RADIUS + bsz_x.x * bsz_x.y - 1) / (bsz_x.x * bsz_x.y);
  gsz_x.y = N_VAR;
  gsz_x.z = 1;

  /* No assumption but need for-loop along nz within the kernel */
  bsz_y.x = dev.maxThreadsPerBlock / 4;
  bsz_y.y = 1;
  bsz_y.z = 1;
  gsz_y.x = (nx * RADIUS + bsz_y.x - 1) / bsz_y.x;
  gsz_y.y = N_VAR;
  gsz_y.z = 1;

  /* No assumption */
  bsz_z.x = dev.maxThreadsPerBlock / 4;
  bsz_z.y = 1;
  bsz_z.z = 1;
  gsz_z.x = (nx * ny * RADIUS + bsz_z.x - 1) / bsz_z.x;
  gsz_z.y = N_VAR;
  gsz_z.z = 1;
}

void update_ghosts(R *f)
{
  periodic_x<<<gsz_x, bsz_x>>>(f, nx, ny, nz);
  periodic_y<<<gsz_y, bsz_y>>>(f, nx, ny, nz);
  periodic_z<<<gsz_z, bsz_z>>>(f, nx, ny, nz);
}
