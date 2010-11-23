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

  const Z i = threadIdx.x + (threadIdx.y ? 0 : RADIUS);
  const Z ii= threadIdx.x + (threadIdx.y ? nx - RADIUS : 0);
  const Z j = threadIdx.z;
  const Z k = blockIdx.x;
  const Z h = blockIdx.y;

  R *data  = f + h * ntotal + offset;
  R *ghost = data + ndata + offset + nx * (2 * RADIUS) * nz;

  ghost[(k * ny + j) * 2 * RADIUS + i] = data[(k * ny + j) * nx + ii];
}

__global__ void periodic_y(R *f, const Z nx, const Z ny, const Z nz)
{
  const Z offset = nx * ny * RADIUS;
  const Z ndata  = nx * ny * nz;
  const Z ntotal = ndata + (nx * ny + ny * nz + nz * nx) * (2 * RADIUS);

  const Z i = threadIdx.x;
  const Z j = threadIdx.y + (threadIdx.z ? 0 : RADIUS);
  const Z jj= threadIdx.y + (threadIdx.z ? ny - RADIUS : 0);
  const Z k = blockIdx.x;
  const Z h = blockIdx.y;

  R *data  = f + h * ntotal + offset;
  R *ghost = data + ndata + offset;

  ghost[(k * 2 * RADIUS + j) * nx + i] = data[(k * ny + jj) * nx + i];
}

__global__ void periodic_z(R *f, const Z nx, const Z ny, const Z nz)
{
  const Z offset = nx * ny * RADIUS;
  const Z ndata  = nx * ny * nz;
  const Z ntotal = ndata + (nx * ny + ny * nz + nz * nx) * (2 * RADIUS);

  const Z l = blockIdx.x * blockDim.x + threadIdx.x;
  const Z ll= l + offset;
  const Z h = blockIdx.y;

  if(l < offset) {
    R *lower = f + h * ntotal;
    R *upper = f + h * ntotal + ndata;
    lower[l ] = upper[l ];
    upper[ll] = lower[ll];
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

  bsz_x = make_uint3(RADIUS, 2, ny);
  gsz_x = make_uint3(nz, N_VAR, 1);

  bsz_y = make_uint3(nx, RADIUS, 2);
  gsz_y = make_uint3(nz, N_VAR, 1);

  bsz_z = make_uint3(dev.maxThreadsPerBlock / 4, 1, 1);
  gsz_z = make_uint3((nx * ny * RADIUS + bsz_z.x - 1) / bsz_z.x, N_VAR, 1);
}

void update_ghosts(R *f)
{
  periodic_x<<<gsz_x, bsz_x>>>(f, nx, ny, nz);
  periodic_y<<<gsz_y, bsz_y>>>(f, nx, ny, nz);
  periodic_z<<<gsz_z, bsz_z>>>(f, nx, ny, nz);
}
