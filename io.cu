#include <stdio.h>
#include "gpen.h"

static R *Host;
static Z Nx, Ny, Nz;

void initialize_io(void *h, const Z nx, const Z ny, const Z nz)
{
  Host = (R *)h;

  Nx = nx;
  Ny = ny;
  Nz = nz;
}

Z output(Z i, const R *f)
{
  cudaError_t err;

  const Z ndata  = Nx * Ny * Nz;
  const Z hghost = (Nx * Ny + Ny * Nz + Nz * Nx) * RADIUS;
  const Z ntotal = ndata + 2 * hghost;

  char  name[256];
  FILE *file;

  Z h;

  for(h = 0; h < N_VAR; ++h) {
    err = cudaMemcpy(Host + h * ndata, f + hghost + h * ntotal,
                     sizeof(R) * ndata, cudaMemcpyDeviceToHost);
    if(cudaSuccess != err) error(cudaGetErrorString(err));
  }

  sprintf(name, "%04d.raw", i);
  file = fopen(name, "wb");
  fwrite(&Nx,  sizeof(Z), 1, file);
  fwrite(&Ny,  sizeof(Z), 1, file);
  fwrite(&Nz,  sizeof(Z), 1, file);
  fwrite(&h,   sizeof(Z), 1, file); /* after the for-loop, h == N_VAR */
  fwrite(Host, sizeof(R), ndata * N_VAR, file);
  fclose(file);

  return i;
}
