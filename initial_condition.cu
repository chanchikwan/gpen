#include "gpen.h"

static R *Host;
static Z Nx, Ny, Nz;
static R Lx, Ly, Lz;

void initialize_initial_condition(void *h, const Z nx, const Z ny, const Z nz,
                                           const R lx, const R ly, const R lz)
{
  Host = (R *)h;

  Nx = nx;
  Ny = ny;
  Nz = nz;

  Lx = lx;
  Ly = ly;
  Lz = lz;
}

void initial_condition(R *f, Q (*f0)(double, double, double))
{
  cudaError_t err;

  const Z ndata  = Nx * Ny * Nz;
  const Z hghost = (Nx * Ny + Ny * Nz + Nz * Nx) * RADIUS;
  const Z ntotal = ndata + 2 * hghost;

  R *lnrho = Host + 0 * ndata;
  R *ux    = Host + 1 * ndata;
  R *uy    = Host + 2 * ndata;
  R *uz    = Host + 3 * ndata;

  Z h, i, j, k;

  for(k = 0; k < Nz; ++k) {
    const double z = (double)Lz * k / Nz;
    for(j = 0; j < Ny; ++j) {
      const double y = (double)Ly * j / Ny;
      for(i = 0; i < Nx; ++i) {
        const double x = (double)Lx * i / Nx;
        const Z l = (k * Ny + j) * Nx + i;
        const Q f = f0(x, y, z);
        lnrho[l] = f.lnrho;
        ux   [l] = f.ux   ;
        uy   [l] = f.uy   ;
        uz   [l] = f.uz   ;
      }
    }
  }

  for(h = 0; h < N_VAR; ++h) {
    err = cudaMemcpy(f + hghost + h * ntotal, Host + h * ndata,
                     sizeof(R) * ndata, cudaMemcpyHostToDevice);
    if(cudaSuccess != err) error(cudaGetErrorString(err));
  }
}
