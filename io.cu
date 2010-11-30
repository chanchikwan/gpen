#include <stdio.h>
#include "gpen.h"

static R *Host;
static Z Nx, Ny, Nz;
static R Lx, Ly, Lz;

#define HEADER  "# vtk DataFile "
#define RESERVE HEADER "\nwriting...\n"
#define STAMP   HEADER "Version 2.0\n"

static void bswap(void *dat, size_t len, size_t cnt)
{
  const int h = len / 2;
  const int l = len - 1;

  size_t i;
  for(i = 0; i < cnt; ++i) {
    char *d = (char *)dat + i * len;

    int j;
    for(j = 0; j < h; ++j) {
      char temp = d[j];
      d[j] = d[l-j];
      d[l-j] = temp;
    }
  }
}

static int stamp(FILE *file, float time)
{
  long int pos = ftell(file);
  if(pos) {
    /* If not at the beginning of a file, go back and place stamp */
    fseek(file, 0L, SEEK_SET);
    fputs(STAMP, file);
    fseek(file, pos, SEEK_SET);
    return 0;
  } else {
    /* If at the beginning of a file, reserve space for the stamp */
    fputs(RESERVE, file);
    fprintf(file, "time = %lf\n", time);
    return 1;
  }
}

void initialize_io(void *h, const Z nx, const Z ny, const Z nz,
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

Z output_raw(Z i, const R *f)
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
  fwrite(&Lx,  sizeof(R), 1, file);
  fwrite(&Ly,  sizeof(R), 1, file);
  fwrite(&Lz,  sizeof(R), 1, file);
  fwrite(&Nx,  sizeof(Z), 1, file);
  fwrite(&Ny,  sizeof(Z), 1, file);
  fwrite(&Nz,  sizeof(Z), 1, file);
  fwrite(&h,   sizeof(Z), 1, file); /* after the for-loop, h == N_VAR */
  fwrite(Host, sizeof(R), ndata * N_VAR, file);
  fclose(file);

  return i;
}

Z output_vtk(Z i, const R *f)
{
  cudaError_t err;

  const Z ndata  = Nx * Ny * Nz;
  const Z hghost = (Nx * Ny + Ny * Nz + Nz * Nx) * RADIUS;
  const Z ntotal = ndata + 2 * hghost;

  char  name[256];
  FILE *file;

  Z h, j, k;

  const R *den= Host + 0 * ndata;
  const R *ux = Host + 1 * ndata;
  const R *uy = Host + 2 * ndata;
  const R *uz = Host + 3 * ndata;

  float *buff = (float *)malloc(sizeof(float) * 3 * Nx * Ny * Nz);

  for(h = 0; h < N_VAR; ++h) {
    err = cudaMemcpy(Host + h * ndata, f + hghost + h * ntotal,
                     sizeof(R) * ndata, cudaMemcpyDeviceToHost);
    if(cudaSuccess != err) error(cudaGetErrorString(err));
  }

  sprintf(name, "gpen-%04d.vtk", i);
  file = fopen(name, "wb");

  { /* Header */
    stamp(file, 0.0);
    fprintf(file, "BINARY\n");
  }

  { /* Grid */
    fprintf(file, "DATASET STRUCTURED_POINTS\n");
    fprintf(file, "DIMENSIONS %d %d %d\n", Nx, Ny, Nz);
    fprintf(file, "ORIGIN %e %e %e \n",    0., 0., 0.);
    fprintf(file, "SPACING %e %e %e \n",   1., 1., 1.);
    fprintf(file, "POINT_DATA %d\n",       Nx* Ny* Nz);
  }

  fprintf(file, "\nSCALARS d float\n");
  fprintf(file, "LOOKUP_TABLE default\n");
  for(i = 0; i < Nx; ++i)
    for(j = 0; j < Ny; ++j)
      for(k = 0; k < Nz; ++k)
        buff[(k * Ny + j) * Nx + i] = den[(i * Ny + j) * Nz + k];

#ifndef WORD_BIGENDIAN
  bswap(buff, sizeof(float), Nx * Ny * Nz);
#endif
  fwrite(buff, sizeof(float), Nx * Ny * Nz, file);

  fprintf(file, "\nSCALARS K float\n");
  fprintf(file, "LOOKUP_TABLE default\n");
  for(i = 0; i < Nx; ++i)
    for(j = 0; j < Ny; ++j)
      for(k = 0; k < Nz; ++k) {
        buff[(k * Ny + j) * Nx + i] =
          ux[(i * Ny + j) * Nz + k] * ux[(i * Ny + j) * Nz + k] +
          uy[(i * Ny + j) * Nz + k] * uy[(i * Ny + j) * Nz + k] +
          uz[(i * Ny + j) * Nz + k] * uz[(i * Ny + j) * Nz + k];
      }

#ifndef WORD_BIGENDIAN
  bswap(buff, sizeof(float), Nx * Ny * Nz);
#endif
  fwrite(buff, sizeof(float), Nx * Ny * Nz, file);

  fprintf(file, "\nVECTORS u float\n");
  for(i = 0; i < Nx; ++i)
    for(j = 0; j < Ny; ++j)
      for(k = 0; k < Nz; ++k) {
        buff[3 * ((k * Ny + j) * Nx + i) + 0] = ux[(i * Ny + j) * Nz + k];
        buff[3 * ((k * Ny + j) * Nx + i) + 1] = uy[(i * Ny + j) * Nz + k];
        buff[3 * ((k * Ny + j) * Nx + i) + 2] = uz[(i * Ny + j) * Nz + k];
      }

#ifndef WORD_BIGENDIAN
  bswap(buff, sizeof(float), 3 * Nx * Ny * Nz);
#endif
  fwrite(buff, sizeof(float), 3 * Nx * Ny * Nz, file);

  free(buff);

  return stamp(file, 0.0);
}
