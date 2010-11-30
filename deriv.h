#define D_X ( Dx[0] * (tile[j][i+1] - tile[j][i-1]) \
            + Dx[1] * (tile[j][i+2] - tile[j][i-2]) \
            + Dx[2] * (tile[j][i+3] - tile[j][i-3]) )

#define D_Y ( Dy[0] * (tile[j+1][i] - tile[j-1][i]) \
            + Dy[1] * (tile[j+2][i] - tile[j-2][i]) \
            + Dy[2] * (tile[j+3][i] - tile[j-3][i]) )

#define D_Z ( Dz[0] * (CACHE(l, 4) - CACHE(l, 2)) \
            + Dz[1] * (CACHE(l, 5) - CACHE(l, 1)) \
            + Dz[2] * (CACHE(l, 6) - CACHE(l, 0)) )

#define D_XX ( Dxx[0] * (tile[j][i  ]               ) \
             + Dxx[1] * (tile[j][i+1] + tile[j][i-1]) \
             + Dxx[2] * (tile[j][i+2] + tile[j][i-2]) \
             + Dxx[3] * (tile[j][i+3] + tile[j][i-3]) )

#define D_YY ( Dyy[0] * (tile[j  ][i]               ) \
             + Dyy[1] * (tile[j+1][i] + tile[j-1][i]) \
             + Dyy[2] * (tile[j+2][i] + tile[j-2][i]) \
             + Dyy[3] * (tile[j+3][i] + tile[j-3][i]) )

#define D_ZZ ( Dzz[0] * (CACHE(l, 3)              ) \
             + Dzz[1] * (CACHE(l, 4) + CACHE(l, 2)) \
             + Dzz[2] * (CACHE(l, 5) + CACHE(l, 1)) \
             + Dzz[3] * (CACHE(l, 6) + CACHE(l, 0)) )
