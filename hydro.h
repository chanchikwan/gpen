#define ADVECTION(I) ( CACHE(1, RADIUS) * fx[I] \
                     + CACHE(2, RADIUS) * fy[I] \
                     + CACHE(3, RADIUS) * fz[I] )
