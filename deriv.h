#define D_X (  (45.f/60) * (tile[j][i+1] - tile[j][i-1]) \
             - ( 9.f/60) * (tile[j][i+2] - tile[j][i-2]) \
             + ( 1.f/60) * (tile[j][i+3] - tile[j][i-3]))

#define D_Y (  (45.f/60) * (tile[j+1][i] - tile[j-1][i]) \
             - ( 9.f/60) * (tile[j+2][i] - tile[j-2][i]) \
             + ( 1.f/60) * (tile[j+3][i] - tile[j-3][i]))

#define D_Z (  (45.f/60) * (GPEN(l, 4) - GPEN(l, 2))     \
             - ( 9.f/60) * (GPEN(l, 5) - GPEN(l, 1))     \
             + ( 1.f/60) * (GPEN(l, 6) - GPEN(l, 0)))
