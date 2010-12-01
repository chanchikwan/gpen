pro vis, i

  name = string(i, format='(i04)') + '.raw'
  print, 'loading: ' + name

  openr, lun, name, /get_lun

    lx = 0d & readu, lun, lx
    ly = 0d & readu, lun, ly
    lz = 0d & readu, lun, lz

    nx = 0L & readu, lun, nx
    ny = 0L & readu, lun, ny
    nz = 0L & readu, lun, nz
    nv = 0L & readu, lun, nv

    print, lx, ly, lz
    print, nx, ny, nz, nv

    f  = dblarr(nx, ny, nz, nv)
    readu, lun, f

  close, lun & free_lun, lun

  x = lx * findgen(nx) / nx
  y = ly * findgen(ny) / ny
  z = lz * findgen(nz) / nz

  ; hard-wire 4 variables for hydro
  rho = exp(f[*,*,*,0])
  ux  = f[*,*,*,1]
  uy  = f[*,*,*,2]
  uz  = f[*,*,*,3]

  print, min(rho), max(rho)
  
  window, 0, xSize=2*nx, ySize=2*ny
  tvscl, rho[*,*,nz/2], 0
  tvscl,  ux[*,*,nz/2], 1
  tvscl,  uy[*,*,nz/2], 2
  tvscl,  uz[*,*,nz/2], 3

end
