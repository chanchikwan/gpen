pro vis, i

  name = string(i, format='(i04)') + '.raw'
  print, 'loading: ' + name

  openr, lun, name, /get_lun
    nx = 0L & readu, lun, nx
    ny = 0L & readu, lun, ny
    nz = 0L & readu, lun, nz
    nv = 0L & readu, lun, nv
    f = fltarr(nx, ny, nz, nv)
    readu, lun, f
  close, lun & free_lun, lun

  x = findgen(nx) / nx
  y = findgen(ny) / ny
  z = findgen(nz) / nz

  ; hard-wire 4 variables for hydro
  rho = exp(f[*,*,*,0])
  ux  = f[*,*,*,1]
  uy  = f[*,*,*,2]
  uz  = f[*,*,*,3]

  print, min(rho), max(rho)
  
  window, 0, xSize=2*nx, ySize=2*ny
  tvscl, rho[*,*,nz/2-1], 0
  tvscl,  ux[*,*,nz/2-1], 1
  tvscl,  uy[*,*,nz/2-1], 2
  tvscl,  uz[*,*,nz/2-1], 3

end
