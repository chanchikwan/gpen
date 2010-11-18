pro vis, i

  nv = 4

  name = string(i, format='(i04)') + '.raw'
  print, 'loading: ' + name

  openr, lun, name, /get_lun
    nx = 0L & readu, lun, nx
    ny = 0L & readu, lun, ny
    nz = 0L & readu, lun, nz
    f = fltarr(nx, ny, nz, nv)
    readu, lun, f
  close, lun & free_lun, lun

  x = findgen(nx) / nx
  y = findgen(ny) / ny
  z = findgen(nz) / nz
  lnrho = f[*,*,*,0]
  
  window, 0, xSize=nx, ySize=ny
  tvscl, exp(lnrho[*,*,nz/2-1])

end
