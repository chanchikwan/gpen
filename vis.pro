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
  rho = f[*,*,*,0]
  ux  = f[*,*,*,1]
  uy  = f[*,*,*,2]
  uz  = f[*,*,*,3]

  !p.multi=[0,4,3,0,1]
  shade_surf, reform(rho[*,*,nz/2])
  shade_surf, reform(rho[*,ny/2,*])
  shade_surf, reform(rho[nx/2,*,*])
  shade_surf, reform( ux[*,*,nz/2])
  shade_surf, reform( ux[*,ny/2,*])
  shade_surf, reform( ux[nx/2,*,*])
  shade_surf, reform( uy[*,*,nz/2])
  shade_surf, reform( uy[*,ny/2,*])
  shade_surf, reform( uy[nx/2,*,*])
  shade_surf, reform( uz[*,*,nz/2])
  shade_surf, reform( uz[*,ny/2,*])
  shade_surf, reform( uz[nx/2,*,*])

end
