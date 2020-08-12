readcol,'kmh.par',lam,cext,csca,kap,g,pol,thet

cext(*)=0.2
csca(*)=0.2
kap(*)=0.2
g(*)=0.0
pol(*)=0.99999999
thet(*)=90.0

writecol,'rayleigh.par',lam,cext,csca,kap,g,pol,thet,FMT='(e12.6,6(1x,f12.8))'

end