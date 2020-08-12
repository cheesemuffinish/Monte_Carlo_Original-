      FUNCTION gammp(a,x)
      REAL a,gammp,x
CU    USES gcf,gser
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.d0)pause 'bad arguments in gammp'
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.d0-gammcf
      end if
      return
      END
