      FUNCTION gammq(a,x)
      REAL a,gammq,x
CU    USES gcf,gser
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.d0)pause 'bad arguments in gammq'
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammq=1.d0-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      end if
      return
      END
