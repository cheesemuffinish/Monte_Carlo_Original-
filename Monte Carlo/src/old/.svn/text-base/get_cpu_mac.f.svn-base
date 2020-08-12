      real function get_cpu()

c returns elapsed CPU time using ETIME system function (Sun, Linux, Mac)
c  need the extra _ for the mac to compile without the -N15 option 
c  cfitsio lib won't work with -N15 option.  what a frikkin pain
c  this compiler is.

      real time(2)

      get_cpu = etime_(time)


      return
      end
	
