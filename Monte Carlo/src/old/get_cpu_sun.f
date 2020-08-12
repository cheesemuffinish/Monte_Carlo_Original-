      real function get_cpu()

c returns elapsed CPU time using ETIME system function (Sun, Linux)

      real time(2)

      get_cpu = etime(time)


      return
      end
	
