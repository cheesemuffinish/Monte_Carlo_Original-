# Call with
#
# $ python make_fine.py input output

import sys
import numpy as np
from hyperion.util.interpolate import interp1d_fast_loglog, \
                                      interp1d_fast_loglin

if len(sys.argv) != 3:
    print "Usage: python make_fine.py input output"
    sys.exit(0)

kmh = np.loadtxt(sys.argv[1])

# Use following to have uniform sampling of wavelengths
# wav = np.logspace(np.log10(kmh[0, 0]), np.log10(kmh[-1, 0]), 1000)

# Use following to refine only between two last wavelengths
wav = np.hstack([kmh[:-2, 0], np.logspace(np.log10(kmh[-2, 0]),
                                          np.log10(kmh[-1, 0]), 50)])

cext = interp1d_fast_loglog(kmh[:, 0], kmh[:, 1], wav)
csca = interp1d_fast_loglog(kmh[:, 0], kmh[:, 2], wav)
chi = interp1d_fast_loglog(kmh[:, 0], kmh[:, 3], wav)
c5 = interp1d_fast_loglin(kmh[:, 0], kmh[:, 4], wav)
c6 = interp1d_fast_loglin(kmh[:, 0], kmh[:, 5], wav)
c7 = interp1d_fast_loglin(kmh[:, 0], kmh[:, 6], wav)

np.savetxt(sys.argv[2], zip(wav, cext, csca, chi, c5, c6, c7), fmt='%11.4e')
