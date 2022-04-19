import sys
#sys.path.append('/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-05-20')
sys.path.append('/home/byu.local/sh747/underwater/scott-hollingsworth/codes/python-general-signal-processing/byuarglib/byuarglib')

from autospec import autospec
from ESAUdata import ESAUdata

ESAUdata

Gxx,f,OASPL = autospec(ch1, fs, ns, N, unitflag, pref)