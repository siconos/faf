from ctypes import cdll, c_float, c_longlong, byref
import numpy as np
try:
    with_papi = True
    papi=cdll.LoadLibrary('/usr/local/lib/libpapi.so')
except:
    try:
        with_papi = True
        papi=cdll.LoadLibrary('/usr/lib/x86_64-linux-gnu/libpapi.so.5.4.3')
    except:
        try:
            with_papi = True
            papi=cdll.LoadLibrary('/home/bremond/faf/install/lib/libpapi.so')
        except:
            with_papi = False


def init_flop():
    if with_papi:
        ireal_time = c_float()
        iproc_time = c_float()
        iflpops = c_longlong()
        imflops = c_float()
        papi.PAPI_flops(byref(ireal_time), byref(iproc_time), byref(iflpops),
                        byref(imflops))


def get_flop(real_time, proc_time, flpops, mflops):
    if with_papi:
        r = papi.PAPI_flops(byref(real_time), byref(proc_time), byref(flpops),
                            byref(mflops))
    else:
        real_time.value = np.nan
        proc_time.value = np.nan
        flpops.value = -1
        mflops.value = np.nan
