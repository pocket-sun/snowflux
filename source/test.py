import ctypes as ctp
snowflux = ctp.CDLL("./src/libsnowflux.so")
paras = (ctp.c_double * 9)(2.5, 2.5, 2.5, 9.5, 12., 15.6, 5., 5., 5.)
snowflux.log_likelihood.restype=ctp.c_double
snowflux.log_likelihood(paras)
