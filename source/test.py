import ctypes as ctp
snowflux = ctp.CDLL("./src/libsnowflux.so")
paras = (ctp.c_double * 9)(2.5, 2.5, 2.5, 9.5, 12., 15.6, 5., 5., 5.)
#dist = ctp.c_double(10.)
#res = (ctp.c_double * 14)(0.)
#res_sz = (ctp.c_ulong)(9)
