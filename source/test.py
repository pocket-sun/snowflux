import ctypes
snowflux = ctypes.CDLL("./src/libsnowflux.so")
snowflux.init()
paras = (ctypes.c_double * 9)(2.5, 2.5, 2.5, 9.5, 12., 15.6, 5., 5., 5.)
dist = ctypes.c_double(10.)
