import ctypes as ctp
resflux = ctp.CDLL("./libSNvC.so")
inputs = (ctp.c_double * 9)(2.5, 2.5, 2.5, 9.5, 12., 15.6, 5., 5., 5.)
outputs = (ctp.c_double * 20)(0.)
resflux.pygetSpec(inputs, outputs)
