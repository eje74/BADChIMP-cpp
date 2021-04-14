#!/usr/bin/env python3

from numpy import fromfile, transpose, save, int8

directory="/cluster/home/janlv/BADChIMP-cpp/"
file2 = "Porer-70kv2.modif.raw"  # void = 1, solid = 0
geo1 = fromfile(file2, dtype=int8)
geo2 = transpose(geo1.reshape((1422, 459, 455)))
save('Porer-70kv2', geo2)

