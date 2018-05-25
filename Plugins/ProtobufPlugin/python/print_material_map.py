#!/usr/bin/env python
from __future__ import print_function

import argparse
import sys
import struct
import os
import math
import math

from google.protobuf.internal.decoder import _DecodeVarint32

sys.path.append(os.path.join(os.path.dirname(__file__), "gen"))
import MaterialMap_pb2

def read_material_maps(f):
    while True:
        size_buf = f.read(2)
        if not size_buf: break
        size, pos = _DecodeVarint32(size_buf, 0)
        assert pos == 2
        payload = f.read(size)
        msg = MaterialMap_pb2.MaterialMap()
        msg.ParseFromString(payload)

        yield msg

def fmt(cols_, numfmt="{:6.4f}", width=9):
    cols = []
    for c in cols_:
        if type(c) == float:
            s = numfmt.format(float(c))
        else:
            s = str(c)
        s = s.rjust(width)
        cols.append(s)

    return cols

def main():
    p = argparse.ArgumentParser()
    p.add_argument("matmap")
    p.add_argument("--check-nan", action="store_true")

    args = p.parse_args()


    nmaps = 0
    nmaps_inactive_bins = 0
    has_nan_inf = False

    with open(args.matmap, "rb") as f:
        for mmap in read_material_maps(f):
            print("ID (vol, lay, app, sen):", mmap.vol_id, mmap.lay_id, mmap.app_id, mmap.sen_id)
            active_bins = filter(lambda bin: bin.valid, mmap.bins)
            nbins = mmap.rows * mmap.cols
            print(len(active_bins), "active out of", nbins, "(", mmap.rows, "x", mmap.cols, ")")
            print()
            if len(active_bins) != nbins:
                nmaps_inactive_bins += 1
            nmaps += 1

            for bin in mmap.bins:
                vals = [bin.X0, bin.L0, bin.thickness, bin.rho, bin.A, bin.Z]
                print(bin.valid, *fmt(vals))
                if any(map(math.isnan, vals)) or any(map(math.isinf, vals)):
                    print("Bin has NaN or INF values!!")
                    has_nan_inf = True

    print("Found", nmaps, "material maps")
    print(nmaps_inactive_bins, "have inactive material bins!")
    if has_nan_inf:
        print("Bins with NaN or INF value were found!")

    if args.check_nan and has_nan_inf:
        sys.exit(1)

if "__main__" == __name__:
    main()
