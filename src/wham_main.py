#!/usr/bin/env python

def Main(Iargs):
    return

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument('-rlow', default=1.5,type=float,
                         help='Low cutoff in angstroms')
    parser.add_argument('-rhigh', default=8,type=float,
                         help='High cutoff in angstroms')
    parser.add_argument('-nbin', default=100, type=int,
                         help='Histogram bins')
    parser.add_argument('-subfile', default="lif.distance", type=str,
                         help="File name for colvar in subdirectory")
    parser.add_argument('-subcol',default=1, type=int, 
                         help="Column of subfile for colvar")
    parser.add_argument('-skip', default=10000, type=int,
                         help="How many elements to skip")
    parser.add_argument('-deriv',default=0, type=int,
                         help="[0] Turn off derivative calculation (default) [1] Turn on derivative calculation")
    parser.add_argument('-enerfile', default="flucts.inp", type=str,
                         help="File name with energy information, default flucts.inp")
    Iargs = parser.parse_args()

    Main(Iargs)