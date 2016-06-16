import hackathon_params as hp
import sys

N = sys.argv[1]
steps = sys.argv[2]
dx = sys.argv[3]
sweeps = sys.argv[4]

hp.hackParams(N, steps, dx, sweeps)
