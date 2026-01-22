# Note: the results were obtained by running simulation, and the time values are pasted here.
# so that the script for plots could easily get them.
# Results are in ms

# # Barnes Hut, theta 0.8
# bh2Dresults082D = {
#     'n=1000000 steps=100': 96829,
#     'n=1000000 steps=10': 9790,
#     'n=1000000 steps=1': 937,
#     'n=100000 steps=100': 5935,
#     'n=100000 steps=10': 620,
#     'n=100000 steps=1': 65,
#     'n=10000 steps=10': 25,
#     'n=10000 steps=100': 275,
#     'n=10000 steps=1': 3,
#     'n=1000 steps=100': 17,
#     'n=1000 steps=10': 1,
#     'n=1000 steps=1': 0,
#     'n=100 steps=100': 1,
#     'n=100 steps=10': 0,
#     'n=100 steps=1': 0,
#     'n=10 steps=100': 0,
#     'n=10 steps=10': 0,
#     'n=10 steps=1': 0
# }

# Barnes Hut theta 0.5
bh2Dresults052D = {
    'n=1000000 steps=100': 97060,
    'n=1000000 steps=10': 9750,
    'n=1000000 steps=1': 933,
    'n=100000 steps=100': 5938,
    'n=100000 steps=10': 629,
    'n=100000 steps=1': 68,
    'n=10000 steps=100': 226,
    'n=10000 steps=10': 25,
    'n=10000 steps=1': 3,
    'n=1000 steps=100': 17,
    'n=1000 steps=10': 1,
    'n=1000 steps=1': 0,
    'n=100 steps=100': 1,
    'n=100 steps=10': 0,
    'n=100 steps=1': 0,
    'n=10 steps=100': 0,
    'n=10 steps=10': 0,
    'n=10 steps=1': 0
}


# Barnes Hut theta 0.5, OMP, 2D
# OMP 4 threads
bh2Dresults052DOMP = {
    'n=1000000 steps=100': 95936,
    'n=1000000 steps=10': 9803,
    'n=1000000 steps=1': 967,
    'n=100000 steps=100': 6217,
    'n=100000 steps=10': 686,
    'n=100000 steps=1': 67,
    'n=10000 steps=100': 217,
    'n=10000 steps=10': 23,
    'n=10000 steps=1': 3,
    'n=1000 steps=100': 17,
    'n=1000 steps=10': 2,
    'n=1000 steps=1': 0,
    'n=100 steps=100': 2,
    'n=100 steps=10': 0,
    'n=100 steps=1': 0,
    'n=10 steps=100': 0,
    'n=10 steps=10': 0,
    'n=10 steps=1': 0
}


# MPI 2D, -n 4
mpi2Dresults_n4 = {
    # 'n=1000000 steps=100': ,
    # 'n=1000000 steps=10': ,
    # 'n=1000000 steps=1': ,
    # 'n=100000 steps=100': ,
    'n=100000 steps=10': 327880,
    'n=100000 steps=1': 29229,
    'n=10000 steps=100': 27806,
    'n=10000 steps=10': 2784,
    'n=10000 steps=1': 279,
    'n=1000 steps=100': 317,
    'n=1000 steps=10': 45,
    'n=1000 steps=1': 4,
    'n=100 steps=100': 3,
    'n=100 steps=10': 0,
    'n=100 steps=1': 0,
    'n=10 steps=100': 0,
    'n=10 steps=10': 0,
    'n=10 steps=1': 0
}


# MPI Reduced 2D, -n 4
mpiReduced2Dresults_n4 = {
    # 'n=1000000 steps=100': ,
    # 'n=1000000 steps=10': ,
    # 'n=1000000 steps=1': ,
    # 'n=100000 steps=100': ,
    'n=100000 steps=10': 43176,
    'n=100000 steps=1': 3952,
    'n=10000 steps=100': 8746,
    'n=10000 steps=10': 372,
    'n=10000 steps=1': 56,
    'n=1000 steps=100': 39,
    'n=1000 steps=10': 4,
    'n=1000 steps=1': 0,
    'n=100 steps=100': 1,
    'n=100 steps=10': 0,
    'n=100 steps=1': 0,
    # 'n=10 steps=100': ,
    # 'n=10 steps=10': ,
    # 'n=10 steps=1': 
}


# Serial 2D
serialResults2D = {
    # 'n=1000000 steps=100': ,
    # 'n=1000000 steps=10': ,
    # 'n=1000000 steps=1': ,
    # 'n=100000 steps=100': ,
    'n=100000 steps=10': 1.05445e+06,
    'n=100000 steps=1': 105668,
    'n=10000 steps=100': 102610,
    'n=10000 steps=10': 10234,
    'n=10000 steps=1': 1027,
    'n=1000 steps=100': 1031,
    'n=1000 steps=10': 105,
    'n=1000 steps=1': 10,
    'n=100 steps=100': 13,
    'n=100 steps=10': 1,
    'n=100 steps=1': 0,
    'n=10 steps=100': 0,
    'n=10 steps=10': 0,
    'n=10 steps=1': 0
}

# SerialReduced 2D
serialReducedResults2D = {
    # 'n=1000000 steps=100': ,
    # 'n=1000000 steps=10': ,
    # 'n=1000000 steps=1': ,
    # 'n=100000 steps=100': ,
    'n=100000 steps=10': 555798,
    'n=100000 steps=1': 55523,
    'n=10000 steps=100': 54014,
    'n=10000 steps=10': 5395,
    'n=10000 steps=1': 538,
    'n=1000 steps=100': 539,
    'n=1000 steps=10': 55,
    'n=1000 steps=1': 5,
    'n=100 steps=100': 5,
    'n=100 steps=10': 0,
    'n=100 steps=1': 0,
    'n=10 steps=100': 0,
    'n=10 steps=10': 0,
    'n=10 steps=1': 0
}




# Serial 2D, OMP 4
serialResults2DOMP4 = {
    # 'n=1000000 steps=100': ,
    # 'n=1000000 steps=10': ,
    # 'n=1000000 steps=1': ,
    # 'n=100000 steps=100': ,
    'n=100000 steps=10': 313199,
    'n=100000 steps=1': 16552,
    # 'n=10000 steps=100': ,
    'n=10000 steps=10': 1546,
    'n=10000 steps=1': 143,
    # 'n=1000 steps=100': ,
    'n=1000 steps=10': 16,
    'n=1000 steps=1': 2,
    # 'n=100 steps=100': ,
    'n=100 steps=10': 0,
    'n=100 steps=1': 0,
    # 'n=10 steps=100': ,
    'n=10 steps=10': 0,
    'n=10 steps=1': 0
}



# SerialReduced 2D, OMP 4
serialReducedResults2DOMP4 = {
    # 'n=1000000 steps=100': ,
    # 'n=1000000 steps=10': ,
    # 'n=1000000 steps=1': ,
    # 'n=100000 steps=100': ,
    'n=100000 steps=10': 164959,
    'n=100000 steps=1': 31170,
    # 'n=10000 steps=100': ,
    'n=10000 steps=10': 2880,
    'n=10000 steps=1': 273,
    # 'n=1000 steps=100': ,
    'n=1000 steps=10': 29,
    'n=1000 steps=1': 3,
    # 'n=100 steps=100': ,
    'n=100 steps=10': 0,
    'n=100 steps=1': 0,
    # 'n=10 steps=100': ,
    'n=10 steps=10': 0,
    'n=10 steps=1': 0
}

# BarnesHut, 3D, theta 0.5
bh3Dresults05 = {
    # 'n=1000000 steps=100': ,
    'n=1000000 steps=10': 9300,
    'n=1000000 steps=1': 937,
    # 'n=100000 steps=100': ,
    'n=100000 steps=10': 625,
    'n=100000 steps=1': 69,
    # 'n=10000 steps=100': ,
    'n=10000 steps=10': 37,
    'n=10000 steps=1': 5,
    # 'n=1000 steps=100': ,
    'n=1000 steps=10': 2,
    'n=1000 steps=1': 0,
    # 'n=100 steps=100': ,
    'n=100 steps=10': 0,
    'n=100 steps=1': 0,
    # 'n=10 steps=100': ,
    'n=10 steps=10': 0,
    'n=10 steps=1': 0
}