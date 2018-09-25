import numpy as np
import time
import sys
banan=5
filc = " ".join(sys.argv[1:]).split('.')[0] + '.npz'
npzfile = np.load(filc)
c = npzfile['c']*-1
b = npzfile['b']
A = npzfile['A']
bix = npzfile['bix'] #KOlumnindex för tillåten startbas
zcheat = npzfile['zcheat']
xcheat = npzfile['xcheat']

bix = bix - 1

t1 = time.time()

[m, n] = np.shape(A)
print ('Rows: ' + repr(m) + ' cols: ' + repr(n))

# Create nix
nix = np.setdiff1d(range(n), bix)

B = A[:, bix]
N = A[:, nix]
cB = c[bix]
cN = c[nix]

# Start iterations
iter = 0
while iter >= 0:
    iter += 1
    in_col=np.argmax(c)
    ut_list=[]

    for a_rows,b_rows in zip(A,b):
        ut_list.append(b_rows/a_rows[in_col])

    ut_row = ut_list.index(min(ut_list))
    bix[ut_row] = in_col
    print(A[ut_row[in_col]])
    A[ut_row] = A[ut_row]/A[ut_row[in_col]]
    for a_rows,b_rows in zip(A,b):
        temp = a_rows[in_col]
        for i in range(len(a_rows)):
            a_rows[i] = a_rows[i] - temp*A[ut_row[i]]




    # calc right-hand-sides and reduced costs
    # --------
    # calc most negative reduced cost, rc_min,
    # and index for entering variable, inkvar
    # --------

    # if rc_min >= -1.0E-12:
    #     print('Ready')
    #     iter = -1
    #
    #     # construct solution, x, and check it
    #     # --------
    #
    #     diffx = np.linalg.norm(x - xcheat)
    #     diffz = z - zcheat
    #     print('xdiff: ' + repr(diffx))
    #     print('zdiff: ' + repr(diffz))
    # else:
    #     # calc entering column, a
    #     # --------
    #
    #     if max(a) <= 0:
    #         # unbounded solution
    #         print('Unbounded solution!')
    #         iter = -1
    #     else:
    #         # calc leaving var, utgvar
    #         # --------
    #
    #         print(' Iter: ' + repr(iter) + ' z: ' + repr(z) + ' rc: ' + repr(rc_min) + ' ink: ' + repr(
    #             inkvar + 1) + ' utg: ' + repr(utgvar + 1))
    #
    #         # make new partition
    #         # --------

#elapsed = time.time() - t1
#print('Elapsed time: ' + repr(elapsed))