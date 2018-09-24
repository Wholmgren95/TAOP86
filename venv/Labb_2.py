import numpy as np
import time
import sys

filc = " ".join(sys.argv[1:]).split('.')[0] + '.npz'
npzfile = np.load(filc)
c = npzfile['c']
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
    in_col=np.argmin(c)
    rc_min = c[in_col]
    ut_list=[]

    for a_rows,b_rows in zip(A,b):
        ut_list.append(b_rows/a_rows[in_col])

    inkvar = ut_list.index(min(ut_list))
    # bix[inkvar] = in_col
    # A[inkvar] = A[inkvar]/A[inkvar][in_col] #hela raden delat på minsta
    #
    # for i in range(len(A[0])):
    #     temp = A[i][in_col]
    #     if i is not inkvar:
    #         A[i] = A[i] - temp*A[inkvar]
    #
    # for a_rows,b_rows in zip(A,b):
    #     temp = a_rows[in_col]
    #     a_rows = a_rows - temp*A[inkvar]




    # calc right-hand-sides and reduced costs
    # --------
    # calc most negative reduced cost, rc_min,
    # and index for entering variable, inkvar
    # --------

    if rc_min >= -1.0E-12:
        print('Ready')
        iter = -1

        # construct solution, x, and check it
        # --------

        diffx = np.linalg.norm(x - xcheat)
        diffz = z - zcheat
        print('xdiff: ' + repr(diffx))
        print('zdiff: ' + repr(diffz))
    else:
        # calc entering column, a
        # --------

        if max(a) <= 0:
            # unbounded solution
            print('Unbounded solution!')
            iter = -1
        else:
            # calc leaving var, utgvar
            # --------

            print(' Iter: ' + repr(iter) + ' z: ' + repr(z) + ' rc: ' + repr(rc_min) + ' ink: ' + repr(
                inkvar + 1) + ' utg: ' + repr(utgvar + 1))

            # make new partition
            # --------

elapsed = time.time() - t1
print('Elapsed time: ' + repr(elapsed))