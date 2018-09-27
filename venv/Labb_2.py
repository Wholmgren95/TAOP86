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
xB = np.transpose(bix)
xN = np.transpose(nix)
x = 2 #xn + xb
# Start iterations
iter = 0
while iter >= 0:
    iter += 1  
    # bix[utgvar] = inkvar
    # A[utgvar] = A[utgvar]/A[utgvar][inkvar] #hela raden delat på minsta
    #
    # for i in range(len(A[0])):
    #     temp = A[i][inkvar]
    #     if i is not utgvar:
    #         A[i] = A[i] - temp*A[utgvar]
    #
    # for a_rows,b_rows in zip(A,b):
    #     temp = a_rows[inkvar]
    #     a_rows = a_rows - temp*A[utgvar]
    
    # calc right-hand-sides and reduced costs
    # --------
    

    B_inv = np.linalg.inv(B)
    xB = np.dot(B_inv, b)
    xN = 0
    y = np.transpose(np.dot(np.transpose(cB),B_inv))
    z = np.dot(np.transpose(b),y)
    c_hattN = cN - (np.dot(np.transpose(N),y)) #np.transpose(N)

    # calc most negative reduced cost, rc_min,
    # and index for entering variable, inkvar
    # --------

    inkvar = np.argmin(c)
    rc_min = c[inkvar]

    
    for e in c_hattN:
        if e <=0:
            break
    if rc_min >= -1.0E-12:
        print('Ready')
        iter = -1

        bix[utgvar] = inkvar
            # construct solution, x, and check it
            # --------
        diffx = np.linalg.norm(x - xcheat)
        diffz = z - zcheat
        print('xdiff: ' + repr(diffx))
        print('zdiff: ' + repr(diffz))
    else:
        # calc entering column, a
        # --------
        # a = inkommande kolumn hatt
        # A[:,inkvar] är att ta ut kolumn ur A på plats inkvar
        a = np.dot(B_inv, A[:, inkvar])

        if max(a) <= 0:
            # unbounded solution
            print('Unbounded solution!')
            iter = -1
        else:
            # calc leaving var, utgvar
            # --------
            ut_list = []

            for a_rows, b_rows in zip(A, b):
                ut_list.append(b_rows / a_rows[inkvar])

            utgvar = ut_list.index(min(ut_list))

            print(' Iter: ' + repr(iter) + ' z: ' + repr(z) + ' rc: ' + repr(rc_min) + ' ink: ' + repr(
                utgvar + 1) + ' utg: ' + repr(utgvar + 1))

            # make new partition
            # --------


elapsed = time.time() - t1
print('Elapsed time: ' + repr(elapsed))