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

diffx=0
diffz=0
# Start iterations
iter = 0
while iter >= 0:
    iter += 1  

    B_inv = np.linalg.inv(B)
    xB = np.dot(B_inv, b)

    y = np.transpose(np.dot(np.transpose(cB),B_inv))
    z = np.dot(np.transpose(b),y)
    c_hattN = cN - (np.dot(np.transpose(N),y)) #np.transpose(N)

    # calc most negative reduced cost, rc_min,
    # and index for entering variable, inkvar
    # --------

    inkvar = np.argmin(c_hattN)
    rc_min = c_hattN[inkvar]

    

    if(c_hattN >=0).all():
        iter=-1
    if rc_min >= -1.0E-12:
        print('Ready')
        iter = -1


            # construct solution, x, and check it
            # --------
#        diffx = np.linalg.norm(x - xcheat)
 #       diffz = z - zcheat
     #   print('xdiff: ' + repr(diffx))
      #  print('zdiff: ' + repr(diffz))
    else:
        # calc entering column, a
        # --------
        # a = inkommande kolumn hatt
        # A[:,inkvar] är att ta ut kolumn ur A på plats inkvar
        a = np.dot(B_inv, N[:, inkvar]) #A ist för N

        if max(a) <= 0:
            # unbounded solution
            print('Unbounded solution!')
            iter = -1
        else:
            # calc leaving var, utgvar
            # --------
            ut_list = []

            for a_rows, b_rows in zip(N, xB):
                if a_rows[inkvar] is not 0:
                    ut_list.append(b_rows / a_rows[inkvar])

            utgvar = ut_list.index(min(ut_list))

            print(' Iter: ' + repr(iter) + ' z: ' + repr(z) + ' rc: ' + repr(rc_min) + ' ink: ' + repr(
                utgvar + 1) + ' utg: ' + repr(utgvar + 1))

            # make new partition
            # --------

            #bix[utgvar] = inkvar


            bix[utgvar], nix[inkvar] = nix[inkvar], bix[utgvar]
            temp = np.copy(B[:, utgvar])
            B[:, utgvar] = N[:, inkvar]
            N[:, inkvar] = temp
            cN[inkvar], cB[utgvar] = cB[utgvar], cN[inkvar]
            # N = np.roll(N,1,1)
            # nix = np.roll(nix,1,0)
            #cN = np.roll(cN,1,0)
          #  c = np.append(cN,cB)


elapsed = time.time() - t1
print('Elapsed time: ' + repr(elapsed))
print(zcheat)
print("z:",z)