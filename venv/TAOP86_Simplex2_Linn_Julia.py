import numpy as np
import time
import sys
from copy import copy

filc = " ".join(sys.argv[1:]).split('.')[0] + '.npz'
npzfile = np.load(filc)
c = npzfile['c'].astype(float)
b = npzfile['b'].astype(float)
A = npzfile['A'].astype(float)
bix = npzfile['bix']
zcheat = npzfile['zcheat']
xcheat = npzfile['xcheat']
bounded = True
bix = bix - 1
t1 = time.time()
[m, n] = np.shape(A)

# Create nix
nix = np.setdiff1d(range(n), bix)
B = A[:, bix]
N = A[:, nix]
cB = c[bix]
cN = c[nix]
x_B= None

#Start iterations
iter = 0
while iter >= 0:
    iter += 1
    # calc right-hand-sides and reduced costs
    # --------
    B_inverse = np.linalg.inv(B)
    # calculate x_B
    x_B = np.dot(B_inverse, b)
    # calculate "skuggpris", dualvariabel y
    cB_T = np.transpose(cB)
    y = np.transpose(np.dot(cB_T, B_inverse))
    b_T = np.transpose(b)
    # målfunktion
    z = np.dot(b_T, y)

    #calculate Cn-hatt
    cN_hat = cN - np.dot(np.transpose(N), y)

    # If Cn-hatt are all positive, we are done, end loop
    rc_min = -0.0000000000001
    if (cN_hat >= -0.000000001).all():
        iter = -1

    else:
        # calculate inkommande variabel, inkvar
        temp_inkvar = -1
        rc_min = -0.0000000000001

        for rc in cN_hat:
            temp_inkvar += 1
            if rc < rc_min:
                rc_min = rc
                inkvar = temp_inkvar


        a = np.dot(B_inverse, N[:, inkvar])

        if (a <= 0).all():
            # if unbounded solution, end loop
            print('Unbounded solution!')
            z = 'Inf'
            bounded = False
            iter = -1

        else:
            # calc leaving var, utgvar
            # --------
            in_col_divided = []

            for q in range(x_B.size):
                if a[q] != 0:
                    in_col_divided.append(x_B[q]/a[q])
                else:
                    in_col_divided.append(100000)

            outgoing_var = 100000000

            temp_index_outgoing = -1

            for d in in_col_divided:
                temp_index_outgoing += 1
                if (d < outgoing_var) and (d > 0):
                    outgoing_var = d
                    utgvar = temp_index_outgoing

            temp1 = copy(nix[inkvar])
            temp2 = copy(bix[utgvar])
            nix[inkvar] = temp2
            bix[utgvar] = temp1

            temp = copy(B[:, utgvar])
            B[:, utgvar] = N[:, inkvar]
            N[:, inkvar] = temp

            temp = copy(cN[inkvar])
            cN[inkvar] = cB[utgvar]
            cB[utgvar] = temp


x = []
for t in range(n):
    if t in bix:
        x.append(x_B[np.where(bix == t)][0])
    else:
        x.append(0)



elapsed = time.time() - t1

print('Elapsed time: ' + repr(elapsed))
print(" ")
print("RESULT(z and x): ")
print(z)
#print(x)
print(" ")

print("CORRECT: ")
print(zcheat)
#print(xcheat)

if bounded:
    print("DIFFERENCE z, x")
    print( np.round(z, 3) - np.round(zcheat, 3))
