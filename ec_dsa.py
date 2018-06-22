import ec_curves
from ec_fields import f_mult, f_inv, f_add
import numpy as np
import importlib
importlib.reload(ec_curves)

#This module covers Signature algorithms, currently just ECDSA
#TODO: Get ECIES in this as well.



#ECDSA: 
#pick curve, pick finite field, such that #E = fr, with r a large
#prime and f a small integer. Pick a base point G in E(Fq) of order
#r. calculate Q = aG for a secret a, and publish Q, G, r, E, Fq

###sign:
#Given the public information, and a secret integer a, sign a message
#m.
###Inputs:
#ec - the elliptic curve
#r - A large prime such that #E = fr, with f small.
#G - a generator point for a subgroup of order r
#a - a secret integer to scale G by.
#m - the message to sign.
###Outputs:
#m - the message again. 
#rpoint - A point on the curve as part of the verification.
#s - the signature. 
def sign(ec, r, G, a, m):
    k = np.random.randint(1, r)
    rpoint = ec.scale(G, k)
    s = f_inv(r, k)
    temp = f_mult(r, a, rpoint.x)
    temp = f_add(r, m, temp)
    s = f_mult(r, s, temp)
    return m, rpoint, s

###verify:
#Do computations with the signature and check if it matches the
#rpoint.
###Inputs:
#ec - the elliptic curve
#r - the same r as above
#G - the same G as above
#Q - the same public Q as above
#rpoint - the point to verify with.
#s - the signature

def verify(ec, r, G, Q, m, rpoint, s):
    u_1 = f_mult(r, f_inv(r, s), m)
    u_2 = f_mult(r, f_inv(r, s), rpoint.x)
    v = ec.scale(G, u_1)
    temp = ec.scale(Q, u_2)
    v = ec.add(v, temp)
    if v.equiv_c(ec, rpoint):
        return True
    else:
        return False
