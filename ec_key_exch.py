import ec_weil
import ec_curves
import numpy as np
from ec_fields import f_mult, f_add, f_inv
import importlib
importlib.reload(ec_weil)
importlib.reload(ec_curves)

#This module will contain all the key exchange algorithms.
#TODO: Implement tate pairing and one round Tri DH with tate pairing.


#This is an example of one round tripartite diffie hellman key
#exchange, following Joux's idea of using pairing based methods.
###Inputs:
#ec - an elliptic curve(an ec_curves.Curve() object)
#P - one of the public points. (ec_curves.Point() object)
#Q - another public point.
#n - P and Q should both be [n]-torsion points.
###Outputs:
#Several strings of debugging:
#The broadcasted points by A, B, and C
#The random integers each of A, B, and C chose
#The outputs of the Weil pairing for each of them.
#The final shared key that each of them computed. These should be the
#same.

def tri_dh(ec, P, Q, n):
    if ec.verify(P) != True or ec.verify(Q) != True:
        return "These points were not on the curve"
    nP = ec.scale(P, n)
    nQ = ec.scale(Q, n)
    if nP.equiv_c(ec) == False:
        return "P is not an n-torsion point"
    if nQ.equiv_c(ec) == False:
        return "Q is not an n-torsion point"
    
    a = np.random.randint(1, 100)
    b = np.random.randint(1, 100)
    c = np.random.randint(1, 100)
    
    P_a = ec.scale(P, a)
    Q_a = ec.scale(Q, a)
    print("A broadcasts P_a:\n" + str(P_a) + "\n and Q_a: " + str(Q_a))
    P_b = ec.scale(P, b)
    Q_b = ec.scale(Q, b)
n    print("B broadcasts P_b:\n" + str(P_b) + "\n and Q_b: " + str(Q_b))
    P_c = ec.scale(P, c)
    Q_c = ec.scale(Q, c)
    print("C broadcasts P_c:\n" + str(P_c) + "\n and Q_c: " + str(Q_c))
    print(str(a) + " A: " + str(F_w(ec, a, P_b, Q_c, n)))
    print(str(b) + " B: " + str(F_w(ec, b, P_a, Q_c, n)))
    print(str(c) + " C: " + str(F_w(ec, c, P_a, Q_b, n)))
    
###F_w: This is the actual computations that each of A, B, and C do in
###order to make the shared secret. It is simply just computing
### e_n(p, q)^x, where e_n is the weil pairing on n-torsion points.

###Inputs:
#ec - an elliptic curve
#x - an integer to scale the weil pairing.
#p, q - n-torsion points.
#n - the integer specifying which torsion subgroup we're working in.
###Outputs:
#w - the value of e_n(p,q)^x

def F_w(ec, x, P, Q, n):
    s = ec_curves.Point()
    s.x = 0
    s.y = 0
    s.z = 1#note that this will only work when 0,0 is on the
    #curve. usually we can just pick random points.
    w = ec_weil.weil(ec, P, Q, n, s)
    w = w**x%(ec.q)
    return w
