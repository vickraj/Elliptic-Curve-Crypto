import ec_curves, ec_fields
import numpy as np
from sympy import sieve
from math import gcd
import importlib
importlib.reload(ec_curves)
importlib.reload(ec_fields)

#This module will be used for both factoring and for primality proving. 

def factor(n):
    #The classic ECM factoring. Pick a random curve mod n and a
    #random point on it, where
    #n is the number to factor, and then scale that point many times,
    #noting that we are doing gcd calculations during the
    #scaling. Then if the curve/point was randomly chosen correctly,
    #we'll get a non-trivial gcd - aka a factor of n.
    #Essentially we are just scaling this point and hoping it becomes
    #the identity modulo a factor of N.

    #Inputs: n, a number to factor
    #Outputs: a prime factor of n
    #TODO: make it exit nicely instead of breaking the program once it
    #finds a factor. 
    x_0 = np.random.randint(0, n)
    y_0 = np.random.randint(0, n)
    a = np.random.randint(0, n)
    b = ((y_0**2)%n - (x_0**3)%n - (a*x_0)%n)%n
    ec = ec_curves.Curve()
    ec.a = 1
    ec.b = a
    ec.c = b
    ec.q = n
    
    p = ec_curves.Point()
    p.x = x_0
    p.y = y_0
    p.z = 1
    
    b = 10
    s = ec_curves.Point()
    s.x = 0
    s.y = 1
    s.z = 0
    for i in range(1, b):
        p = ec.scale(p, i, 1)
        if p == True:
            return 1
        if p.equiv_c(ec, s):
            return 0
    print("iteration over")
    return 0

#This is the classic elliptic curve primality proving, based on the 
#Goldwasser-Kilian Theorem. If E is an elliptic curve over the
#rationals, with M and N positive and M > (N^1/4 + 1)^2 and N prime to
#the discriminant to E. IF there is a point P such that MP is zero mod
#N, and (M/l)P is strongly nonzero mod N for every prime l|M, then N
#is prime.
#Let us simplify for M=q to be prime.
#Then a primality certificate for p is (p, A, B, x_1, y_1, q), where P
#= (x_1, y_1, 1) is a point on E: y^2 = x^3 + Ax + B over Q, p is
#prime to disc(E), and q > (p^{1/4} + 1)^2 such that qP is zero mod p.
#This reduces p's primality to q's primality recursively.

#Inputs: p - a number to test primality.
# bnd - a bound for how "small" we allow c to be in cq = m to be.
#Outputs: prime or composite.

#TODO: Implement schoof's algorithm for counting points, and find an
#easy way to factor m to cq (or pick a q and then random curves until
#we get q to factor into m)
def isPrime(p, bnd):
    a = np.random.randint(0, p)
    x_0 = np.random.randint(0,p)
    y_0 = np.random.randint(0, p)
    b = y_0**2 - x_0**3 - a*x_0
    while gcd(4*(a**3) + 27*(b**2), p) != 1:
        a = np.random.randint(0, p)
        x_0 = np.random.randint(0,p)
        y_0 = np.random.randint(0, p)
        b = y_0**2 - x_0**3 - a*x_0
    x = ec_curves.Point()
    x.x = x_0
    x.y = y_0
    x.z = 1
    ec = ec_curves.Curve()
    ec.a = 1
    ec.b = a
    ec.c = b
    ec.q = p
    #m = schoof(ec)
    if m == False:
        return "composite"
    else:
        c = 1
        q = 1
        #wre m = cq, with c a small integer and q a probable prime.
        #if(miller_rabin(q)):
        #    return restart.
        point = ec.scale(x, c)
        if gcd(point.z, p) != 1:
            return "restart"
        else:
            x_0 = f_mult(p, point.x, f_inv(p, point.z))
            y_0 = f_mult(p, point.y, f_inv(p, point.z))
            qpoint = ec_curves.Point()
            qpoint.x = x_0
            qpoint.y = y_0
            qpoint.z = 1
            qpoint = ec.scale(qpoint, q)
            if qpoint.z%p != 0:
                return "composite"
            if q > b:
                ret = isPrime(q, b)
                while ret == "restart":
                    ret = isPrime(q, b)
                if ret == "composite":
                    return "restart"
                if ret == "prime":
                    return "prime"
            if q <= b:
                for i in range (2, q):
                    if q%i == 0:
                        return "composite"
                return "" + str(q) + "is prime"

###
#This is the miller-rabin algorithm for probabilistic primality
#testing.
#Essentially it tries to find a nontrivial square root of 1 as proof
#of compositeness, since only 1 and -1 are square roots of 1 in a
#prime field.

#Inputs: 
def miller_rabin(n):
    a = np.random.randint(1, n)
    t = n-1
    s = 0
    while t%2 == 0: #Write n-1 = 2^S*t, with t odd.
        t = (int) (t/2)
        s += 1
    print(t)
    print(s)
    b = 1
    for i in range (0, t): #calculate a^t
        b = b*a
        b = b%n
    print(b)
    if b%n == 1 or b%n == n-1: #if this is already a square root of 1,
        #it can't be a witness. 
        return True
    for i in range (1, s): #Otherwise check if it squared is -1 (maybe
        #b^2 can be a witness
        b = (b**2)%n
        if b%n == n-1:
            return True
    print("False: Witness: " + str(a) + " for " + str(n))#if neither
    #of these, than it is a witness for compositeness. 
    return False
