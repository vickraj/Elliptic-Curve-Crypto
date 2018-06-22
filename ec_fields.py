import numpy as np
import sympy as sp

def poly_find(q):#TODO: hardcode the polynomials we'll use for
    #field extensions for nonprime q's.
    return "";

#f_add: Add two numbers in a prime field. 
def f_add(q, a, b):
    if q <= 0:
        su = a+b
    else:
        su = (int) ((a%q + b%q)%q)
    return su

#multiply two numbers naively in a prime field. 
def f_mult(q, a, b): #TODO: Implement karatsuba or FFT for faster integer multiplication. 
    if q <= 0:
        ret = a*b
        return ret
    else:
        ret = (a%q * b%q)%q
        return (int)(ret)
    return ret

#Get a random element in a prime power field. 
def f_rand(q): #TODO: when prime power field support happens, get a
    #random element. 
    t = sp.symbols('t')
    idea = sp.poly(t**2 + 1, t)
    r_1 = np.random.randint(0,3)
    r_2 = np.random.randint(0,3)
    return r_1 + r_2*t

#invert an element in a prime field. 
def f_inv(q, a, fac=0): #extended euc alg. while we are working with just
    #primes TODO: extended poly euclidean algo when we allow for
    #nonprimes as well.
    if q <= 0:
        if a != 0:
            return 1/a
        else:
            return "Error: Inverting 0"
    else:#q is prime right now.
        if a == 0:
            return "Error: Inverting 0"
        q_mat = np.matrix([[0, 1],[1, 0]])
        r_mat = np.matrix([[0],[0]])
        s_mat = np.matrix([[1],[0]])
        t_mat = np.matrix([[0],[1]])
        r_mat[0] = q
        r_mat[1] = a
        q_mat[1,1] = (int) (-r_mat[0]/r_mat[1])
        while r_mat[1] != 0:
            r_mat = q_mat*r_mat
            s_mat = q_mat*s_mat
            t_mat = q_mat*t_mat
            if r_mat[1] != 0:    
                q_mat[1,1] = (int)(-r_mat[0]/r_mat[1])
                #print(q_mat)
        gcd = s_mat[0]*q + t_mat[0]*a
        if gcd != 1 and gcd != q and fac == 1:
            print("Factor found: "+ str(gcd) + " of " + str(q))
            quit()
            
        return (int) (t_mat[0]%q)

def f_sqrt(p, n):
    #Calculate a solution to x^2 cong n mod p, if a solution exists.
    #Uses the Tonelli-Shanks algorithm, which is polynomial if the
    #generalized riemann hypothesis is true.
    
    leg = pow(n, (int) ((p-1)/2), p) 
    if leg != 1:
        return False
    if n == 0:
        return 0
    if p == 2:
        return p
    if p%4 == 3:
        return pow(n, (int) ((p+1)/4), p)

    s = 0
    t = p-1
    while t%2 == 0:
        s += 1
        t = (int) (t/2)
    leg = 1
    while leg == 1:
        a = np.random.randint(0, p)
        leg = pow(a, int((p-1)/2), p)
    M = s
    c = pow(a, t, p)
    r = pow(n, (int) ((t + 1)/2), p)
    t = pow(n, t, p)
    while t != 0 and t != 1:
        su = t
        i = 0
        while su != 1:
            su = pow(su, 2, p)
            i += 1
        b = pow(c, pow(2, M - i - 1), p)
        M = i
        c = pow(b, 2, p)
        t = t*c%p
        r = r*b%p
        
    if t == 0:
        return 0
    if t == 1:
        return r
        
