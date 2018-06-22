import ec_curves
from ec_fields import f_mult, f_add, f_inv
import sympy as sp
import importlib

importlib.reload(ec_curves)

#This module will cover most of the methods necessary to compute the
#Weil Pairing on a given elliptic curve. Currently it only covers
#curves over the reals or prime fields. May also add some code to
#compute functions at a divisor to show weil reciprocity as well. 




#Given P, Q as m-torsion points on an elliptic curve, let f_P, f_Q be
#rational functions on the curve such that div(f_P) = m(P) - m(O), and
#div(f_Q) = m(Q) - m(O) - these functions must exist since the RHS is
#a principal divisor (degree is 0 and sum is O).

#The Weil pairing is a bilinear pairing defined as [f_P(Q +
#S)/f_P(S)]/[f_Q(P - S)/f_Q(-S)], where S is some point on the curve
#to make this defined and non-zero. 
def weil(ec, p, q, m, s): #s is any point on the curve not P, -Q, P-Q, O

    #Verify p, q are m-torsion points, and on curve
    if ec.verify(p) == False or ec.verify(q) == False:
        return "At least one point is not on the curve"
    
    mp = ec.scale(p, m)
    mq = ec.scale(q, m)
    if mp.equiv_c(ec) != True:
        return "P is not an m-torsion point"
    if mq.equiv_c(ec) != True:
        return "Q is not an m-torsion point"
    
    x,y = sp.symbols('x y')
    f_pn, f_pd = millers_algorithm(ec, p, m)
    f_qn, f_qd = millers_algorithm(ec, q, m)
    found = False
    while found == False:#Find an s such that s notin {O, P, -Q, P-Q} to prevent
        #degeneracy of the weil computations.
        s = ec.rand()
        q.neg()
        if s.equiv_c(ec):
            found = False
        elif s.equiv_c(ec, p):
            found = False
        elif s.equiv_c(ec, q):
            
            found = False
        elif s.equiv_c(ec, ec.add(p, q)):
            found = False
        else:
            found = True
        q.neg()
    
    ps = ec.add(q, s)
    neg = ec_curves.Point()
    neg.x = s.x
    neg.y = -s.y
    neg.z = s.z
    qs = ec.add(p, neg)
    
    num = f_pn.subs([(x, ps.plane()[0]), (y, ps.plane()[1])])
    
    num = f_mult(ec.q, num, f_inv(ec.q, f_pd.subs([(x, ps.plane()[0]),
    (y, ps.plane()[1])])))

    den = f_pn.subs([(x, s.plane()[0]), (y, s.plane()[1])])
    den = f_mult(ec.q, den, f_inv(ec.q, f_pd.subs([(x, s.plane()[0]),
    (y, s.plane()[1])])))

    
    num2 = f_qn.subs([(x, qs.plane()[0]), (y, qs.plane()[1])])
    num2 = f_mult(ec.q, num2, f_inv(ec.q, f_qd.subs([(x,
    qs.plane()[0]), (y, qs.plane()[1])])))

    den2 = f_qn.subs([(x, neg.plane()[0]), (y, neg.plane()[1])])
    den2 = f_mult(ec.q, den2, f_inv(ec.q, f_qd.subs([(x,
    neg.plane()[0]), (y, neg.plane()[1])])))
    
    num = f_mult(ec.q, num, den2)
    den = f_mult(ec.q, den, num2)
    ret = f_mult(ec.q, num, f_inv(ec.q, den))
    
    print("Weil: " + str(ret))
    return(ret)


###h_func:
#A helper function used by Millers Algorithm.
###Inputs:
#p, q - two points 
###Outputs:
#a rational function h_p, denoted in terms of numerators and denominators.
#So it returns h_np and h_dp, where h_np is the numerator and h_dp is the
#denominator of the rational function.

def h_func(ec, p, q):

    p_x = p.plane()[0]
    p_y = p.plane()[1]
    q_x = q.plane()[0]
    q_y = q.plane()[1]
    x, y = sp.symbols('x y')
    if p_x == q_x and p_y != q_y:
        h_pqn = x - p_x
        h_pqd = 1
    else:
        if p_x == q_x and p_y == q_y and p_y == 0:
            h_pqn = x - p_x
            h_pqd = 1
        else:
            if p_x == q_x and p_y == q_y:
                m = (3*p_x*p_x + ec.b)/(2*p_y)
                m = f_mult(ec.q, p_x, p_x)
                temp = f_add(ec.q, m, m)
                m = f_add(ec.q, m, temp)
                m = f_add(ec.q, m, ec.b)
                temp = f_add(ec.q, p_y, p_y)
                m = f_mult(ec.q, m, f_inv(ec.q, temp))
                
                #print("yay")
            else:
                m = f_add(ec.q, q_y, -p_y)
                temp = f_add(ec.q, q_x, -p_x)
                #print("this is m")
                #print(m)
                #print("this is inv")
                #print(f_inv(ec.q, temp))
                m = f_mult(ec.q, m, f_inv(ec.q, temp))
            temp = f_add(ec.q, p_x, q_x)
            temp1 = f_mult(ec.q, m, m)
            temp = f_add(ec.q, temp, -temp1)
            h_pqn = (y - p_y - m*(x -p_x))
            h_pqd = (x + temp)
    
    return h_pqn, h_pqd

###millers:
#This is millers algorithm to come up with a function whose divisor is m(P) - m(O).
###Inputs:
#p - an m-torsion point
#m - the integer determining which torsion subgroup we're in.
###Outputs:
#f_p - a function whose divisor is m(P) - m(O)
def millers_algorithm(ec, p, m):
    t = p
    f_n = 1
    f_d = 1
    x,y = sp.symbols('x y')
    #verify p is an m torsion point on ec
    bin_m = bin(m)[2:]
    #print(bin_m)
    i = len(bin_m) - 2
    #print(i)
    while i >= 0:
        f_pn, f_pd = h_func(ec, t, t)
        f_n = f_n**2 * f_pn
        f_d = f_d**2 * f_pd
        t = ec.add(t, t)
        if bin_m[i] == '1':
            f_pn, f_pd = h_func(ec, t, p)
            f_n = f_n * f_pn
            f_d = f_d *f_pd
        
            t = ec.add(t, p)
        f_n = f_n%(ec.q)
        f_d = f_d%(ec.q)
        i = i - 1
    return f_n, f_d
    
