import ec_weil, ec_curves, ec_fields
import importlib
import sympy as sp
importlib.reload(ec_weil)
importlib.reload(ec_curves)
importlib.reload(ec_fields)


#This module will be used for general attacks on discrete logs.
#In particular, it covers Baby Steps Giant steps as a general method
#of computing the discrete log on elliptic curves, and the MOV attack as an example,
#though there is currently not support for solving the discrete log in
#prime power fields (index calculus methods)

#TODO: Cover index calculus methods.
#TODO: Cover Pollards Rho and Lambda Algorithms
#TODO: Cover Pohlig Hellman method
#Solve the dl problem q = kp, with n the group order size.
def bsgs(ec, p, q, n): 
    m = (int) (n**(1/2))+ 1
    li = []
    ip = ec_curves.Point()
    ip.x = 0
    ip.y = 1
    ip.z = 0
    li.append(ip)
    mp = ec_curves.Point()
    mp = ec.scale(p, m)
    print(mp)
    for i in range (0, m):
        ip = ec.add(p, ip)
        li.append(ip)
    for i in range(0, m):
        #print("i:" + str(i))
        #print(li[i])
        if li[i].equiv_c(ec, q):
            print("found i! " + str(i))
    qp = ec_curves.Point()
    qp = q
    mp.neg()
    for j in range(1, m):
        #print(qp)
        qp = ec.add(qp, mp)
        for i in range(0, m):
            if li[i].equiv_c(ec, qp):
                #print(qp)
                print("found i! " + str(i) + " j:" + str(j) + " m:" + str(m))
                return i + m*j
    return "can't find discrete log"


def mov(ec, p, q):
    #Only going to do in the case of supersingular curves.
    # Examples the curves y^2 = x^3 - ax over Fp, with p > 3 a prime a a
    # quadratic nonresidue in F_p, and p cong 3 mod 4.
    #for p == 13, the quadratic residues are 1, 3, 4, 9, 10, 12
    #TODO: Requires computations in F_p^2, which hasn't been implemented
    #yet.
    #This follows from the original MOV paper, which is not public afaik. 
    a = ec_curves.Point()
    a.x = 0
    a.y = 1
    a.z = 0
    temp = ec.add(ec, p, a)
    n = 1
    while not a.equiv_c(temp):
        temp = ec.add(ec, temp, p)
        n += 1
    #n is the order of the point P. 
    if quad_res(ec.b, ec.q):
        k = 2
        c = 2
        n_1 = ec.q + 1
    else:
        k = 2
        c = 1
        n_1 = (int) ((ec.q + 1)/2)
    ec.q = ec.q**k
    rand_p = ec.rand() #TODO 
    rand_p = ec.scale(rand_p, (c*n_1)/n)
    s= ec_curves.Point() #The point 0, 0 will work for the class 1 and
    #2 supersingular curves we're working with here. 
    s.x = 0
    s.y = 0
    s.z = 1
    a = ec_weil.weil(ec, p, rand_p, n, s)
    b = ec_weil.weil(ec, q, rand_p, n, s)
    #finite field discrete log of al' = b in F_q^k
    print("The discrete log problem in F_q^k is " + str(a) + "x = " +
          str(b))
    
def quad_res(a, p): #check whether a number is a quad res in F_p.
    #essentially it just calculates the legendre symbol for F_p, and will
    #be in more detail when this expands to F_p^d
    euler =(int) ( a**((p-1)/2))
    if euler%p == 1:
        return True
    else:
        return False
    
    
