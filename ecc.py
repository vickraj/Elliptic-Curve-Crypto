import ec_curves
import ec_fields
import ec_weil
import ec_disc
import ec_primes
import ec_dsa
import ec_key_exch
import importlib
import sympy as sp
import numpy as np
importlib.reload(ec_fields)
importlib.reload(ec_curves)
importlib.reload(ec_weil)
importlib.reload(ec_disc)
importlib.reload(ec_primes)
importlib.reload(ec_dsa)
importlib.reload(ec_key_exch)
def main():
    ###Some examples below, and then the terminal code. 
    ec = ec_curves.Curve();
    ec.a = 1
    ec.b = 37
    ec.c = 0
    ec.q = 1009
    print("Welcome to Elliptic Curve Crypto Suite!")
    print(ec);
    x = ec_curves.Point();
    x.x = 417
    x.y = 952
    x.z = 1
    y = ec_curves.Point();
    y.x = 561
    y.y = 153
    y.z = 1
    s = ec_curves.Point()
    s.x = 0
    s.y = 0
    s.z = 1
    su = ec.add(x, y)
    print(su)
    m = 7
    a,b = sp.symbols('a b')
    
    
    print((int) (ec_weil.weil(ec, x, y, m, s)))
    ec.a = 1
    ec.b = 0
    ec.c = 1
    ec.q = 599
    s.x = 60
    s.y = 19
    s.z = 1
    y.x = 277
    y.y = 239
    y.z = 1
    print(ec_disc.bsgs(ec, s, y, 600))
    ec.a = 1
    ec.b = 0
    ec.c = 1
    ec.q = 9
    ret = 0
    #while ret == 0:
        #ret = ec_primes.factor(38861*38867)
    i = 0
    ret = True
    n = 3886123
    while i < 10 and ret == True:
        ret = ec_primes.miller_rabin(n)
        i += 1
    if ret == False:
        print(str(n) + " is definitely composite")
    else:
        print(str(n) + " probably prime")

    print("ECDSA:")
    ec.b = 2
    ec.c = 2
    ec.q = 17
    G = ec_curves.Point()
    G.x = 5
    G.y = 1
    G.z = 1
    priv_k = 3
    pub_k = ec.scale(G, priv_k)
    Q = pub_k
    #r = schoof(ec)
    r = 19
    m = 10935823
    m, rpoint, s = ec_dsa.sign(ec, r, G, priv_k, m)
    print(str(m) + " verification: " + str(rpoint) + " "+str(s))
    m += 2
    print(ec_dsa.verify(ec, r, G, pub_k, m, rpoint, s))

    print("One Round Tripartite Diffie Hellman Key Exchange:")
    
    m = 7
    x = ec_curves.Point();
    x.x = 417
    x.y = 952
    x.z = 1
    y = ec_curves.Point();
    y.x = 561
    y.y = 153
    y.z = 1
    ec.b = 37
    ec.c = 0
    ec.q = 1009
    ec_key_exch.tri_dh(ec, x, y, m)

    ###Everything above is just examples. 
    print("Welcome to Elliptic Curve Crypto Suite:")
    p = ec_curves.Point();
    q = ec_curves.Point();
    while True:
        print(str(ec))
        st = input("You can [c]hange curve \n" +
               "[f]actor an integer \n" +
              "check [p]rimality \n" +
              "Calculate [d]iscrete logs on this curve \n" +
              "Do the [w]eil pairing of two points \n" +
              "see an example of one-round tripartite diffie-hellman key e[x]change \n" +
              "Do [e]cdsa \n" + 
              "[a]dd points \n" +
              "[s]cale points \n" +
              "get a [r]andom point \n >>>")
        if st == 'c':
            print("Weierstrass form: y^2 = x^3 + ax + b")
            ec.b = int(input("a?"))
            ec.c = int(input("c?"))
            ec.q = int(input("What field is it over? (only prime fields)"))
            i = 0n
            ret = True
            prime = False
            while prime == False:
                while i < 10 and ret == True:
                    ret = ec_primes.miller_rabin(ec.q)
                    i +=1 
                if ret == True:
                    prime = True
                else:
                    ec.q = input("That was definitely not prime. Try again.")
        if st == 'f':
            n = input("Number to factor?")
            while n != 'b':
                ret = True
                i = 0
                while i < 10 and ret == True:
                    ret = ec_primes.miller_rabin(int(n))
                    i += 1
                if ret == True:
                    n = input("This number is probably prime. Try again, or [b]ack")
                else:
                    ec_primes.factor(int(n))
        if st == 'p':
            n = input("Number to test primality?")
            i = 0
            ret = True
            while i < 10 and ret == True:
                ret = ec_primes.miller_rabin(int(n))
                i += 1
            if ret == True:
                print("Probably Prime")
            else:
                print("Definitely composite")
        if st == 'd':
            p = ec_curves.Point()
            q = ec_curves.Point()
            p.x = int(input("Discrete log is aP = Q. What's the x coordinate of P?"))
            p.y = int(input("y coordinate of P?"))
            p.z = 1
            q.x = int(input("x coordinate of Q?"))
            q.y = int(input("y coordinate of Q?"))
            q.z = 1
            n = input("Group order of elliptic curve?")
            #n = schoof(ec)
            print(ec_disc.bsgs(ec, p, q, n))
        if st == 'w':
            p = ec_curves.Point()
            q = ec_curves.Point()
            p.x = int(input("The Weil Pairing requires two [n]-torsion points, P and Q." +
                            "\nWhat's the x-coordinate of P?"))
            p.y = int(input("y coordinate of P?"))
            p.z = 1
            q.x = int(input("x coordinate of Q?"))
            q.y = int(input("y coordinate of Q?"))
            q.z = 1
            n = int(input("Value of n?"))
            s = ec_curves.Point()
            s.x = 0
            s.y = 0
            s.z = 1
            print(ec_weil.weil(ec, p, q, n, s))
        if st == 'x':
            p.x = int(input("Tri DH requires two independent [n]-torsion points, P and Q." +
                            "\n What's the x-coordinate of P?"))
            p.y = int(input("y coordinate of P?"))
            p.z = 1
            q.x = int(input("x coordinate of Q?"))
            q.y = int(input("y coordinate of Q?"))
            q.z = 1
            n = int(input("[n]?"))
            print(ec_key_exch.tri_dh(ec, p, q, n))
        if st == 'e':
            g = ec_curves.Point()
            g.x = int(input("ECDSA requires a generator point. x coordinate?"))
            g.y = int(input("y coordinate?"))
            g.z = 1
            priv_k = int(input("ECDSA requires a private key to scale this point."))
            r = int(input("r such that #E = fr, with r a large prime and f small."))
            pub_k = ec.scale(g, priv_k)
            Q = pub_k
            #r = schoof(ec)
            r = 19
            sig = input("[s]ign or [v]erify?")
            if sig == 's':
                m = int(input("message to sign?"))
                m, rpoint, s = ec_dsa.sign(ec, r, G, priv_k, m)
                print(str(m) + " verification: " + str(rpoint) + " "+str(s))
            if sig == 'v':
                m = int(input("message?"))
                rpoint = ec_curves.Point()
                rpoint.x = int(input("x coordinate of rpoint?"))
                rpoint.y = int(input("y coordinate of rpoint?"))
                rpoint.z = 1
                s = int(input("s?"))
                print(ec_dsa.verify(ec, r, G, pub_k, m, rpoint, s))
        if st == 'a':
            p.x = int(input("Adding requires two points, P and Q." +
                            "\n What's the x-coordinate of P?"))
            p.y = int(input("y coordinate of P?"))
            p.z = 1
            q.x = int(input("x coordinate of Q?"))
            q.y = int(input("y coordinate of Q?"))
            q.z = 1
            print(ec.add(p, q))
        if st == 's':
            
            p.x = int(input("scaling requires one point, P and an integer k." +
                            "\n What's the x-coordinate of P?"))
            p.y = int(input("y coordinate of P?"))
            p.z = 1y
            k = int(input("How much to scale this by?"))
            print(ec.scale(p, k))
        if st == 'r':
            print("A random point on this curve:")
            p = ec.rand()
            print(p)
            print(ec.verify(p))
if __name__== '__main__':
    main()
