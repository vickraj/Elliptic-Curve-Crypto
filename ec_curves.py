from ec_fields import f_add, f_mult, f_inv, f_sqrt, f_rand
from sympy import *
import numpy as np


###A module containing the two classes, Curves and Points, and related methods.

class Curve: #all curves initialize as the trivial elliptic curve.
    def __init__(self):
        self.q = 0; #q is the size of the finite field, 0 for rationals, -1 for reals.
        self.a = 0;
        self.b = 0;
        self.c = 0;
        
    ###Add: add two points on this elliptic curve. If the fac flag is on, then stop
    #as soon as a nontrivial gcd is found in any inversions.
    #The addition is done in the standard chord and flip method.
    ###Inputs:
    #x, y - points to add, should be on the elliptic curve.
    #fac - a flag determining whether we should be caring about gcd's.
    ###Outputs:
    #0 on failure.
    #p - a point that is the sum of x and y on this curve on success. 
    def add(self, x, y, fac =0):
        if ec.verify(x) != True:
            print("This is not on the curve.")
            return 0
        if ec.verify(y) != True:
            print("This is not on the curve.")
            return 0
        
        if y.z == 0:
            return x
        elif x.z == 0:
            return y
        else:
            q = self.q
            x_1 = x.plane()[0]
            y_1 = x.plane()[1]
            x_2 = y.plane()[0]
            y_2 = y.plane()[1]
            ret = Point()
            if x_1 != x_2:
                m = f_add(self.q, y_2, -y_1)
                m = f_mult(self.q, m, f_inv(self.q, f_add(self.q, x_2, -x_1), fac))
                x_3 = f_mult(q, m, m)
                x_3 = f_add(q, x_3, -x_1)
                x_3 = f_add(q, x_3, -x_2)
                y_3 = f_add(q, x_1, -x_3)
                y_3 = f_mult(q,m,y_3)
                y_3 = f_add(q, y_3, -y_1)
            elif y_1 != y_2:
                ret.x = 0
                ret.y = 1
                ret.z = 0
                return ret
            elif y_1 != 0: #P_1 = P_2, y_1 != 0
                m = f_mult(q, 3, x_1)
                m = f_mult(q, m, x_1)
                m = f_add(q, m, self.b)
                div = f_mult(q, 2, y_1)
                m = f_mult(q, m, f_inv(q, div, fac))
                #print(m)
                x_3 = f_mult(q, m, m)
                x_3 = f_add(q, x_3, -x_1)
                x_3 = f_add(q, x_3, -x_1)
                y_3 = f_add(q, x_1, -x_3)
                y_3 = f_mult(q,m,y_3)
                y_3 = f_add(q, y_3, -y_1)
                
            elif y_1 == 0:
                ret.x = 0
                ret.y = 1
                ret.z = 0
                return ret
            ret.x = x_3
            ret.y = y_3
            ret.z = 1
            return ret


    ###scale:
    #A way to easily scale a point on an elliptic curve by an integer.
    #We'll do double and add scaling.
    ###Inputs:
    #point - the point to scale
    #n - the integer to scale it by.
    #fac - an optional flag determining whether we should care about addition.
    ###Outputs:
    #0 on failure
    #ret - the scaled point. 
    def scale(self, point, n, fac = 0):#We'll do double and add scaling.
        if n == 0:
            return 0
        elif n == 1:
            return point
        elif n%2 == 1:
            return self.add(point, self.scale(point, n-1, fac), fac)
        else:
            return self.scale(self.add(point, point, fac), n/2, fac)
    ###rand
    #A method that returns a random point on this elliptic curve.
    ###Inputs: just the elliptic curve
    ###Outputs:
    #ret - a random point on this elliptic curve.
    
    def rand(self):#pick a random point on this ec. TODO: Prime power support.
        
        y_0 = False
        while(y_0 == False):
            x_0 = np.random.randint(0, self.q)
            y_square = f_mult(self.q, x_0, x_0)
            y_square = f_mult(self.q, x_0, y_square)
            temp = f_mult(self.q, x_0, self.b)
            y_square = f_add(self.q, temp, y_square)
            y_square = f_add(self.q, y_square, self.c)
            y_0 = f_sqrt(self.q, y_square)
        r = Point()
        r.x = x_0
        r.y = y_0
        r.z = 1
        return r

    
    ###verify:
    #A method that determines whether a point is actually on an elliptic curve.
    ###Inputs:
    #p - a point() object
    ###Outputs:
    #A boolean determining whether it is or is not on the curve.
    
    def verify(self, p):
        x = p.plane()[0]
        y = p.plane()[1]
        if self.q > 0:
            x = int(x)
            y = int(y)
            y = pow(y, 2, self.q)
            rhs = pow(x, 3, self.q)
            rhs += (self.b*x)%(self.q)
            rhs += self.c
            rhs = rhs%(self.q)
            if y == rhs:
                return True
            else:
                return False
        else:
            if y**2 == (x**3 + self.b*x + self.c):
                return True
            else:
                return False
        
    ###__str__
    #Debugging information that says exactly what this elliptic curve is.
    
    def __str__(self):
        stra = str(self.a);
        strb = str(self.b);
        if self.b < 0:
            strb = str(-1*self.b);
        strc = str(self.c);
        if self.c < 0:
            strc = str(-1*self.c)
        field = "";
        if self.q == 0:
            field = "rationals"
        elif self.q == -1:
            field = "reals"
        else:
            field = "finite field with "+str(self.q)+" elements"
        if self.a == 0:
            stra = ""
        elif self.a == 1:
            stra = "x^3"
        elif self.a == -1:
            stra = "-x^3"
        else:
            stra += "x^3"
        if self.b == 0:
            strb = ""
        elif self.b == 1 or self.b == -1:
            strb = "x"
        else:
            strb = strb +"x"
        string = stra;
        if self.a == 0:
            string = "";
        if self.b < 0 and self.a != 0:
            string += " - " + strb
        elif self.b > 0 and self.a != 0:
            string += " + " + strb
        elif self.a == 0 and self.b < 0:
            string += "-" + strb
        elif self.a == 0 and self.b > 0:
            string += strb
        if string == "":
            if self.c < 0:
                string += "-"+strc;
            else:
                string += strc;
        else:
            if self.c < 0:
                string += " - " + strc
            elif self.c > 0:
                string += " + " + strc
        return("This curve is over the " + field + " with equation" +
               " y^2 = "+string)
    
class Point: #all points initialize as the point at infinity, we'll
    #use projective coordinates.
    def __init__(self):
        self.x = 0;
        self.y = 1;
        self.z = 0;

    ###neg:
    #this will negate the point (as long as the curve is in weierstrass form.

    def neg(self):
        self.y = -self.y;
        if self.z == 0:
            self.y = 1
    ###plane:
    #This gives the points in plane coordinates, when otherwise they would be in
    #projective coordinates.
    ###Inputs: None
    ###Outputs:
    #the points in plane coordinates as a tuple. 
    def plane(self):
        if self.z == 0:
            return (0, 1)
        else:
            return (self.x/self.z, self.y/self.z)
    ###equiv_c
    #this allows points to check if they are equivalent to another point on a curve.
    ###Inputs:
    #ec - an elliptic curve.
    #m - the other point - default is the other point is the point at infinity.
    ###Outputs:
    #True or False, depending if the points are equivalent on this curve.
    def equiv_c(self, ec, m = 0):
        if m == 0:
            p = Point()
            p.x = 0
            p.y = 1
            p.z = 0
            return(self.equiv_c(ec, p))
        if self.z == 0 or m.z == 0:
            if self.z == 0 and m.z == 0:
                return True
            else:
                return False
        temp = f_inv(ec.q, m.z)
        s_temp = f_inv(ec.q, self.z)    
        self_x = f_mult(ec.q, (int) (self.x), s_temp)
        self_y = f_mult(ec.q, (int) (self.y), s_temp)
        m_x = f_mult(ec.q, m.x, m.z)
        m_y = f_mult(ec.q, m.y, m.z)
        if ec.q > 0:
            if m_x%(ec.q) == self_x%(ec.q) and (m_y%(ec.q)) == self_y%ec.q:
                return True
            else:
                return False
        else:
            if m_x == self_x and m_y == self_y:
                return True
            else:
                return False
        
    ###str:
    #the string method that gives the point in both plane coordinates and
    #projective coordinates. 
    def __str__(self):
        if self.z != 0:
            return "Plane: " + str((self.x/self.z, self.y/self.z)) + "\n" + "Projective: " + str((self.x, self.y, self.z))
        else:
            return "Point at infinity: (0, 1, 0)"
