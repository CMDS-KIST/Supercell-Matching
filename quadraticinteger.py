from eisenstein.eisenstein import EisensteinInt
from gaussint import GaussInt
from sympy.solvers.diophantine.diophantine import diop_DN, transformation_to_DN, find_DN
from sympy import Matrix
from sympy.ntheory import factorint
from sympy.core import integer_nthroot
from itertools import chain, combinations
import cmath, math

class EisensteinInt_hex(EisensteinInt):
    def __init__(self, real, complex):
        super(EisensteinInt_hex, self).__init__(real, complex)

    def __str__(self):
        a = self.real
        b = self.imaginary

        if b == 0:
            return str(a)
        elif (b > 0):
            s = ""
            op = "+"
        else:
            s = "-"
            op = "-"

        if (abs(b) == 1):
            b_ = "\\omega_{{hex}}"
        else:
            b_ = "{}\\omega_{{hex}}".format(str(abs(b)))

        if a == 0:
            result = "{} {}".format(s, b_)
        else:
            result = "{} {} {}".format(str(a), op, b_)

        return result.strip()
    
    def arg(self):
        return cmath.phase(self.complex_form())

class GaussInt_square(GaussInt):
    def __init__(self, a, b):
        super(GaussInt_square, self).__init__(a, b)
    
    def __str__(self):  # Overload string conversion used by print
        return "" + str(self.r) + ((" + "+str(self.i)) if (self.i >= 0) else (" - "+str(-self.i))) + "\\omega_{{sq}}"
    
    def arg(self):
        return cmath.phase(self.r + self.i*1j)

    def is_unit(self):
        flag_unit = False
        if self.r**2 + self.i**2 == 1:
            flag_unit = True
        return flag_unit

    def gcd(self, b):
        v = super(GaussInt_square, self).gcd(b)
        return GaussInt_square(v.r, v.i)

from sympy import isprime, primefactors, factorint, solve, symbols, Symbol, Eq
from math import sqrt, pi, sin, cos, atan
import matplotlib.pyplot as plt
import numpy as np
import mpmath

class EuclideanDomainInt:
    """
    Stores the quadratic integer of Euclidean Doomain in the form a + bω

    Create Euclidean Integer by:
	     a = EuclideanDomainInt(5,7)  # Create (5 + 7ω)
	     a = EuclideanDomainInt(13)  # Create (13)

    Functions implemented
         Basic functions: init(), ==, hash(),  str(), <, >, <=, >=
	     Arithmetic functions: abs(), +, divmod(), //, %, *, -

         a.debug_str() - Returns a simple debug string to describe the
            data members.

         a.associates() - Returns a list of the product of the number and
            each of the units.
         a.canonical() - Returns the associate in the first sextant.
         a.complex_form() - Returns the complex form.
         a.conjugate() - Returns an EisensteinInt representing he conjugate.
         a.dot_product(b) - Returns the dot product.
         a.norm() - Returns an integer representing the norm.
         a.polar_form() - Returns a tuple of the radius and angle.

         a.gcd(b) - Compute the greatest common divisor of a and b.
    """

    def __init__(self, real=0, imaginary=0, d=-1):
        assert isinstance(real, int)
        assert isinstance(imaginary, int)
        assert isinstance(d, int)

        self.real = real
        self.imaginary = imaginary
        self.d = d
        if d % 4 != 1:
            self.omega = cmath.sqrt(d)
        else:
            self.omega = (-1.0 + cmath.sqrt(d))/2.0

    def __hash__(self):
        return hash((self.real, self.imaginary))

    def __eq__(self, other):
        if isinstance(other, EuclideanDomainInt):
            return (self.real == other.real) and (self.imaginary == other.imaginary) and (self.d == other.d)
        else:
            return False

    def __str__(self):
        a = self.real
        b = self.imaginary

        if b == 0:
            return str(a)
        elif (b > 0):
            s = ""
            op = "+"
        else:
            s = "-"
            op = "-"

        if (abs(b) == 1):
            b_ = "ω"
        else:
            b_ = f"{str(abs(b))}ω"

        if a == 0:
            result = f"{s} {b_}"
        else:
            result = f"{str(a)} {op} {b_}"

        return result.strip()
    
    def __repr__(self):
        return f"EuclideanDomainInt({self.real}, {self.imaginary}, {self.d})"

    def debug_str(self):
        result = f"EuclideanDomeinInt({self.real}, {self.imaginary}, {self.d})"
        return result

    def __abs__(self):
        return sqrt(self.norm())

    def __add__(self, other):
        if isinstance(other, int):
            other = EuclideanDomainInt(other)

        a = self.real
        b = self.imaginary
        c = other.real
        d = other.imaginary

        sum_real = a + c
        sum_imaginary = b + d

        return EuclideanDomainInt(sum_real, sum_imaginary)

    def __divmod__(self, other):
        if isinstance(other, int):
            other = EuclideanDomainInt(other, 0, self.d)

        assert(other != EuclideanDomainInt(0, 0, self.d))
        assert(self.d == other.d)

        a = self
        b = other

        q = a // b
        r = a - (q*b)

        return q,r

    def __floordiv__(self, other):

        if isinstance(other, int):
            other = EuclideanDomainInt(other, 0, self.d)

        assert(other != EuclideanDomainInt(0, 0, self.d))
        assert(self.d == other.d)

        a = self
        b = other

        numerator = a * (b.conjugate())
        denominator = b.norm()

        nr = numerator.real
        ni = numerator.imaginary

        qr = nr // denominator
        qi = ni // denominator

        if (2*qr+1)*denominator < 2*nr:
            qr += 1

        if (2*qi+1)*denominator < 2*ni:
            qi += 1

        q = EuclideanDomainInt(qr, qi, self.d)

        return q

    def divmod_brute_force(self, other):
        if isinstance(other, int):
            other = EuclideanDomainInt(other, 0, self.d)

        assert(other != EuclideanDomainInt(0, 0, self.d))
        assert(self.d == other.d)

        a = self
        b = other

        q = a.floor_div_brute_force(b)
        r = a - (q*b)

        return q,r

    def floor_div_brute_force(self, other):
        if isinstance(other, int):
            other = EuclideanDomainInt(other, 0, self.d)

        assert(other != EuclideanDomainInt(0, 0, self.d))
        assert(self.d == other.d)

        a = self
        b = other

        numerator = a * (b.conjugate())
        denominator = b.norm()

        nr = numerator.real
        ni = numerator.imaginary

        qr = nr // denominator
        qi = ni // denominator

        candidates = []
        q = EuclideanDomainInt(qr, qi, self.d)
        r = a - (q*b)
        candidates.append((q,r))

        q = EuclideanDomainInt(qr + 1, qi, self.d)
        r = a - (q*b)
        candidates.append((q,r))

        q = EuclideanDomainInt(qr, qi + 1, self.d)
        r = a - (q*b)
        candidates.append((q,r))

        q = EuclideanDomainInt(qr + 1, qi + 1, self.d)
        r = a - (q*b)
        candidates.append((q,r))

        return EuclideanDomainInt.best_candidate(candidates)[0]

    @staticmethod
    def best_candidate(candidates):
        min_r = float("inf")
        min_c = ()
        for c in candidates:
            c_r = c[1].norm()
            if c_r <= min_r:
                min_r = c_r
                min_c = c
        return min_c

    def __lt__(self, other):
        if isinstance(other, int):
            other = EuclideanDomainInt(other, 0, self.d)

        assert(self.d == other.d)

        return self.norm() < other.norm()

    def __le__(self, other):
        assert(self.d == other.d)
        return self < other or self == other

    def __gt__(self, other):
        if isinstance(other, int):
            other = EuclideanDomainInt(other, 0, self.d)

        assert(self.d == other.d)
        return self.norm() > other.norm()

    def __ge__(self, other):
        assert(self.d == other.d)
        return self > other or self == other

    def __mod__(self, other):
        assert(self.d == other.d)
        q, r = divmod(self, other)
        return r

    def __mul__(self, other):
        # (a+bω)*(m+nω)=(am-bnd)+(an+bm)ω if d % 4 != 1
        # (a+bω)*(m+nω)=(am-bn(d-1)/4)+(an+b(m-n))ω if d % 4 = 1

        if isinstance(other, int):
            other = EuclideanDomainInt(other, 0, self.d)

        assert(self.d == other.d)
        
        a = self.real
        b = self.imaginary
        m = other.real
        n = other.imaginary
        if self.d % 4 != 1:
            product_real = (a*m) + (b*n*self.d)
            product_imaginary = (b*m + a*n)
        else:
            product_real = (a*m) + (b*n)*(self.d - 1)//4
            product_imaginary = (a*n) + b*(m-n)

        return EuclideanDomainInt(product_real, product_imaginary, self.d)

    def __sub__(self, other):
        if isinstance(other, int):
            other = EuclideanDomainInt(other, 0, self.d)

        assert(self.d == other.d)

        a = self.real
        b = self.imaginary
        m = other.real
        n = other.imaginary

        difference_real = a - m
        difference_imaginary = b - n

        return EuclideanDomainInt(difference_real, difference_imaginary, self.d)

    def associates(self):
        units = EuclideanDomainInt.units(self.d)
        associates = list(map(lambda x: x * self, units))
        return associates

    def canonical(self):
        associates = self.associates()

        for a in associates:
            # In first sextant
            if (a.sextant() == EisensteinInt(1)):
                return a

    def complex_form(self):
        r = self.real
        i = self.imaginary

        return r + i*self.omega

    def conjugate(self):

        a = self.real
        b = self.imaginary

        if self.d % 4 != 1:
            conjugate_real = a
        else:
            conjugate_real = a-b
        conjugate_imag = -b
        return EuclideanDomainInt(conjugate_real, conjugate_imag, self.d)

    def dot_product(self, other):
        a = self.real
        b = self.imaginary
        m = other.real
        n = other.imaginary

        dot_product = a*m + b*n - (b*m + a*n) /2
        return dot_product

    def norm(self):
        # Norm(a) = a^2 - db^2 if d % 4 != 1
        # Norm(a) = a^2 - ab - b^2(d - 1)/4 if d % 4 = 1

        a = self.real
        b = self.imaginary
        if self.d % 4 != 1:
            norm = a**2 - b**2*self.d
        else:
            norm = a**2 - a*b - b**2*(self.d - 1)//4
        return norm

    def polar_form(self):
        r = abs(self)

        a = self.real
        b = self.imaginary

        z = a + b*self.omega
        return cmath.phase(z)

    def sextant(self):
        return self.signum()[1]

    def signum(self):
        a = self.real
        b = self.imaginary
        if (a == 0 and b==0):
            return (self, EisensteinInt(0)) # origin
        elif (a > b and b >=0):
            return (self, EisensteinInt(1)) # first
        elif (b >= a and a>0):
            return (self * EisensteinInt(0,-1), EisensteinInt(1,1)) # second
        elif (b > 0 and 0>=a):
            return (EisensteinInt(-1,-1) * self, EisensteinInt(0,1)) # third
        elif (a < b and b<=0):
            return (self * EisensteinInt(-1,0), EisensteinInt(-1)) # fourth
        elif (b <= a and a<0):
            return (self * EisensteinInt(0,1), EisensteinInt(-1,-1)) # fifth
        else:
            return (self * EisensteinInt(1,1), EisensteinInt(0,-1)) # sixth

    @staticmethod
    def quadratic_form(c, d):
        r = c.real
        i = c.imag

        if d % 4 != 1:
            re = r
            ie = i / sqrt(-d)
        else:
            ie = 2* i / sqrt(-d)
            re = r + ((1/2) * ie)

        return EuclideanDomainInt(round(re), round(ie), d)

    @staticmethod
    def units(d):
        units = []
        if d == -1:
            unit = EuclideanDomainInt(1, 0, d)
            for i in range(4):
                units.append(unit)
                unit = unit * EuclideanDomainInt(0, 1, d)
            return units
        elif d == -3:
            unit = EuclideanDomainInt(1,0,d)
            for i in range(4):
                units.append(unit)
                unit = unit * EuclideanDomainInt(1, 1, d)
            return units
        else:
            units = [EuclideanDomainInt(1, 0, d), EuclideanDomainInt(-1, 0, d)]
            return units

    def is_even(self):
        # is_even iff. a+b is congruent to 0 mod 3

        a = self.real
        b = self.imaginary

        if (((a+b) % 3) == 0):
            return True
        else:
            return False

    def is_prime(self):
        # n is an eisenstein prime if
        # 1) b=0 and a=p, prime and p≡2 mod 3
        # 2) a=0 and b=p, prime and p≡2 mod 3
        # 3) Norm(n) = p is prime where p=3 or p≡1 mod 3

        a = self.real
        b = self.imaginary

        if (b == 0 and isprime(a) and (a % 3 == 2)):
            return True
        if (a == 0 and isprime(b) and (b % 3 == 2)):
            return True

        norm = self.norm()
        is_prime = isprime(norm)

        if (is_prime):
            # Norm == 3 is 1-ω
            return (norm == 3) or (norm % 3 == 1)
        else:
            return False

    def is_unit(self):
        if (self.norm() == 1):
            return True
        else:
            return False

    def gcd(self, other):

        a = self
        b = other

        if (a.norm() < b.norm()):
            return b.gcd(a)

        while (b.norm() > 0):
            q, r = divmod(a, b)
            a, b = b, r

        return a
    
class DedekindDomainInt:
    """
    Stores the quadratic integer of Dedekind Doomain in the form a + bω

    Create Dedekind Integer by:
	     a = DedekindDomainInt(5,7)  # Create (5 + 7ω)
	     a = DedekindDomainInt(13)  # Create (13)

    Functions implemented
         Basic functions: init(), ==, hash(),  str(), <, >, <=, >=
	     Arithmetic functions: abs(), +, *, -

         a.debug_str() - Returns a simple debug string to describe the
            data members.

         a.associates() - Returns a list of the product of the number and
            each of the units.
         a.canonical() - Returns the associate in the first sextant.
         a.complex_form() - Returns the complex form.
         a.conjugate() - Returns an EisensteinInt representing he conjugate.
         a.dot_product(b) - Returns the dot product.
         a.norm() - Returns an integer representing the norm.
         a.polar_form() - Returns a tuple of the radius and angle.

         a.gcd(b) - Compute the greatest common divisor of a and b.
    """

    def __init__(self, real=0, imaginary=0, d=-1):
        assert isinstance(real, int)
        assert isinstance(imaginary, int)
        assert isinstance(d, int)

        self.real = real
        self.imaginary = imaginary
        self.d = d
        if d % 4 != 1:
            self.omega = cmath.sqrt(d)
        else:
            self.omega = (-1.0 + cmath.sqrt(d))/2.0

    def __hash__(self):
        return hash((self.real, self.imaginary))

    def __eq__(self, other):
        if isinstance(other, DedekindDomainInt):
            return (self.real == other.real) and (self.imaginary == other.imaginary) and (self.d == other.d)
        else:
            return False

    def __str__(self):
        a = self.real
        b = self.imaginary

        if b == 0:
            return str(a)
        elif (b > 0):
            s = ""
            op = "+"
        else:
            s = "-"
            op = "-"

        if (abs(b) == 1):
            b_ = "ω"
        else:
            b_ = f"{str(abs(b))}ω"

        if a == 0:
            result = f"{s} {b_}"
        else:
            result = f"{str(a)} {op} {b_}"

        return result.strip()
    
    def __repr__(self):
        return f"DedekindDomainInt({self.real}, {self.imaginary}, {self.d})"

    def debug_str(self):
        result = f"DedekindDomeinInt({self.real}, {self.imaginary}, {self.d})"
        return result

    def __abs__(self):
        return sqrt(self.norm())

    def __add__(self, other):
        if isinstance(other, int):
            other = DedekindDomainInt(other)

        a = self.real
        b = self.imaginary
        c = other.real
        d = other.imaginary

        sum_real = a + c
        sum_imaginary = b + d

        return DedekindDomainInt(sum_real, sum_imaginary)

    @staticmethod
    def best_candidate(candidates):
        min_r = float("inf")
        min_c = ()
        for c in candidates:
            c_r = c[1].norm()
            if c_r <= min_r:
                min_r = c_r
                min_c = c
        return min_c

    def __lt__(self, other):
        if isinstance(other, int):
            other = DedekindDomainInt(other, 0, self.d)

        assert(self.d == other.d)

        return self.norm() < other.norm()

    def __le__(self, other):
        assert(self.d == other.d)
        return self < other or self == other

    def __gt__(self, other):
        if isinstance(other, int):
            other = DedekindDomainInt(other, 0, self.d)

        assert(self.d == other.d)
        return self.norm() > other.norm()

    def __ge__(self, other):
        assert(self.d == other.d)
        return self > other or self == other

    def __mod__(self, other):
        assert(self.d == other.d)
        q, r = divmod(self, other)
        return r

    def __mul__(self, other):
        # (a+bω)*(m+nω)=(am-bnd)+(an+bm)ω if d % 4 != 1
        # (a+bω)*(m+nω)=(am-bn(d-1)/4)+(an+b(m-n))ω if d % 4 = 1

        if isinstance(other, int):
            other = DedekindDomainInt(other, 0, self.d)

        assert(self.d == other.d)
        
        a = self.real
        b = self.imaginary
        m = other.real
        n = other.imaginary
        if self.d % 4 != 1:
            product_real = (a*m) + (b*n*self.d)
            product_imaginary = (b*m + a*n)
        else:
            product_real = (a*m) + (b*n)*(self.d - 1)//4
            product_imaginary = (a*n) + b*(m-n)

        return DedekindDomainInt(product_real, product_imaginary, self.d)

    def __sub__(self, other):
        if isinstance(other, int):
            other = DedekindDomainInt(other, 0, self.d)

        assert(self.d == other.d)

        a = self.real
        b = self.imaginary
        m = other.real
        n = other.imaginary

        difference_real = a - m
        difference_imaginary = b - n

        return DedekindDomainInt(difference_real, difference_imaginary, self.d)

    def associates(self):
        units = DedekindDomainInt.units(self.d)
        associates = list(map(lambda x: x * self, units))
        return associates

    def canonical(self):
        associates = self.associates()

        for a in associates:
            # In first sextant
            if (a.sextant() == EisensteinInt(1)):
                return a

    def complex_form(self):
        r = self.real
        i = self.imaginary

        return r + i*self.omega

    def conjugate(self):

        a = self.real
        b = self.imaginary

        if self.d % 4 != 1:
            conjugate_real = a
        else:
            conjugate_real = a-b
        conjugate_imag = -b
        return DedekindDomainInt(conjugate_real, conjugate_imag, self.d)

    def dot_product(self, other):
        a = self.real
        b = self.imaginary
        m = other.real
        n = other.imaginary

        dot_product = a*m + b*n - (b*m + a*n) /2
        return dot_product

    def norm(self):
        # Norm(a) = a^2 - db^2 if d % 4 != 1
        # Norm(a) = a^2 - ab - b^2(d - 1)/4 if d % 4 = 1

        a = self.real
        b = self.imaginary
        if self.d % 4 != 1:
            norm = a**2 - b**2*self.d
        else:
            norm = a**2 - a*b - b**2*(self.d - 1)//4
        return norm

    def polar_form(self):
        r = abs(self)

        a = self.real
        b = self.imaginary

        z = a + b*self.omega
        return cmath.phase(z)

    @staticmethod
    def quadratic_form(c, d):
        r = c.real
        i = c.imag

        if d % 4 != 1:
            re = r
            ie = i / sqrt(-d)
        else:
            ie = 2* i / sqrt(-d)
            re = r + ((1/2) * ie)

        return DedekindDomainInt(round(re), round(ie), d)

    @staticmethod
    def units(d):
        units = []
        if d == -1:
            unit = DedekindDomainInt(1, 0, d)
            for i in range(4):
                units.append(unit)
                unit = unit * DedekindDomainInt(0, 1, d)
            return units
        elif d == -3:
            unit = DedekindDomainInt(1,0,d)
            for i in range(4):
                units.append(unit)
                unit = unit * DedekindDomainInt(1, 1, d)
            return units
        else:
            units = [DedekindDomainInt(1, 0, d), DedekindDomainInt(-1, 0, d)]
            return units

    def is_even(self):
        # is_even iff. a+b is congruent to 0 mod 3

        a = self.real
        b = self.imaginary

        if (((a+b) % 3) == 0):
            return True
        else:
            return False
        
    def factor(self):
        def diop_DN_with_trivial(D, N):
            trivial_solution = []
            sqrtN, issqr = integer_nthroot(N, 2)
            if issqr and sqrtN != 0:
                trivial_solution.append((sqrtN, 0))
                trivial_solution.append((-sqrtN, 0))
            if N % abs(D) == 0:
                sqrtND, issqr = integer_nthroot(N//abs(D), 2)
                if issqr and sqrtND != 0:
                    trivial_solution.append((0, sqrtND))
                    trivial_solution.append((0, -sqrtND))

            if D % 4 != 1:
                solutions = diop_DN(D, N)
            else:
                x, y = symbols("x, y", integer=True)
                D1, N1 = find_DN(x**2 - x*y - (D - 1)//4*y**2 - N)
                A, B = transformation_to_DN(x**2 - x*y - (D - 1)//4*y**2 - N)
                transformed_solutions = diop_DN(D1, N1)
                solutions = []
                if len(transformed_solutions) != 0:
                    for transformed_solution in transformed_solutions:
                        solution = A*Matrix((transformed_solution[0], transformed_solution[1])) + B
                        solutions.append((solution[0], solution[1]))
            
            solutions_2nd_quadrant = []
            solutions_3rd_quadrant = []
            solutions_4th_quadrant = []
            for solution in solutions:
                solutions_2nd_quadrant.append((-solution[0], solution[1]))
                solutions_3rd_quadrant.append((-solution[0], -solution[1]))
                solutions_4th_quadrant.append((solution[0], -solution[1]))
            return trivial_solution + solutions + solutions_2nd_quadrant + solutions_3rd_quadrant + solutions_4th_quadrant
        
        factors = []
        n = self.norm()
        pf = factorint(n, multiple=True)
        f = set([math.prod(x) for x in list(chain.from_iterable(combinations(pf, r) for r in range(0, len(pf)+1)))])
        for x in sorted(list(f))[:math.ceil(len(f)/2)]:
            y = n//x
            x_solutions = diop_DN_with_trivial(self.d, x)
            if len(x_solutions) > 0:
                y_solutions = diop_DN_with_trivial(self.d, y)
                if len(y_solutions) > 0:
                    for x_solution in x_solutions:
                        for y_solution in y_solutions:
                            cand_x = DedekindDomainInt(int(x_solution[0]), int(x_solution[1]), self.d)
                            cand_y = DedekindDomainInt(int(y_solution[0]), int(y_solution[1]), self.d)
                            cand_mul = cand_x*cand_y
                            if cand_mul == self:
                                factors.append(cand_x)
                                factors.append(cand_y)
        return factors

    def is_prime(self):
        # n is an eisenstein prime if
        # 1) b=0 and a=p, prime and p≡2 mod 3
        # 2) a=0 and b=p, prime and p≡2 mod 3
        # 3) Norm(n) = p is prime where p=3 or p≡1 mod 3
        
        if self.d == -3:
            a = self.real
            b = self.imaginary

            if (b == 0 and isprime(a) and (a % 3 == 2)):
                return True
            if (a == 0 and isprime(b) and (b % 3 == 2)):
                return True

            norm = self.norm()
            is_prime = isprime(norm)

            if (is_prime):
                # Norm == 3 is 1-ω
                return (norm == 3) or (norm % 3 == 1)
            else:
                return False
        else:
            factors = self.factor()
            norms = set([x.norm() for x in factors])
            if len(norms) == 1 and norms[0] == 1:
                return True
            else:
                return False

    def is_unit(self):
        if (self.norm() == 1):
            return True
        else:
            return False

    def gcd(self, other):
        factor_self = set(self.factor())
        factor_other = set(other.factor())
        factor_cd = factor_self & factor_other
        return max(factor_cd)        

if __name__ == "__main__":
    z1 = EuclideanDomainInt(1, 0, -1)
    z2 = EuclideanDomainInt(0, 1, -1)

    z = z1.gcd(z2)
    print(z)

    z1 = DedekindDomainInt(6, 0, -5)
    z2 = DedekindDomainInt(1, -1, -5)

    z = z1.gcd(z2)
    print(z.norm())