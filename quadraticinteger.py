from eisenstein.eisenstein import EisensteinInt
from NumberTheoryPython.gaussint import GaussInt
import cmath

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