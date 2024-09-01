import cmath, math
from quadraticinteger import EisensteinInt_hex, GaussInt_square, DedekindDomainInt

def wood_notation_quad_latex_style(a, b, D):
    int_tol = 1e-5
    quad_int = DedekindDomainInt(a, b, D)
    norm = quad_int.norm()
    sqrt_norm = math.sqrt(norm)
    angle = cmath.phase(quad_int.complex_form())
    if abs(sqrt_norm - math.floor(sqrt_norm)) > int_tol:
        wood_str = f'(\\sqrt{{{norm:.0f}}}\\times\\sqrt{{{norm:.0f}}})\\mathrm{{R}}{angle/cmath.pi*180:.1f}^{{\\circ}}'
    else:
        wood_str = f'({sqrt_norm:.0f}\\times{sqrt_norm:.0f})\\mathrm{{R}}{angle/cmath.pi*180:.1f}^{{\\circ}}'
    return wood_str

def matrix_notation_quad_latex_style(a, b, quad_d):
    if quad_d % 4 != 1:
        matrix_str = f'\\begin{{pmatrix}} {a} & {b} \\\\ {b*quad_d} & {a}\\end{{pmatrix}}'
    elif quad_d % 4 == 1:
        matrix_str = f'\\begin{{pmatrix}} {a} & {b} \\\\ {b*(quad_d - 1)//4} & {a - b}\\end{{pmatrix}}'
    return matrix_str

def common_supercell_wood_notations(eisenint_supercells, error, D):
    common_supercell_str = ""
    if len(eisenint_supercells) == 3:
        if D == -3:
            sym_id = 'hex'
        elif D == -1:
            sym_id = 'sq'
        elif D % 4 != 1:
            sym_id = 'rec'
        elif D % 4 == 1:
            sym_id = 'obl'
        quadints = eisenint_supercells[0:2]
        angle_diff = eisenint_supercells[2]
        wood_notations = [wood_notation_quad_latex_style(quadint.real, quadint.imaginary, D) for quadint in quadints]
        common_supercell_str = f'{angle_diff:.2f} & ${wood_notations[0]}$ & $ {matrix_notation_quad_latex_style(quadints[0].real, quadints[0].imaginary, D)}\\cong {quadints[0].real} + {quadints[0].imaginary}\\omega_{{{sym_id}}}$ & ${wood_notations[1]}$ & $ {matrix_notation_quad_latex_style(quadints[1].real, quadints[1].imaginary, D)} \\cong {quadints[1].real} + {quadints[1].imaginary}\\omega_{{{sym_id}}}$ & {quadints[0].norm() + quadints[1].norm()} & {error*100:.3f}\\% \\\\'

    return common_supercell_str

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 3:
        real_part = int(sys.argv[1])
        imag_part = int(sys.argv[2])
    else:
        real_part = 1
        imag_part = 1
    wood_notation = wood_notation_quad_latex_style(real_part, imag_part, -3)
    print(wood_notation)

    common_str = common_supercell_wood_notations([[1, -15], [-1, -16]], 0, -3)
    print(common_str)