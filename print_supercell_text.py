import cmath, math
from quadraticinteger import EisensteinInt_hex, GaussInt_square

def wood_notation_hex_eigen(real, complex):
    omega = (-1 + cmath.sqrt(-3.0))/2.0
    eisen_number = real + complex*omega
    norm = abs(eisen_number)**2
    angle = cmath.phase(eisen_number)
    wood_str = f'sqrt({norm:.1f})×sqrt({norm:.1f})R({angle/cmath.pi*180:.1f} degree)'
    return wood_str

def wood_notation_quad_rec_latex_style(a, b, D):
    quad_int = a + b*cmath.sqrt(D)
    norm = (quad_int*(quad_int.conjugate())).real
    angle = cmath.phase(quad_int)
    wood_str = f'(\\sqrt{{{norm:.0f}}}\\times\\sqrt{{{norm:.0f}}})\\mathrm{{R}}{angle/cmath.pi*180:.1f}^{{\\circ}}'
    return wood_str

def wood_notation_hex_eisen_latex_style(real, complex):
    int_tol = 1e-5
    omega = (-1 + cmath.sqrt(-3.0))/2.0
    eisen_number = real + complex*omega
    abs_eisen_number = abs(eisen_number)
    norm = abs_eisen_number**2
    angle = cmath.phase(eisen_number)
    if abs(abs_eisen_number - math.floor(abs_eisen_number)) > int_tol:
        wood_str = f'(\\sqrt{{{norm:.0f}}}\\times\\sqrt{{{norm:.0f}}})\\mathrm{{R}}{angle/cmath.pi*180:.1f}^{{\\circ}}'
    else:
        wood_str = f'{abs_eisen_number:.0f}\\times{abs_eisen_number:.0f}\\mathrm{{R}}{angle/cmath.pi*180:.1f}^{{\\circ}}'
    return wood_str

def wood_notation_square_gauss_latex_style(real, complex):
    int_tol = 1e-5
    omega = 1j
    eisen_number = real + complex*omega
    abs_eisen_number = abs(eisen_number)
    norm = abs_eisen_number**2
    angle = cmath.phase(eisen_number)
    if abs(abs_eisen_number - math.floor(abs_eisen_number)) > int_tol:
        wood_str = f'(\\sqrt{{{norm:.0f}}}\\times\\sqrt{{{norm:.0f}}})\\mathrm{{R}}{angle/cmath.pi*180:.1f}^{{\\circ}}'
    else:
        wood_str = f'{abs_eisen_number:.0f}\\times{abs_eisen_number:.0f}\\mathrm{{R}}{angle/cmath.pi*180:.1f}^{{\\circ}}'
    return wood_str

def matrix_notation_square_gauss_latex_style(real, complex):
    matrix_str = f'\\begin{{pmatrix}} {real} & {complex} \\\\ {-complex} & {real}\\end{{pmatrix}}'
    return matrix_str

def matrix_notation_hex_eisen_latex_style(real, complex):
    matrix_str = f'\\begin{{pmatrix}} {real} & {complex} \\\\ {-complex} & {real - complex}\\end{{pmatrix}}'
    return matrix_str

def matrix_notation_quad_rec_latex_style(a, b, D):
    matrix_str = f'\\begin{{pmatrix}} {a} & {b} \\\\ {b*D} & {a}\\end{{pmatrix}}'
    return matrix_str

def common_supercell_wood_notations(eisenint_supercells, error, D):
    common_supercell_str = ""
    if len(eisenint_supercells) == 2:
        if D == -3:
            eiseninttype_supercells = [EisensteinInt_hex(eisenint.real, eisenint.imaginary) for eisenint in eisenint_supercells]
            wood_notations = [wood_notation_hex_eisen_latex_style(eisenint.real, eisenint.imaginary) for eisenint in eiseninttype_supercells]
            angle_diff = cmath.phase(eiseninttype_supercells[0].complex_form()) - cmath.phase(eiseninttype_supercells[1].complex_form())
            common_supercell_str = f'{angle_diff/cmath.pi*180:.2f} & ${wood_notations[0]}$ & $ {matrix_notation_hex_eisen_latex_style(eisenint_supercells[0].real, eisenint_supercells[0].imaginary)}\\cong {eiseninttype_supercells[0]}$ & ${wood_notations[1]}$ & $ {matrix_notation_hex_eisen_latex_style(eisenint_supercells[1].real, eisenint_supercells[1].imaginary)} \\cong {eiseninttype_supercells[1]}$ & {eiseninttype_supercells[0].norm() + eiseninttype_supercells[1].norm()} & {error*100:.3f}\\% \\\\'
        elif D == -1:
            eiseninttype_supercells = [GaussInt_square(eisenint.i, eisenint.r) for eisenint in eisenint_supercells]
            wood_notations = [wood_notation_square_gauss_latex_style(gaussint.r, gaussint.i) for gaussint in eisenint_supercells]
            angle_diff = eisenint_supercells[0].arg() - eisenint_supercells[1].arg()
            common_supercell_str = f'{angle_diff/cmath.pi*180:.2f} & ${wood_notations[0]}$ & $ {matrix_notation_square_gauss_latex_style(eisenint_supercells[0].r, eisenint_supercells[0].i)}\\cong {eisenint_supercells[0]}$ & ${wood_notations[1]}$ & $ {matrix_notation_square_gauss_latex_style(eisenint_supercells[1].r, eisenint_supercells[1].i)} \\cong {eisenint_supercells[1]}$ & {eisenint_supercells[0].norm() + eisenint_supercells[1].norm()} & {error*100:.3f}\\% \\\\'
        else:
            quadints = [[int(x.real), int(x.imag/math.sqrt(-D))] for x in eisenint_supercells]
            angle_diff = cmath.phase(eisenint_supercells[0]) - cmath.phase(eisenint_supercells[1])
            wood_notations = [wood_notation_quad_rec_latex_style(quadint[0], quadint[1], D) for quadint in quadints]
            common_supercell_str = f'{angle_diff/cmath.pi*180:.2f} & ${wood_notations[0]}$ & $ {matrix_notation_quad_rec_latex_style(quadints[0][0], quadints[0][1], D)}\\cong {quadints[0][0]} + {quadints[0][1]}\\omega_{{rec}}$ & ${wood_notations[1]}$ & $ {matrix_notation_quad_rec_latex_style(quadints[1][0], quadints[1][1], D)} \\cong {quadints[1][0]} + {quadints[1][1]}\\omega_{{rec}}$ & {round(abs(eisenint_supercells[0])**2 + abs(eisenint_supercells[1])**2):.0f} & {error*100:.3f}\\% \\\\'

    return common_supercell_str

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 3:
        real_part = int(sys.argv[1])
        imag_part = int(sys.argv[2])
    else:
        real_part = 1
        imag_part = 1
    wood_notation = wood_notation_hex_eisen_latex_style(real_part, imag_part)
    print(wood_notation)

    common_str = common_supercell_wood_notations([[1, -15], [-1, -16]], 0, -3)
    print(common_str)