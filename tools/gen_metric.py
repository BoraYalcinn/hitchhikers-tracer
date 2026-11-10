import argparse
import re
import sys
from pathlib import Path

import sympy as sp


# ---------------------------------------------------------------------------
# Spacetimes:
# ---------------------------------------------------------------------------

def schwarzschild():
    """Schwarzschild metric, signature (-,+,+,+), M = 1."""
    t, r, theta, phi = sp.symbols('t r theta phi', real=True)
    coords = [t, r, theta, phi]

    f = 1 - 2 / r
    g = sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(theta)**2)
    return coords, g


# Add new spacetimes here, then register them below.
METRICS = {
    'schwarzschild': schwarzschild,
}


# ---------------------------------------------------------------------------
# Processing: 
# ---------------------------------------------------------------------------

def inverse_and_derivatives(coords, g):
    """Return the inverse metric and its derivatives.

    ginv  : 4x4 matrix, g^{mu nu}
    dginv : list of 4 matrices, dginv[a][mu, nu] = d(g^{mu nu}) / d(x^a)
    """
    ginv = g.inv()
    ginv = ginv.applyfunc(sp.simplify)
    dginv = [ginv.diff(c).applyfunc(sp.simplify) for c in coords]
    return ginv, dginv


def geodesic_rhs_entries(ginv, dginv):
    """Right-hand side of the geodesic equation in Hamiltonian form.

    Returns the momentum symbols and the list of (Fortran target, expression):
      dx^mu/dlambda =  sum_nu g^{mu nu} p_nu
      dp_mu/dlambda = -1/2 sum_{a,b} (d_mu g^{ab}) p_a p_b
    Expanding the sums symbolically drops every term that is identically
    zero, so the generated code only evaluates the non-zero ones.
    """
    p = sp.symbols('p1 p2 p3 p4', real=True)
    entries = []
    for mu in range(4):
        e = sum(ginv[mu, nu] * p[nu] for nu in range(4))
        entries.append((f'dx({mu + 1})', sp.simplify(e)))
    for mu in range(4):
        e = -sp.Rational(1, 2) * sum(dginv[mu][a, b] * p[a] * p[b]
                                     for a in range(4) for b in range(4))
        entries.append((f'dp({mu + 1})', sp.simplify(e)))
    return p, entries


def print_report(coords, ginv, dginv):
    """Print the inverse metric and the non-zero derivatives."""
    print('Inverse metric:')
    sp.pprint(ginv)
    print()
    print('Non-zero derivatives of the inverse metric:')
    for a in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                d = dginv[a][mu, nu]
                if d != 0:
                    print(f'd/d{coords[a]} g^({mu},{nu}) =', d)


# ---------------------------------------------------------------------------
# Fortran output
# ---------------------------------------------------------------------------

def to_fortran(expr, kind='wp'):
    """Convert a sympy expression to a Fortran expression string.

    kind : name of the Fortran kind parameter used for real literals
    """
    text = sp.fcode(expr, source_format='free', standard=95)
    # sympy writes literals like 2.0d0; use the kind parameter instead
    text = re.sub(r'(\d+\.\d*)d0', rf'\1_{kind}', text)
    # whole-number literals such as 1d0 become 1.0_wp
    return re.sub(r'\b(\d+)d0\b', rf'\1.0_{kind}', text)


def symmetric_entries(name, m, lead=''):
    """Non-zero entries of a symmetric 4x4 matrix as (Fortran target, expression).

    name : Fortran array name, e.g. 'ginv'
    lead : text put before the two indices, e.g. '2,' gives name(2,mu,nu)
    """
    entries = []
    for mu in range(4):
        for nu in range(mu, 4):
            if m[mu, nu] != 0:
                entries.append((f'{name}({lead}{mu + 1},{nu + 1})', m[mu, nu]))
                if mu != nu:
                    entries.append((f'{name}({lead}{nu + 1},{mu + 1})', m[mu, nu]))
    return entries


def subroutine_lines(sub_name, outputs, coords, entries, momenta=None, zero_outputs=True,
                     kind='wp'):
    """Fortran lines of one subroutine that fills its output arrays.

    sub_name     : name of the subroutine, e.g. 'metric_inv'
    outputs      : list of (array name, shape text), e.g. [('ginv', '(4,4)')]
    entries      : list of (target text, expression) to assign
    momenta      : symbols p1..p4 if the routine also takes p(4) as input
    zero_outputs : zero the output arrays first (not needed when every
                   element is assigned explicitly)
    kind         : Fortran kind parameter of every real in the routine
                   ('wp' double precision, 'sp' single precision)
    """
    # common subexpression elimination: shared parts are computed once
    exprs = [e for _, e in entries]
    temps, reduced = sp.cse(exprs, symbols=sp.numbered_symbols('c'))

    # which coordinates are really used?
    used = set()
    for _, e in temps:
        used |= e.free_symbols
    for e in reduced:
        used |= e.free_symbols
    used_coords = [(i, c) for i, c in enumerate(coords) if c in used]
    used_momenta = [(i, m) for i, m in enumerate(momenta or ()) if m in used]

    arg_names = ['x'] + (['p'] if momenta else []) + [var for var, _ in outputs]
    lines = []
    lines.append(f'  attributes(host,device) pure subroutine {sub_name}({", ".join(arg_names)})')
    lines.append(f'    real({kind}), intent(in)  :: x(4)')
    if momenta:
        lines.append(f'    real({kind}), intent(in)  :: p(4)')
    for var, shape in outputs:
        lines.append(f'    real({kind}), intent(out) :: {var}{shape}')
    names = ([str(c) for _, c in used_coords] + [str(m) for _, m in used_momenta]
             + [str(s) for s, _ in temps])
    if names:
        lines.append(f'    real({kind}) :: ' + ', '.join(names))
    lines.append('')
    for i, c in used_coords:
        lines.append(f'    {c} = x({i + 1})')
    for i, m in used_momenta:
        lines.append(f'    {m} = p({i + 1})')
    if used_coords or used_momenta:
        lines.append('')
    if zero_outputs:
        for var, _ in outputs:
            lines.append(f'    {var} = 0.0_{kind}')
        lines.append('')
    for sym, e in temps:
        lines.append(f'    {sym} = {to_fortran(e, kind)}')
    if temps:
        lines.append('')
    for (target, _), e in zip(entries, reduced):
        lines.append(f'    {target} = {to_fortran(e, kind)}')
    lines.append(f'  end subroutine {sub_name}')
    return lines


def generate_fortran(name, coords, g, ginv, dginv):
    """Build the text of the Fortran module for one spacetime."""
    # inverse metric and its derivatives
    entries_inv = symmetric_entries('ginv', ginv)
    for a in range(4):
        entries_inv += symmetric_entries('dginv', dginv[a], lead=f'{a + 1},')
    # metric itself (needed to lower indices)
    entries_cov = symmetric_entries('gcov', g)
    # geodesic equation with the zero terms already dropped
    momenta, entries_rhs = geodesic_rhs_entries(ginv, dginv)

    lines = []
    lines.append('! AUTO-GENERATED by tools/gen_metric.py -- DO NOT EDIT.')
    lines.append(f'! Metric: {name}. Regenerate with: python tools/gen_metric.py {name}')
    lines.append('!')
    lines.append('! x(1:4)          coordinates (t, r, theta, phi)')
    lines.append('! gcov(mu,nu)     metric g_{mu nu}')
    lines.append('! ginv(mu,nu)     inverse metric g^{mu nu}')
    lines.append('! dginv(a,mu,nu)  derivative d g^{mu nu} / d x^a')
    lines.append('! metric_rhs      geodesic equation dx/dlambda, dp/dlambda for momentum p,')
    lines.append('!                 with only the non-zero metric terms evaluated')
    lines.append('!')
    lines.append('! The *_sp routines are the same formulas in single precision for the')
    lines.append('! GPU renderer; consumer GPUs run double precision at 1/64 of the')
    lines.append('! single precision rate. The double precision routines stay the')
    lines.append('! reference that the tests check against.')
    lines.append('module metric')
    lines.append('  implicit none')
    lines.append('  private')
    lines.append('  public :: wp, sp, metric_inv, metric_cov, metric_rhs')
    lines.append('  public :: metric_cov_sp, metric_rhs_sp')
    lines.append('')
    lines.append('  integer, parameter :: wp = kind(1.0d0)')
    lines.append('  integer, parameter :: sp = kind(1.0e0)')
    lines.append('')
    lines.append('contains')
    lines.append('')
    lines += subroutine_lines('metric_inv', [('ginv', '(4,4)'), ('dginv', '(4,4,4)')],
                              coords, entries_inv)
    lines.append('')
    lines += subroutine_lines('metric_cov', [('gcov', '(4,4)')], coords, entries_cov)
    lines.append('')
    lines += subroutine_lines('metric_rhs', [('dx', '(4)'), ('dp', '(4)')], coords,
                              entries_rhs, momenta=momenta, zero_outputs=False)
    lines.append('')
    lines += subroutine_lines('metric_cov_sp', [('gcov', '(4,4)')], coords, entries_cov,
                              kind='sp')
    lines.append('')
    lines += subroutine_lines('metric_rhs_sp', [('dx', '(4)'), ('dp', '(4)')], coords,
                              entries_rhs, momenta=momenta, zero_outputs=False, kind='sp')
    lines.append('')
    lines.append('end module metric')
    return '\n'.join(lines) + '\n'


# ---------------------------------------------------------------------------
# Entry point: 
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description='Generate the Fortran metric module.')
    parser.add_argument('name', nargs='?', default='schwarzschild',
                        help=f'spacetime, one of: {", ".join(METRICS)}')
    parser.add_argument('--report', action='store_true',
                        help='print the inverse metric and derivatives')
    parser.add_argument('--out', default=None,
                        help='output file (default: src/generated/metric.f90)')
    args = parser.parse_args()

    if args.name not in METRICS:
        print(f'Unknown metric "{args.name}". Available: {", ".join(METRICS)}')
        sys.exit(1)

    coords, g = METRICS[args.name]()
    assert g.shape == (4, 4)
    ginv, dginv = inverse_and_derivatives(coords, g)

    if args.report:
        print_report(coords, ginv, dginv)

    if args.out:
        out = Path(args.out)
    else:
        out = Path(__file__).resolve().parent.parent / 'src' / 'generated' / 'metric.f90'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(generate_fortran(args.name, coords, g, ginv, dginv))
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()