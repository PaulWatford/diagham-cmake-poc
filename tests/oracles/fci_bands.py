#!/usr/bin/env python3
"""Band structures of the FCI lattice models from their published Bloch Hamiltonians (no DiagHam).

    python3 tests/oracles/fci_bands.py --write tests/data/fci_bands
    python3 tests/oracles/fci_bands.py --check tests/data/fci_bands      (selftest.fci_bands_oracle)

Further models, Bloch Hamiltonians as published, momenta k = 2 pi (i/Nx, j/Ny[, l/Nz]), eigenvalues of the
dense matrix sorted per momentum (numpy, nothing from DiagHam):
    haldane      Haldane (PRL 61, 2015): B1 = t1 (1 + e^{i(x+y)} + e^{iy}), d0 = 2 t2 cos(phi) (cos x + cos y + cos(x+y)),
                 d3 = 2 t2 sin(phi) (sin x + sin y - sin(x+y)) + mu_s, phi in radians
    kagome       three-band Kagome, nearest and next-nearest hopping with spin-orbit terms (Tang, Mei & Wen, PRL 106, 236802)
    zhangqi      Zhang-Qi two-orbital square lattice (Wang, Lu & Zhang, PRL 107, 146803)
    bhz          Bernevig-Hughes-Zhang, four bands with an inversion-breaking Delta
    simple_ti3d  3D simple topological insulator, d = (sin kx, sin ky, sin kz, M - cos kx - cos ky - cos kz)

Checkerboard lattice (Sun, Gu, Katsura & Das Sarma, PRL 106, 236803 (2011)), two bands:
    h(k) = d1(k) 1 + Re B1(k) sigma_x + Im B1(k) sigma_y + d3(k) sigma_z
    B1 = 4 t1 [cos(kx/2) cos(ky/2) cos(phi) + i sin(kx/2) sin(ky/2) sin(phi)],  phi = pi/4
    d1 = 4 tpp cos kx cos ky,   d3 = mu_s + 2 t2 (cos kx - cos ky)
    E_{0,1}(k) = d1 -/+ sqrt(|B1|^2 + d3^2),  k = 2 pi (i/Nx, j/Ny).
The program's default hoppings t1 = 1, t2 = 1/(2+sqrt 2), tpp = 1/(2+2 sqrt 2) are the paper's
flat-band point. The expected file lists kx ky E_0 E_1 per momentum, in the order the program's
--export-onebodytext writes them; the tests compare columns E_0 and E_1 (check_spectrum spectrum).
"""
import argparse, itertools, math, sys
from pathlib import Path

import numpy as np


def checkerboard(nx, ny, t1=1.0, t2=1.0 / (2 + math.sqrt(2)), tpp=1.0 / (2 + 2 * math.sqrt(2)), mu_s=0.0):
    rows = []
    for i in range(nx):
        for j in range(ny):
            kx, ky = 2 * math.pi * i / nx, 2 * math.pi * j / ny
            b1 = complex(4 * t1 * math.cos(kx / 2) * math.cos(ky / 2) * math.cos(math.pi / 4),
                         4 * t1 * math.sin(kx / 2) * math.sin(ky / 2) * math.sin(math.pi / 4))
            d1 = 4 * tpp * math.cos(kx) * math.cos(ky)
            d3 = mu_s + 2 * t2 * (math.cos(kx) - math.cos(ky))
            g = math.sqrt(abs(b1) ** 2 + d3 * d3)
            rows.append((i, j, d1 - g, d1 + g))
    return rows


def _bands(build, dims):
    rows = []
    for idx in itertools.product(*[range(n) for n in dims]):
        k = [2 * math.pi * i / n for i, n in zip(idx, dims)]
        rows.append((*idx, *[float(x) for x in np.linalg.eigvalsh(build(*k))]))
    return rows


def haldane(x, y, t1=1.0, t2=1.0, phi=1.0 / 3, mu_s=0.0):
    b1 = t1 * complex(1 + math.cos(x + y) + math.cos(y), math.sin(x + y) + math.sin(y))
    d0 = 2 * t2 * math.cos(phi) * (math.cos(x) + math.cos(y) + math.cos(x + y))
    d3 = 2 * t2 * math.sin(phi) * (math.sin(x) + math.sin(y) - math.sin(x + y)) + mu_s
    return np.array([[d0 + d3, b1], [np.conj(b1), d0 - d3]])


def kagome(kx, ky, t1=1.0, t2=-0.3, l1=0.28, l2=0.2):
    hab = complex(-2 * t1, -2 * l1) * math.cos(kx / 2) + complex(-2 * t2, 2 * l2) * math.cos((kx - 2 * ky) / 2)
    hac = complex(-2 * t1, 2 * l1) * math.cos(ky / 2) + complex(-2 * t2, -2 * l2) * math.cos((2 * kx - ky) / 2)
    hbc = complex(-2 * t1, -2 * l1) * math.cos((kx - ky) / 2) + complex(-2 * t2, 2 * l2) * math.cos((kx + ky) / 2)
    h = np.zeros((3, 3), dtype=complex)
    h[0, 1], h[0, 2], h[1, 2] = hab, hac, hbc
    return h + np.triu(h, 1).conj().T


def zhangqi(kx, ky, theta=math.pi / 4, mu_s=0.0):
    haa = 2 * math.sqrt(2) * (math.sin(theta) * math.cos(kx) + math.cos(theta) * math.cos(ky))
    hab = math.sqrt(2) * (math.cos(theta) * complex(math.cos(kx), math.sin(kx)) - math.sin(theta) * complex(math.cos(ky), math.sin(ky)))
    return np.array([[mu_s + haa, hab], [np.conj(hab), -haa]])


def bhz4(kx, ky, a=-13.68, b=-16.9, c=-0.0263, d=-0.514, m=-2.058, delta=1.2):
    hd = c - 2 * d * (2 - math.cos(kx) - math.cos(ky))
    h3 = -2 * b * (2 - math.cos(kx) - math.cos(ky)) + m
    h = np.zeros((4, 4), dtype=complex)
    h[0, 0], h[1, 1], h[2, 2], h[3, 3] = hd + h3, hd - h3, hd + h3, hd - h3
    h[0, 1], h[0, 3], h[1, 2], h[2, 3] = complex(a * math.sin(kx), a * math.sin(ky)), delta, -delta, complex(-a * math.sin(kx), a * math.sin(ky))
    return h + np.triu(h, 1).conj().T


def simple_ti3d(kx, ky, kz, mass=1.5):
    d2 = complex(math.sin(ky), -math.sin(kz))
    d1 = math.sin(kx)
    d3 = mass - math.cos(kx) - math.cos(ky) - math.cos(kz)
    h = np.zeros((4, 4), dtype=complex)
    h[0, 0], h[1, 1], h[2, 2], h[3, 3] = d3, -d3, d3, -d3
    h[0, 1], h[2, 3], h[0, 3], h[1, 2] = d1, d1, d2, d2
    return h + np.triu(h, 1).conj().T


CASES = [
    ("checkerboard_3x3", lambda: checkerboard(3, 3), "default hoppings (flat-band point), 3x3 unit cells"),
    ("checkerboard_4x4", lambda: checkerboard(4, 4), "default hoppings, 4x4"),
    ("checkerboard_4x3_t2_0.2_tpp_0.1", lambda: checkerboard(4, 3, t2=0.2, tpp=0.1), "away from the flat-band point"),
    ("haldane_3x3", lambda: _bands(haldane, (3, 3)), "Haldane, t1 = t2 = 1, phi = 1/3 rad, 3x3"),
    ("haldane_4x3_t2_0.4_phi_0.9_mus_0.3", lambda: _bands(lambda x, y: haldane(x, y, 1.0, 0.4, 0.9, 0.3), (4, 3)), "Haldane t2 = 0.4, phi = 0.9, mu_s = 0.3, 4x3"),
    ("kagome_3x3", lambda: _bands(kagome, (3, 3)), "Kagome three bands, defaults, 3x3"),
    ("kagome_4x3_t2_0.1_l1_0.5", lambda: _bands(lambda x, y: kagome(x, y, 1.0, 0.1, 0.5, 0.0), (4, 3)), "Kagome t2 = 0.1, l1 = 0.5, l2 = 0, 4x3"),
    ("zhangqi_3x3", lambda: _bands(zhangqi, (3, 3)), "Zhang-Qi, theta = pi/4, 3x3"),
    ("zhangqi_4x3_theta_0.15_mus_0.2", lambda: _bands(lambda x, y: zhangqi(x, y, 0.15 * math.pi, 0.2), (4, 3)), "Zhang-Qi theta = 0.15 pi, mu_s = 0.2, 4x3"),
    ("bhz_3x3", lambda: _bands(bhz4, (3, 3)), "BHZ four bands, defaults, 3x3"),
    ("simple_ti3d_3x3x2_m1.5", lambda: _bands(simple_ti3d, (3, 3, 2)), "3D simple TI, M = 1.5, 3x3x2"),
    ("simple_ti3d_3x2x3_m2.5", lambda: _bands(lambda x, y, z: simple_ti3d(x, y, z, 2.5), (3, 2, 3)), "3D simple TI, M = 2.5, 3x2x3"),
]


def fukui_hatsugai(build, n, band=0):
    # Lattice Chern number of one band (Fukui, Hatsugai & Suzuki 2005): the sum of the gauge-invariant
    # plaquette Berry phases arg(<u1|u2><u2|u3><u3|u4><u4|u1>) / 2 pi, links along +kx then +ky. It needs a
    # Bloch matrix periodic under k -> k + 2 pi, which the Haldane matrix is.
    u = {}
    for i in range(n):
        for j in range(n):
            u[i, j] = np.linalg.eigh(build(2 * math.pi * i / n, 2 * math.pi * j / n))[1][:, band]
    tot = 0.0
    for i in range(n):
        for j in range(n):
            a, b, c, d = u[i, j], u[(i + 1) % n, j], u[(i + 1) % n, (j + 1) % n], u[i, (j + 1) % n]
            tot += float(np.angle(np.vdot(a, b) * np.vdot(b, c) * np.vdot(c, d) * np.vdot(d, a)))
    return round(tot / (2 * math.pi))


# DiagHam's --singleparticle-chernnumber has the opposite orientation to the link order above (its value
# is minus the Fukui-Hatsugai number for Haldane) and is a small-angle approximation that reaches the
# integer only as 1/N^2, so the tests run it on a 64x64 grid with tolerance 0.02.
def _haldane_chern(phi, mu):
    return -fukui_hatsugai(lambda x, y: haldane(x, y, 1.0, 1.0, phi, mu), 24)


CHERN_CASES = [
    ("chern_haldane_phi0.5_mus0.3", lambda: _haldane_chern(0.5, 0.3), "Haldane phi = 0.5 rad, mu_s = 0.3, lowest band"),
    ("chern_haldane_phi-0.5_mus0.3", lambda: _haldane_chern(-0.5, 0.3), "Haldane phi = -0.5 rad, mu_s = 0.3"),
    ("chern_haldane_phi0.5_mus6", lambda: _haldane_chern(0.5, 6.0), "Haldane phi = 0.5 rad, mu_s = 6, trivial insulator"),
]


def render_all():
    return {f"{name}.dat": f"# momentum indices then sorted energies -- {note}; from the Bloch Hamiltonian in tests/oracles/fci_bands.py, do not edit\n"
            + "".join(" ".join(str(v) if isinstance(v, int) else f"{v:.15g}" for v in r) + "\n" for r in f()) for name, f, note in CASES} | {
        f"{name}.txt": f"# lattice Chern number, {note}; from tests/oracles/fci_bands.py, do not edit\n{f()}\n" for name, f, note in CHERN_CASES}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", metavar="DIR")
    ap.add_argument("--check", metavar="DIR")
    a = ap.parse_args()
    files = render_all()
    if a.write:
        Path(a.write).mkdir(parents=True, exist_ok=True)
        for n, t in files.items():
            (Path(a.write) / n).write_text(t)
        print(f"wrote {len(files)} files to {a.write}")
    if a.check:
        bad = 0
        for n, t in files.items():
            p = Path(a.check) / n
            if not p.exists():
                bad += 1; print(f"MISSING {p}"); continue
            av = [[float(x) for x in l.split()] for l in p.read_text().splitlines() if l and not l.startswith("#")]
            bv = [[float(x) for x in l.split()] for l in t.splitlines() if l and not l.startswith("#")]
            if len(av) != len(bv) or max(abs(x - y) for ra, rb in zip(av, bv) for x, y in zip(ra, rb)) > 1e-12:
                bad += 1; print(f"MISMATCH {p}")
        print(f"{len(files)} files, {bad} mismatches")
        sys.exit(1 if bad else 0)
    if not a.write and not a.check:
        for n, t in files.items():
            print(f"== {n}\n" + "\n".join(t.splitlines()[:4]))


if __name__ == "__main__":
    main()
