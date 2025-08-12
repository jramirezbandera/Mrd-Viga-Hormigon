import math
from dataclasses import dataclass
from typing import List, Tuple, Literal


@dataclass
class SteelLayer:
    y_mm: float   # cota desde el borde superior [mm]
    As_mm2: float  # área de la capa [mm²]


@dataclass
class Materials:
    fck_MPa: float = 30.0
    fyk_MPa: float = 500.0
    gamma_c: float = 1.5
    gamma_s: float = 1.15
    alpha_cc: float = 0.85
    Es_MPa: float = 200_000.0
    # Diagrama parabola-rectangulo EC-2 (fck <= 50 MPa):
    eps_c2: float = 0.0020
    eps_cu2: float = 0.0035

    @property
    def fcd(self) -> float:
        return self.alpha_cc * self.fck_MPa / self.gamma_c

    @property
    def fyd(self) -> float:
        return self.fyk_MPa / self.gamma_s


def area_barras(diams_mm: List[float]) -> float:
    return sum(math.pi * (d**2) / 4.0 for d in diams_mm)


def sigma_c_ec2(eps_c: float, fcd: float, eps_c2: float, eps_cu2: float) -> float:
    """
    Tensión de compresión del hormigón (MPa) para εc >= 0.
    EC-2 diagrama parabólico-rectangular con n=2. Sin tracción (para εc <= 0 → 0).
    """
    if eps_c <= 0.0:
        return 0.0
    if eps_c <= eps_c2:
        eta = eps_c / eps_c2
        return fcd * (1.0 - (1.0 - eta)**2)  # n = 2
    if eps_c <= eps_cu2:
        return fcd
    # fuera del dominio último
    return 0.0


def sigma_s_bilineal(eps_s: float, Es: float, fyd: float) -> float:
    """Tensión en MPa (signada). Elasto-perfectamente plástico, simétrico."""
    sig = Es * eps_s
    if sig > fyd:
        return fyd
    if sig < -fyd:
        return -fyd
    return sig


def mrd_fibers_rect(
    b_mm: float,
    h_mm: float,
    steel_layers: List[SteelLayer],
    materials: Materials = Materials(),
    nlayers_concrete: int = 400,
    sign: Literal["positivo", "negativo"] = "positivo",
    tol_N: float = 1e-3,
    itmax: int = 120
) -> Tuple[float, float, float]:
    """
    Calcula MRd [kN·m] para una sección rectangular por fibras, EC-2.
    - sign="positivo": comprime la fibra superior (momento hogging/positivo).
    - sign="negativo": comprime la fibra inferior (sagging/negativo).
    Devuelve (MRd_kNm, c_mm, z_equiv_mm) siendo c la profundidad del eje neutro
    medida desde la fibra comprimida.
    """
    mat = materials
    fcd, fyd, Es = mat.fcd, mat.fyd, mat.Es_MPa
    eps_c2, eps_cu2 = mat.eps_c2, mat.eps_cu2

    # Para momento negativo, reflejamos la geometría respecto al borde comprimido
    # y trabajamos siempre "desde arriba" (y=0 fibra comprimida).
    def map_y(y):
        return y if sign == "positivo" else (h_mm - y)

    layers = [SteelLayer(y_mm=map_y(sl.y_mm), As_mm2=sl.As_mm2)
              for sl in steel_layers]

    # Discretización del hormigón en láminas paralelas al borde comprimido
    dy = h_mm / nlayers_concrete
    # 0 (arriba) → h (abajo)
    y_centers = [(i+0.5)*dy for i in range(nlayers_concrete)]

    def section_resultants(c_mm: float):
        """
        Para un c dado: calcula N (positivo a compresión) y M respecto a la fibra comprimida (y=0).
        Supone ε_top = εcu2; eje neutro en y=c.
        """
        kappa = mat.eps_cu2 / max(c_mm, 1e-6)         # curvatura [1/mm]
        # Hormigón
        Nc = 0.0  # N (MPa*mm² = N)
        Mc = 0.0  # Momento respecto a y=0 (N*mm)
        for y in y_centers:
            eps = kappa * (c_mm - y)  # >0 compresión, <0 tracción
            # MPa, solo compresión
            sig = sigma_c_ec2(eps, fcd, eps_c2, eps_cu2)
            Fc = sig * b_mm * dy       # N
            Nc += Fc
            Mc += Fc * (y + 0.0)       # brazo y
        # Acero
        Ns = 0.0
        Ms = 0.0
        for sl in layers:
            eps_s = kappa * (c_mm - sl.y_mm)
            sig_s = sigma_s_bilineal(eps_s, Es, fyd)  # MPa (signado)
            Fs = sig_s * sl.As_mm2                    # N
            Ns += Fs
            Ms += Fs * sl.y_mm
        N = Nc + Ns
        M = Mc + Ms
        return N, M, kappa

    # Búsqueda de c tal que N≈0
    # El rango físico: c in (muy pequeño, ~varias veces h). Si c->0, casi todo tracción en acero;
    # si c grande, casi todo compresión en hormigón y acero superior.
    c_low = 1e-3
    c_high = 3.0 * h_mm
    N_low, _, _ = section_resultants(c_low)
    N_high, _, _ = section_resultants(c_high)
    # Asegurar cambio de signo ampliando por si acaso
    expand = 0
    while N_low * N_high > 0 and expand < 8:
        c_high *= 1.8
        N_high, _, _ = section_resultants(c_high)
        expand += 1

    c = None
    for _ in range(itmax):
        c_mid = 0.5 * (c_low + c_high)
        N_mid, _, _ = section_resultants(c_mid)
        if abs(N_mid) < tol_N:
            c = c_mid
            break
        if N_low * N_mid <= 0:
            c_high = c_mid
            N_high = N_mid
        else:
            c_low = c_mid
            N_low = N_mid
        c = c_mid

    # Con c encontrado, el momento resistente de cálculo es el valor de M (N*mm)
    N_eq, M_eq, kappa = section_resultants(c)
    # (N_eq puede no ser exactamente cero por la tolerancia)
    MRd_kNm = abs(M_eq) / 1e6  # N*mm → kN·m

    # Brazo equivalente aproximado z = M / T_tracción (informativo)
    # Tomamos T como suma de tracciones (Fs<0) en valor absoluto
    T_abs = 0.0
    for sl in layers:
        eps_s = (mat.eps_cu2 / max(c, 1e-6)) * (c - sl.y_mm)
        sig_s = sigma_s_bilineal(eps_s, Es, fyd)
        Fs = sig_s * sl.As_mm2
        if Fs < 0:
            T_abs += -Fs
    z_eq = (abs(M_eq) / T_abs) if T_abs > 1e-9 else float('nan')

    # Nota: si sign="negativo", MRd es el mismo valor (la orientación ya se invirtió).
    return MRd_kNm, c, z_eq


# ========== EJEMPLO DE USO ==========
if __name__ == "__main__":
    # Sección 300x600 mm, recubrimiento nominal 35 mm, capas:
    b, h = 300.0, 600.0
    cnom = 35.0
    # Profundidad de centros de barras (desde borde superior):
    # Superior: 2Ø16 a y = cnom + 16/2
    y_sup = cnom + 16/2
    As_sup = area_barras([16, 16])
    # Inferior: 4Ø16 a y = h - (cnom + 16/2)
    y_inf = h - (cnom + 16/2)
    As_inf = area_barras([16, 16, 16, 16])

    layers = [
        SteelLayer(y_mm=y_sup, As_mm2=As_sup),
        SteelLayer(y_mm=y_inf, As_mm2=As_inf),
    ]

    mats = Materials(fck_MPa=30.0, fyk_MPa=500.0)
    MRd_pos, c_pos, z_pos = mrd_fibers_rect(
        b, h, layers, mats, sign="positivo")
    MRd_neg, c_neg, z_neg = mrd_fibers_rect(
        b, h, layers, mats, sign="negativo")

    print(f"MRd+ = {MRd_pos:.2f} kN·m (c={c_pos:.1f} mm, z≈{z_pos:.1f} mm)")
    print(f"MRd- = {MRd_neg:.2f} kN·m (c={c_neg:.1f} mm, z≈{z_neg:.1f} mm)")
