import numpy as np
import sympy as sp
def pa_to_db_spl(p_rms, p_ref=2e-5):
    """Convertit une pression efficace (Pa) en dB SPL."""
    return 20 * np.log10(p_rms / p_ref)

def compute_spl(frequency_hz, x_peak_m, Sd=0.0309, rho=1.2041, r_m=1.0, half_space=True,p_ref=2e-5):
    """
    Calcule le SPL à r_m mètres.
    half_space=True -> rayonnement dans un demi-espace (baffle ou sol)
    half_space=False -> espace libre (4π)
    """
    omega = 2 * np.pi * frequency_hz
    a_peak = omega**2 * x_peak_m
    a_rms = a_peak / np.sqrt(2)
    if half_space:
        p_rms = (rho * Sd * a_rms) / (2 * np.pi * r_m)
    else:
        p_rms = (rho * Sd * a_rms) / (4 * np.pi * r_m)
    spl = 20 * np.log10(p_rms / p_ref)
    return spl


def normalize_second_order(H, s):
    """
    Met une fonction de transfert rationnelle H(s) sous forme canonique
    pour un système du second ordre passe-bas :
        H(s) = G / (s**2 + (omega0/Q)*s + omega0**2)
    ou, si le numérateur a des zéros, sous la forme générale :
        H(s) = G * (s**2 + ...) / (s**2 + (omega0/Q)*s + omega0**2)

    Paramètres:
        H : expression sympy rationnelle en s
        s : symbole de Laplace

    Retourne:
        un tuple (G, omega0, Q, num_factors) où num_factors est une liste
        de termes (gain, zéro) pour les zéros éventuels.
    """
    # Mettre sous forme de fraction
    H_num, H_den = sp.fraction(sp.simplify(H))

    # Assurer que le dénominateur est un polynôme en s
    if not H_den.is_polynomial(s):
        raise ValueError("Le dénominateur n'est pas un polynôme en s")

    # Obtenir les coefficients du dénominateur (ordre décroissant)
    coeffs = sp.Poly(H_den, s).all_coeffs()
    order = len(coeffs) - 1
    if order != 2:
        raise ValueError(f"Le dénominateur n'est pas du second ordre (ordre {order})")

    # Normaliser le dénominateur pour que le coefficient de s^2 soit 1
    a2, a1, a0 = coeffs  # a2 * s^2 + a1 * s + a0
    if a2 == 0:
        raise ValueError("Coefficient de s^2 nul")

    den_norm = sp.Poly(H_den / a2, s)  # s^2 + (a1/a2) s + (a0/a2)
    # Extraire les coefficients normalisés
    _, b1, b0 = den_norm.all_coeffs()

    # Calculer omega0 et Q
    # b0 = omega0^2, b1 = omega0 / Q
    omega0 = sp.sqrt(b0)
    Q = omega0 / b1 if b1 != 0 else sp.oo

    # Gain global : G = H_num / a2
    G = sp.simplify(H_num / a2)

    # Gérer les zéros : factoriser le numérateur
    num_poly = sp.Poly(H_num, s)
    num_factors = []
    if num_poly.degree() > 0:
        # Factoriser le numérateur pour extraire les zéros
        factors = sp.factor(num_poly).as_terms()
        # Simplification : on peut extraire les racines
        roots = sp.roots(num_poly, s)
        for r, mult in roots.items():
            num_factors.append((r, mult))

    return G, omega0, Q, num_factors