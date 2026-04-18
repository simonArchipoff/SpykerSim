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

def parametres_second_ordre(H, s):
    """
    Calcule la pulsation propre omega et le facteur de qualité Q
    à partir d'une fonction de transfert H(s) du second ordre.

    Paramètres:
    H : expression sympy rationnelle (numérateur/dénominateur)
    s : variable symbolique (souvent sp.Symbol('s'))

    Retourne:
    (omega, Q) : deux expressions sympy
    """
    # Extraire le dénominateur
    num, den = sp.fraction(H)
    poly = sp.Poly(den, s)

    if poly.degree() != 2:
        raise ValueError("Le dénominateur n'est pas du second ordre")

    # Coefficients : c2*s^2 + c1*s + c0
    coeffs = poly.all_coeffs()  # [c2, c1, c0]
    c2, c1, c0 = coeffs

    # Pulsation propre
    omega = sp.sqrt(c0 / c2)

    # Facteur de qualité
    # Q = sqrt(c0*c2) / c1, avec gestion du cas c1=0 (Q infini)
    if c1 == 0:
        Q = sp.oo
    else:
        Q = sp.sqrt(c0 * c2) / c1

    return omega, Q