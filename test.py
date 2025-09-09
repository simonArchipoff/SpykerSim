
import sympy as sp
from lcapy import transfer
import matplotlib.pyplot as plt
import numpy as np
# Définition des variables symboliques
s = sp.symbols('s')  # Variable de Laplace
X = sp.Function('X')(s)  # Déplacement de la membrane en Laplace
V = sp.Function('V')(s)  # Vitesse de la membrane en Laplace (s*X)
A = sp.Function('A')(s)  # Accélération de la membrane en Laplace (s^2*X)
P = sp.Function('P')(s)  # Pression acoustique en Laplace
F = sp.Function('F')(s)  # Force appliquée en Laplace

# Paramètres mécaniques et acoustiques
Mm, Rm, Cm = sp.symbols('Mm Rm Cm')  # Masse, résistance, compliance mécaniques
Sd, Ra, Ca = sp.symbols('Sd Ra Ca')  # Surface, résistance, compliance acoustiques

# Relation entre déplacement, vitesse et accélération
V = s * X
A = s**2 * X

# Équation mécanique de la membrane dans le domaine de Laplace
F_mech = Mm * A + Rm * V + (1 / Cm) * X
F_acous = Sd * P
eq_meca = sp.Eq(F, F_mech - F_acous)

# Équation acoustique dans le domaine de Laplace
P_eq = sp.Eq(F, Ra * Sd * V + (Sd / Ca) * X)




matrix, constants = sp.linear_eq_to_matrix([P_eq,
                                                     eq_meca], [ P, X])



inverse_matrix = matrix.inv()
solutions = inverse_matrix @ constants

transfer_function = solutions[0] / F
print(transfer_function)


d={ #beyma 12br70
    'Bl' : 12.1,
    'Fs' : 31,
    'Re' : 5.6,
    'Le' : 0.8e-3,
    'Qms' : 4.4,
    'Qes' : 0.44,
    'Qts' : 0.50,
    'Vas' : 142 / 1e3,
    'Cms' : 345e-3,
    'Rms' : 3.3,
    'Mms' : 0.074,
    'Sd' : 0.054,

    'c_air' : 343.21,
    'rho' : 1.2041,
    #box
    'Vb' : 50 / 1e3,
}

d_ = {
    "Re" :d['Re'],
    "Le" : d['Le'],
    "Cm":d['Cms'],
    "Ca":d['Vb'] / ( d['rho'] * d['c_air'] * d['c_air']),
    "Ma": 0.0 * d['Sd'] * 0.1 * d['rho'],
    "Mm" :d['Mms'],
    "Ra" : 1,

}

tf = transfer_function.subs(d)
tf = tf.subs(d_)
Rms_expr = (d['Fs'] * d['rho'] * d['Sd']) / (sp.pi * d['Qms'] * d['Mms'])
tf = tf.subs(Rm, Rms_expr)
tf = tf.subs(sp.pi,np.pi)


# Fonction pour évaluer la fonction de transfert en jω (fréquence complexe)
def evaluate_transfer_function(H_s, f):
    H_lambdified = sp.lambdify(s, H_s, 'numpy')
    omega = f / (2 * np.pi)
    s_value = 1j * omega
    return H_lambdified(s_value)

# Générer un ensemble de fréquences pour le diagramme de Bode
frequency_vals = np.logspace(0.0001, 4, 1000)  # Fréquences de 0.1 à 1000 rad/s
magnitude = []
phase = []


H_value = evaluate_transfer_function(tf.simplify(), frequency_vals)
magnitude = np.abs(H_value)  # Magnitude
phase = np.angle(H_value) # Phase en degrés

# Convertir la magnitude en dB
magnitude_dB = 20 * np.log10(magnitude)

# Créer le diagramme de Bode
fig, axs = plt.subplots(2, 1, figsize=(10, 6))
# Tracer la magnitude
axs[0].semilogx(frequency_vals, magnitude_dB)
axs[0].set_ylabel('Magnitude')
axs[0].set_title('Diagramme de Bode')

# Tracer la phase
axs[1].semilogx(frequency_vals, np.unwrap(phase))
axs[1].set_ylabel('Phase (rad)')
axs[1].set_xlabel('Fréquence (rad/s)')



# Affichage du graphique
plt.tight_layout()
plt.show()
exit(0)

import sympy as sp
import matplotlib.pyplot as plt
import numpy as np

# Variables symboliques pour le temps et la transformée de Laplace
s = sp.symbols('s')

# Variables électriques
F = sp.Function('F')(s)
P = sp.Function('P')(s)
Mm, Rm, Cm = sp.symbols('M_m R_m C_m')  # Paramètres TS mécaniques (masse perte, suspension
x = sp.Function('x')(s)
rho, c, Sd, R_a, C_a = sp.symbols('rho c Sd R_a C_a')  # Paramètres acoustiques (densité, vitesse, surface membrane, perte visqueuses, compliance)
v = x * s
a = v * s

a = F / Mm
v = a /s
x = v / s

F_mech = Mm * a + Rm * v + (1 / Cm) * x
F_acous = Sd * P
F_membrane = F_mech + F_acous + F
eq_meca = sp.Eq(F_membrane, Mm * a)

# Équations acoustiques (impédances acoustiques)
eq_acous = sp.Eq(P, (R_a * v + v / C_a))

# Construire la matrice des équations
matrix, constants = sp.linear_eq_to_matrix([eq_acous,
                                                     #eq_elec,
                                                     eq_meca], [ P, F])

inverse_matrix = matrix.inv()
solutions = inverse_matrix @ constants

# Fonction de transfert : P(s) / V(s)
transfer_function = (solutions[0] / F)

print(transfer_function)


