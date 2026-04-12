
import sympy as sp
import numpy as np

# Variables symboliques pour le temps et la transformée de Laplace
t = sp.symbols('t', real=True)
s = sp.symbols('s')

# Variables électriques
V = sp.Function('V')(s)
I = sp.Function('I')(s)
Re, Le, Bl = sp.symbols('R_e L_e Bl')  # Paramètres TS électriques (resistance, inductance, moteur)
Mm, Rm, Cm = sp.symbols('M_m R_m C_m')  # Paramètres TS mécaniques (masse perte, suspension
x = sp.Function('x')(s) # position membrane
P = sp.Function('P')(s)  # Pression acoustique
rho, c, Sd, R_a, C_a = sp.symbols('rho c Sd R_a C_a')  # Paramètres acoustiques (densité, vitesse, surface membrane, perte visqueuses, compliance)
v = x * s # vitesse membrane
a = v * s # acceleration membrane

F_motor = Bl * I
F_mech = Mm * a + Rm * v + (1 / Cm) * x
F_acous = Sd * P

# Équation électrique
eq_elec = sp.Eq(V, Le * s * I + Re * I  - Bl * v)
# Équation globale de la membrane (somme des forces = M_m * a)
eq_meca = sp.Eq(0, F_motor + F_mech - F_acous)
# Équations acoustiques (impédances acoustiques)
eq_acous = sp.Eq(P, (R_a * v + v / C_a))

symbols = [P, x, I, V]

def compute_transfer(symbol_p,symbol_d):
    # Construire la matrice des équations
    s = [s for s in symbols if s != symbol_d]
    matrix, constants = sp.linear_eq_to_matrix([eq_acous,
                                                         eq_elec,
                                                         eq_meca], s)

    inverse_matrix = matrix.inv()
    solutions = inverse_matrix @ constants

    solutions_dict = dict(zip(s, solutions))
    # Fonction de transfert : P(s) / V(s)
    transfer_function = (solutions_dict[symbol_p] / symbol_d)

    #impedance =  V/solutions_dict[I]

    print("\nFonction de transfert (P(s) / V(s)) :")
    sp.pprint(transfer_function.simplify(),num_columns=200)
    return transfer_function.simplify()

transfer_function = compute_transfer(P,I)

d={ #beyma 12br70
    'Bl' : 12.1,
    'Fs' : 31,
    'Re' : 5.6,
    'Le' : 0.8e-3,
    'Qms' : 4.4,
    'Qes' : 0.44,
    'Qts' : 0.50,
    'Vas' : 142 / 1e3,
    'Cms' : 345e-6,
    'Rms' : 3.3,
    'Mms' : 0.074,
    'Sd' : 0.054,

    'c_air' : 343.21,
    'rho' : 1.2041,
    #box
    'Vb' : 50 / 1e3,
}

#Vas = Cms × d × c² × A²
Cms = d['Vas'] / ( d['Sd'] ** 2 * d['rho'] * d['c_air'] ** 2 )
f = 205e-3 / (0.1255**2 * d['rho'] * d['c_air'] ** 2)

d_ = {
    "R_e" :d['Re'],
    "L_e" : d['Le'],
    "C_m":d['Cms'],
    "C_a":d['Vb'] / ( d['rho'] * d['c_air'] * d['c_air']),
    "M_a":  d['Sd'] * 0.1 * d['rho'],
    "M_m" :d['Mms'],
    "R_a" : 10,

}

print(transfer_function)

#transfer_function = transfer_function.subs(Rm, (1 / d['Qms']) * sp.sqrt(Mm / Cm) )
Rms_expr = (d['Fs'] * d['rho'] * d['Sd']) / ( d['Qms'] * Mm)
Rms_expr = Rms_expr.subs(d_)

transfer_function = transfer_function.subs(Rm, Rms_expr)
transfer_function = transfer_function.subs(d)
transfer_function = transfer_function.subs(d_)
print(transfer_function.simplify())

H_s = transfer_function
 # Transforme en fonction pour numpy

# Fonction pour évaluer la fonction de transfert en jω (fréquence complexe)
def evaluate_transfer_function(H_s, f):
    H_lambdified = sp.lambdify(s, H_s, 'numpy')
    omega = f / (2 * np.pi)
    s_value = 1j * omega
    return H_lambdified(s_value)

# Générer un ensemble de fréquences pour le diagramme de Bode
frequency_vals = np.logspace(0.00000000001, 4, 1000)  # Fréquences de 0.1 à 1000 rad/s
magnitude = []
phase = []


H_value = evaluate_transfer_function(transfer_function, frequency_vals)
magnitude = np.abs(H_value)  # Magnitude
phase = np.angle(H_value) # Phase en degrés

# Convertir la magnitude en dB
magnitude_dB = 20 * np.log10(magnitude)


if 0:
    import matplotlib.pyplot as plt
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