import sympy as sp
import numpy as np

# variable de Laplace
s = sp.symbols('s')

# Variables électriques
V = sp.Function('V')(s)
I = sp.Function('I')(s)

Re, Le, Bl = sp.symbols('R_e L_e Bl')  # Paramètres TS électriques (resistance, inductance, moteur)
Mm, Rm, Cm = sp.symbols('M_m R_m C_m')  # Paramètres TS mécaniques (masse perte, compliance suspension)

x = sp.Function('x')(s) # position membrane
v = x * s # vitesse membrane
a = v * s # acceleration membrane

Pf = sp.Function('Pf')(s)  # Pression acoustique frontal
Pb = sp.Function('Pb')(s) #pression accoustique arrière

# Paramètres acoustiques (densité, vitesse, surface membrane, perte visqueuses, compliance et masse d'air)
rho, c, Sd = sp.symbols('rho c Sd')
R_af, C_af, L_af, R_ab, C_ab, L_ab = sp.symbols('R_af C_af L_af R_ab C_ab L_ab')  

F_motor = Bl * I
F_mech = Mm * a + Rm * v + (1/Cm) * x
F_acous = Sd * (Pf - Pb) #pas sûr du signe

# Équation électrique
eq_elec = sp.Eq(V, Le * s * I + Re * I  - Bl * v)
eq_meca = sp.Eq(0, F_motor + F_mech - F_acous)

# Équations acoustiques (impédances acoustiques)
eq_acousf = sp.Eq(Pf,R_af * v + L_af * a + (1/C_af) * x)
eq_acousb = sp.Eq(Pb,R_ab * v + L_ab * a + (1/C_ab) * x)

symbols = [Pf, Pb, x, I, V]

def compute_transfer(symbol_p,symbol_d):
    # Construire la matrice des équations
    s = [s for s in symbols if s != symbol_d]
    matrix, constants = sp.linear_eq_to_matrix([eq_acousf,
                                                         eq_acousb,
                                                         eq_elec,
                                                         eq_meca], s)
    inverse_matrix = matrix.inv()
    solutions = inverse_matrix @ constants

    solutions_dict = dict(zip(s, solutions))
    transfer_function = (solutions_dict[symbol_p] / symbol_d)

    return transfer_function.simplify()

transfer_function = compute_transfer(Pf,V)

d={ #beyma 12br70
    'Bl' : 12.1,
    'Fs' : 31,
    'Re' : 5.6,
    'Le' : 0.8e-3,
    'Cms' : 345e-6,
    'Rms' : 3.3,
    'Mms' : 0.074,
    'Sd' : 0.054,

    'c_air' : 343.21,
    'rho' : 1.2041,
    #box
    'Vb' : 50 / 1e3,
}


#pas sûr de ces valeurs
a = (Sd / sp.pi)**0.5  # rayon effectif de la membrane (m)
l_af = 0.6 * a  # correction d’extrémité typique
l_ab = 0.6 * a  # longueur du volume arrière (exemple)
L_af = rho * l_af / Sd  # inertance frontale
L_ab = rho * l_ab / Sd  # inertance arrière

d_ = {
    "R_e" : d['Re'],
    "L_e" : d['Le'],
    "C_m":  d['Cms'],

    "C_ab":d['Vb'] / ( d['rho'] * d['c_air'] **2 ),
    #"C_af":  sp.oo # grosse valeur à substituer avec sp.limit, elasticité volume extérieur

    "L_af":L_af,
    "L_ab":L_ab,

    "R_af":rho*d['c_air']/(sp.pi*a**2),

    "R_ab":0,#pour simplifier
    "M_m" :d['Mms'], #masse mécanique
    "R_m" :d['Rms']
}

transfer_function = transfer_function.subs(d_)
#for k,v in d_.items():
#   print(f"substition de {k} = {v}")
#   transfer_function = transfer_function.subs({k:v})
#   sp.pprint(transfer_function.simplify(),num_columns=200)

transfer_function = sp.simplify(sp.limit(transfer_function,C_af,sp.oo))

transfer_function = transfer_function.subs(d)
sp.pprint(transfer_function,num_columns=200)

#print(transfer_function.simplify())

H_s = transfer_function
 # Transforme en fonction pour numpy

# Fonction pour évaluer la fonction de transfert en jω (fréquence complexe)
def evaluate_transfer_function(H_s, f):
    H_lambdified = sp.lambdify(s, H_s, 'numpy')
    omega = f * (2 * np.pi)
    s_value = 1j * omega
    return H_lambdified(s_value)

# Générer un ensemble de fréquences pour le diagramme de Bode
frequency_vals = np.logspace(1, 4, 1000)
magnitude = []
phase = []


H_value = evaluate_transfer_function(transfer_function, frequency_vals)
magnitude = 20 * np.log10(np.abs(H_value))  # Magnitude
phase = np.angle(H_value) # Phase en degrés


if 1:
    import matplotlib.pyplot as plt
    # Créer le diagramme de Bode
    fig, axs = plt.subplots(2, 1, figsize=(10, 6))
    # Tracer la magnitude
    axs[0].semilogx(frequency_vals, magnitude)
    axs[0].set_ylabel('amplitude (db)')
    axs[0].set_title('Diagramme de Bode')
    axs[0].xaxis.set_major_formatter(plt.ScalarFormatter())
    axs[0].grid()

    max_mag = np.max(magnitude)  # valeur maximale en dB
    axs[0].axhline(y=max_mag, color='r', linestyle='--', linewidth=1.5, label=f'Max = {max_mag:.1f} dB')
    axs[0].axhline(y=max_mag - 3, color='g', linestyle='--', linewidth=1.5, label=f'Max - 3 dB = {max_mag - 3:.1f} dB')

    # Tracer la phase
    axs[1].semilogx(frequency_vals, np.unwrap(phase))
    axs[1].set_ylabel('Phase (rad)')
    axs[1].set_xlabel('Fréquence (rad/s)')

    # Affichage du graphique
    plt.tight_layout()
    plt.show()