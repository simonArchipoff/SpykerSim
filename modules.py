import sympy as sp


s = sp.Symbol('s')
rho, c = sp.symbols('rho_air c_air')


class Module:
    def __init__(self,unique_id):
        self._id = unique_id
        self._local_symbols = []
        self._local_equations = []
        self._global_functions = []
        self._values = {}

    def get_local_global_map(self):
        return {sp.symbol(x.name): sp.symbol(self.g(x.name)) for x in self._symbols}

    def get_symbols_global(self):
        return self.get_local_global_map().values()

    def get_subs_values(self):
        return {f"{g.name}": self._values[l.name] for l,g in self.get_local_global_map().items() if l.name in self._values  }

    def get_global_equations(self):
        d = self.get_local_global_map()
        return [x.subs[d] for x in self._local_equations]

    def g(self,s):
        return f"{s}_{self._id}"

def ClosedBox(Module):
    def __init__(self, unique_id):
        super().__init__(self.unique_id)
        R_a, C_a, L_a = sp.symbols('R_a C_a L_a')
        self.P = sp.symbols("P")
        self.z  = (R_a * s + L_a * s * s  + (1 / C_a))


class Speaker(Module):
    def __init__(self,unique_id):
        super().__init__(unique_id)
        # Variables électriques
        self.V = sp.symbols(self.g("V"))
        self.I = sp.symbols(self.g("I"))

        Re, Le, Bl = sp.symbols('R_e L_e Bl')  # Paramètres TS électriques (resistance, inductance, moteur)
        Mm, Rm, Cm = sp.symbols('M_m R_m C_m')  # Paramètres TS mécaniques (masse perte, compliance suspension)
        Sd = sp.symbols('Sd') # surface piston
        self._local_symbols += [Re, Le, Bl, Mm, Rm, Cm,Sd]

        self.x = sp.symbols(self.g("x"))  # position membrane

        v = self.x * s  # vitesse membrane
        a = v * s  # acceleration membrane

        Pf = sp.symbols(self.g("Pf"))  # Pression acoustique frontal
        Pb = sp.symbols(self.g("Pb"))  # pression accoustique arrière

        # Paramètres acoustiques (densité, vitesse de l'air, surface membrane, perte visqueuses, compliance et masse d'air)
        F_motor = Bl * self.I
        F_mech = Mm * a + Rm * v + (1 / Cm) * self.x
        F_acous = Sd * (Pf - Pb)

        eq_elec = sp.Eq(self.V, Le * s * self.I + Re * self.I  + Bl * v)
        eq_meca = sp.Eq(F_motor ,  F_mech + F_acous)

        self.zaf, self.zab = (sp.symbols(self.g("zaf")), sp.symbols(self.g("zab")))
        eq_acousf = sp.Eq(Pf, Sd * self.x * self.zaf)
        eq_acousb = sp.Eq(Pb, Sd * self.x * self.zab)

        self._equations = [eq_elec, eq_meca, eq_acousf, eq_acousb]



if __name__ == "__main__":
    s = Speaker(1)
    pass