from lcapy import *
import matplotlib.pyplot as plt

a = Circuit("""
W in 1; right
R1 1 2 {Rb} ; right
C1 2 0 {Cb}; down
R2 2 3 {Rb}; right


W 0 0_1; right
W 0_1 0_2 ; right
W 0_2 0_3 ; right

W  1 1_u ; up
R3 1_u 4 {Ra} ; right
C2 4 3_u {Ca} ; right
W 3_u 3 ; down

W 3 3_1 ; right
W 3_1 3_2 ; right


E 12 0 opamp 0_3 3_2 A; right, mirror
W 12 out; right


W 3_u 6 ; up
R4 6 7 {Rc} ; right
C3 7 0_c3 {Cc} ; down 
R5 7 8 {Rc} ; right
W 8 8_out ; right
W 8_out 12; down
W 0_c3 0_2 ; down

W 6 9 ; up
R6 9 10 {Ra}; right
C4 10 11 {Ca} ; right
W 11 8 ; down
""")

fi = 122
qi = 1.53

fo = 20
qo = 0.8
ideal = (s**2 + (fi/qi) * s + fi**2) / (s**2 + (fo/qo) * s + fo**2)
k = ((fi / fo) - (qi / qo)) / ((qi / qo) - (fo / fi))
ca = 0.01e-6
rb = 1 / (2 * pi * fi * ca * (2 * qi * (1 + k)))
ra = 2 * k * rb
rc = rb * (fi / fo)**2
cb = ca * (2 * qi * (1+k))**2
cc = cb * (fo / fi)**2


a = a.subs({"Rc":rc.evalf(),"Ra":ra.evalf(),"Rb":rb.evalf(),"Ca":float(ca), "Cb":float(cb), "Cc":float(cc)})

a.draw()
plt.show()


h = -a.transfer(0,1,0, 12).limit("A",oo).simplify()
f = ideal - h(s*2*pi)
print(h)
#ideal(jw).plot((20, 500), log_frequency=True)
h(jw*2*pi).plot((20,500), log_frequency=True)

# I dont know if this difference (about 1e-14) comes from the numerical computation or the circuit itself
f(jw).plot((20,500),log_frequency=True)

plt.show()
