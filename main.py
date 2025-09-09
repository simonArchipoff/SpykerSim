import lcapy as lc
import matplotlib.pyplot as plt
from lcapy import Circuit
import numpy as np

#https://jahonen.kapsi.fi/Audio/Papers/enclosuremodelling.pdf
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

a = lc.Circuit("""
    Re in+ 3 ; right
    Le in- 2 ; right
    P1 in+ in-  ; down
    #Fm 5 4  {Bl}  ; up
    #Gm 3 2 {1/Bl} ; up 
    TF1   4 5 3 2 {1/Bl} ; right 
    #GY1 4 5 3 2 {Bl}; right

    Rms 5 9 {Rms} ; right
    Lms 4 7 {Mms} ; right
    Cms 7 8 {Cms} ; right
    TF2 10 11 8 9 {1/Sd} ; right
    #W 10 11 ; down
    # W 8 9_0 ; down
    # W 9 9_0 ; right   
    Ral 10 11 ; down 
    
    W 11 11_0 ; down
    W 0 11_0 ; right
    
    Rab 10 14 {1176} ; right
    Ccab 14 15 {Vb/(rho*c_air*c_air)}; right
    

    #Ral 11 13 ; right
    W 15 16 ; down
    W 16 11 ; left
    
    
    #connect 0
    W in- 0 ; down
    W 5 0_5 ; down
    #W 11 0_11 ; down
    
    W 0 0_5; right
    #W  0_5 0_11 ; right
    
    #W 110 15 ; down
    #W 17 13_0 ; up
    
    #LMap 15 16 {LMap}; right
    #Rap 16 17 {Rap} ; right
    
    Haccoustic_p 0 pCcab Ccab 1 ; down
    #Hcone 0 cone_excur Lms 1 ; down
    #W 0 0_e ; right
    #Hport_p 0_e pLMap LMap 1 ; down
   """)
a.draw()
plt.show()

#print(a.transfer("P1",("pCcab","0")))
#print(a.transfer("P1",("pLMap","0")))


fsc =d['Fs'] * np.sqrt(1+(d['Vas'] / d['Vb']))
wsc = fsc * 2 * np.pi
V = (d['Vb'] * d['Vas']) / (d['Vb'] + d['Vas'])
vCcab = d['Vb'] / (d['rho'] * d['c_air']**2)

#valeur calculée à partir de Q
vRal = 10 / (wsc * (V / (d['rho'] * d['c_air']**2)))
vRab = 1  / (wsc  * 100 * (V / (d['rho'] * d['c_air']**2)) )

print(f"{V=} {vCcab=},{fsc=},{vRal=}, {vRab=}")
G = (a.transfer("P1",("pCcab","0"))) * lc.s
#E = (a.transfer("P1",("cone_excur","0")))
print(G.simplify())
G = G.subs(d)
G = G.subs({"Ral":vRal,"Rab":vRab})
G(lc.jf).dB.plot((30,200))



#plt.plot(np.abs(E.evaluate(range(2000))))
#plt.plot( 20 * np.log10( np.abs(G.evaluate(range(2000)))))
plt.show()
