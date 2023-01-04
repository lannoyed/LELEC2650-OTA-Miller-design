
"""
 Python sizing

50nm CMOS technology (TSMC)

Miller Design algorithm" of project of ELEC2650 (2022-23)

Sophie Alicandro, Julien Giunta, Diego Lannoye
From Sylvain Favresse, 2022
From Marco Gonzalez's work on MATLAB (2021)

The code will not run until the design algorithm is complete
"""

import numpy as np
import matplotlib.pyplot as plt
import settings
from plot_techno_curves import plot_techno_curves
import scipy.interpolate as scint
import scipy.signal as ss
import control


# # ---------------------------------------------------------------------------
# # ------------------------- Library extraction ------------------------------
# # ---------------------------------------------------------------------------

nch_lvt = settings.read_txt('data/nch_lvt.txt')
nch = settings.read_txt('data/nch.txt')
nch_hvt = settings.read_txt('data/nch_hvt.txt')
pch_lvt = settings.read_txt('data/pch_lvt.txt')
pch = settings.read_txt('data/pch.txt')
pch_hvt = settings.read_txt('data/pch_hvt.txt')

# # Technological curves
# plot_techno_curves([nch_lvt, nch, nch_hvt],['nch_lvt','nch','nch_hvt'])
# plot_techno_curves([pch_lvt, pch, pch_hvt],['pch_lvt','pch','pch_hvt'])
# plt.show()

# # ---------------------------------------------------------------------------
# # ----------------------------- Specifications ------------------------------
# # ---------------------------------------------------------------------------




#ANNONCE

#Miller
M1 = pch_lvt; M2 = pch_lvt
M3 = nch_lvt; M4 = nch_lvt
M5 = pch_lvt; M6 = pch_lvt
M7 = nch_lvt; M8 = pch_lvt


#Boost
M17 = pch_lvt
M16 = nch_lvt
M10=pch_lvt

#FF

M11=pch_lvt
M12=nch_lvt
M13=nch_lvt
M14=nch_lvt
M15=pch_lvt

#Je sais plus pourquoi mais m'en demande pas trop : ok si c'est pcq boucle
global gm1,gm2,gm3,gm4,gm5,gm6,gm7,gm8
global g1,g2,g3,g4,g5,g6,g7,g8
global id1,id2,id3,id4,id5,id6,id7,id8
global idn1,idn2,idn3,idn4,idn5,idn6,idn7,idn8
global vsg1,vsg2,vgs3,vgs4,vsg5,vsg8,vsg6,vgs7


global VIN

def abs_sum(x):
    return np.sum(np.abs(x))

def Result_is_physical(Cc1, Cc2, CL, ibias, VIN, W_array) : 
    if Cc1 < 0 or Cc2 < 0 or CL < 0 or ibias < 0 or VIN < 0 < 0:
        return 0
    if VIN > 1.2:
        return 0
    if ibias > 1e-2:
        return 0
    for W in W_array:
        if W < 0 or W > 1000:
            return 0
    return 1 

def transitor_are_satured (vds_list, vdd):
    # faire -1 au numero pour les transistore 1 à 8
    # faire -2 au numero pour les transistore 10 à 17
    
    # M1 M3 M5
    if abs_sum(np.array([vds_list[0], vds_list[2], vds_list[4]])) > vdd:
        return 0
    # M2 M4 M5
    if abs_sum(np.array([vds_list[1], vds_list[3], vds_list[4]])) > vdd:
        return 0
    # M6 M7
    if abs_sum(np.array([vds_list[5], vds_list[6]])) > vdd:
        return 0
    """ 
    # M16 M17
    if abs_sum(np.array([vds_list[14], vds_list[15]])) > vdd:
        return False
     #M10 M12
    if abs_sum(np.array([vgs_list[8], vgs_list[10]])) > VIN:
        return False
    #M11 M12
    if abs_sum(np.array([vgs_list[9], vgs_list[10]])) > VIN:
        return False
    #M15 et Max M12, M13, M14
    if abs_sum(np.array([vds_list[13], np.max( [vds_list[10], vds_list[11], vds_list[12]]) ] )) > vdd:
        return False
    """
    
    return 1

def Miller_Sizing (fT, CL, Cf, Avout, gmid1, gmid3, gmid5, gmid6, gmid12, gmid16, gmid17, B, B1213, B314) : 

    Flag_Correct_Saturation = 0
    Flag_Correct_Design = 0
    Flag_Mathematical_integrity = 1
    
    
    gmid2 = gmid1
    gmid4 = gmid3
    gmid7 = gmid3
    
    gmid8 = gmid5
    
    Cfavant = Cf/1000.0
    Cfchoix = Cf
    Cfavantchoix = Cfavant
    i = 0
    L = 1e-6 #pour toutes les longueurs

    L1 = L; L2 = L1;
    L3 = L; L4 = L3; 
    L5 = L; L6 = L5;
    L7 = L; 
    L8 = L;
    L9=L10=L11=L12=L13=L14=L15=L16=L17=L; 
    
    while np.abs((Cf-Cfavant)/Cf) > 0.000001: #until Cf converges
        i+=1
        Cfavant=Cf
    # # ---------------------------------------------------------------------------
    # # -------------------- Design algorithm - Miller ----------------------------
    # # ---------------------------------------------------------------------------

    # Differential pair M1, M2
        gm1 = 2*np.pi * fT * Cf
        id1 = gm1 /gmid1
        idn1 = float(scint.interp1d(M1['GMID'], M1['IN'])(gmid1))
        vsg1 = float(scint.interp1d(M1['GMID'], M1['VGS'])(gmid1))
        W_over_L1 = id1 / idn1
        W1 = W_over_L1 * L1
        W2 = W1
        id2 = id1
        vsg2 = vsg1

    #Loading from M3 and M4
        id3= id1
        id4= id1
        idn3 = float(scint.interp1d(M3['GMID'], M3['IN'])(gmid3))
        vgs3 = float(scint.interp1d(M3['GMID'], M3['VGS'])(gmid3))
        vgs4=vgs3
        W_over_L3 = id3 / idn3
        W3 = W_over_L3 * L3
        W4 = W3

    #Biasing M5
        id5 = id1 + id2

    #source en miroir M5 M6 M8
        id8=id5 #B
        ibias = id8
        idn5 = float(scint.interp1d(M5['GMID'], M5['IN'])(gmid5))
        idn8 = float(scint.interp1d(M8['GMID'], M8['IN'])(gmid8))
        vsg5 = float(scint.interp1d(M5['GMID'], M5['VGS'])(gmid5))
        vsg8=vsg5
        W_over_L5 = id5 / idn5
        W5 = W_over_L5 * L5
        W8_over_L8 = id8/  idn8
        W8 = W8_over_L8 * L8


    #transi de Miller M7 + son bias
        gm7 = 10 *gm1
        id7 = gm7/gmid7
        id6=id7
        idn7 = float(scint.interp1d(M7['GMID'], M7['IN'])(gmid7))
        idn6 = float(scint.interp1d(M6['GMID'], M6['IN'])(gmid6))
        vgs7 = float(scint.interp1d(M7['GMID'], M7['VGS'])(gmid7))
        vsg6 = float(scint.interp1d(M6['GMID'], M6['VGS'])(gmid6))
        W_over_L7 = id7 / idn7
        W7 = W_over_L7 * L7
        W_over_L6 = id6 / idn6
        W6 = W_over_L6 * L6

        vea1 = float(scint.interp1d(M1['GMID'], M1['VEA'])(gmid1))
        vea7 = float(scint.interp1d(M7['GMID'], M7['VEA'])(gmid7)) 

    #besoin pour ma fonction de transfert
        gm4=gmid4*id4
        gm1 = gmid1*id1
        g1 = id1/vea1
        g7 = id7/vea7
        gm7 = gmid7*id7



    # # ---------------------------------------------------------------------------
    # # ------------------------- Gain, poles and zeros ---------------------------
    # # ---------------------------------------------------------------------------

        Cgs4 = 2/3 * float(scint.interp1d(M4['GMID'], M4['CGS'])(gmid4)) * W4 * L4
        Cgs6 = 2/3 * float(scint.interp1d(M6['GMID'], M6['CGS'])(gmid6)) * W6 * L6
        Cgs7 = 2/3 * float(scint.interp1d(M7['GMID'], M7['CGS'])(gmid7)) * W7 * L7
        Cgs8 = 2/3 * float(scint.interp1d(M8['GMID'], M8['CGS'])(gmid8)) * W8 * L8
        
        Cgso4 = float(scint.interp1d(M4['GMID'], M4['CGS0'])(gmid4)) * W4
        Cgso6 = float(scint.interp1d(M6['GMID'], M6['CGS0'])(gmid6)) * W6
        Cgso7 = float(scint.interp1d(M7['GMID'], M7['CGS0'])(gmid7)) * W7
        Cgso8 = float(scint.interp1d(M8['GMID'], M8['CGS0'])(gmid8)) * W8
        
        Cgdo2 = float(scint.interp1d(M2['GMID'],M2['CGD0'])(gmid2)) * W2
        Cgdo4 = float(scint.interp1d(M4['GMID'],M4['CGD0'])(gmid4)) * W4
        Cgdo5 = float(scint.interp1d(M5['GMID'],M5['CGD0'])(gmid5)) * W5
        Cgdo6 = float(scint.interp1d(M6['GMID'],M6['CGD0'])(gmid6)) * W6
        Cgdo8 = float(scint.interp1d(M8['GMID'],M8['CGD0'])(gmid8)) * W8
        
        Cbd2 = float(scint.interp1d(M2['GMID'], M2['CBD'])(gmid2)) * W2
        Cbd4 = float(scint.interp1d(M4['GMID'], M4['CBD'])(gmid4)) * W4
        Cbd5 = float(scint.interp1d(M5['GMID'], M5['CBD'])(gmid5)) * W5
        Cbd7 = float(scint.interp1d(M7['GMID'], M7['CBD'])(gmid7)) * W7

        cpp = Cgs4+Cgs6+Cgso4+Cgso6+Cgdo2+Cgdo6+Cbd2+Cbd4
        cpn = Cgs7+Cgs8+Cgso7+Cgso8+Cgdo5+Cgdo8+Cbd7+Cbd5 
        
        
        
        vea6 = float(scint.interp1d(M6['GMID'], M6['VEA'])(gmid6))
        vea8 = float(scint.interp1d(M8['GMID'], M8['VEA'])(gmid8))
        vea_eq = 1 / (1 / vea6 + 1 / vea8)
        
        C1 = Cgs7+Cgso7+Cbd2+Cgdo2+Cbd4+Cgdo4
        Cf = 0.11*(C1+CL+ np.sqrt((C1+CL)**2+2*C1*(CL/0.11)))
        

    #print('Cf = ', Cf, 'F')
        
        
        
    # print(Cf)

    #fonctin de transfert Miller ( formule slides Miller OTA 6/42)
    Gain = (gm1*gm7)/(g1*g7); 
    Gain_dB = 20*np.log10(Gain) #trouver formule du gain
    Poled = (g1*g7)/(gm7*Cf)
    Polend = gm7 *(Cf/(C1*CL + (C1+CL)*Cf))
    Zeron = gm7/Cf 

    # # ---------------------------------------------------------------------------
    # # ----------------- Bode plots - Miller -------------------------------------
    # # ---------------------------------------------------------------------------

    plot = True

    # Transfer function
    fmin = 1e0; fmax = 1e10
    Nw = 1000
    w = np.logspace(np.log10(2*np.pi*fmin), np.log10(2*np.pi*fmax), Nw)
    f = w/2/np.pi

    num = Gain * (1 + 1j * w / Zeron)
    den = (1 + 1j * w / Poled) * (1 + 1j * w / Polend) 

    # Open-loop
    H = num / den

    # Closed-loop
    G = H / (1 + H)

    def get_tf_angle(tf):
        angle = np.angle(tf,deg=True)
        for i in range(len(angle)):
            if angle[i] >= 0:
                angle[i] -= 360
        return angle

    mag_H = 20 * np.log10(abs(H))
    mag_G = 20 * np.log10(abs(G))
    angle_H = get_tf_angle(H)
    angle_G = get_tf_angle(G)

    min_dB = -120; max_dB = 120
    min_angle = -360; max_angle = 180



    # # ---------------------------------------------------------------------------
    # # -------------------- Design algorithm - Boost -----------------------------
    # # ---------------------------------------------------------------------------


    #recuperation from Miller ( eldo pour Avin)
    Avin=62
    alpha=Gain_dB/Avin #poru avoir la differnece entre gain prevu par python et elui trouvé dans eldo
    #print(alpha)
    Iin=ibias #puisque B=1 dans MIller le courant de sortie c'est toujorus Ibias

    #formulas dimensions
    B1=Avout/Avin #va etre le rapport de dimensions du miroir de courant du boost
    id17=Iin
    id16=id17
    idn17 = float(scint.interp1d(M17['GMID'], M17['IN'])(gmid17))
    idn16 = float(scint.interp1d(M16['GMID'], M16['IN'])(gmid16))
    W_over_L17 = id17/idn17
    W_over_L16 = id16/idn16
    W_over_L10 = B1*W_over_L17
    W17=W_over_L17*L17
    W16=W_over_L16*L16
    W10=W_over_L10*L10


    #graphs

    gm16 = gmid16*id16
    Gain_b=Gain/(B1*gm16) #POURQUOI DIVISER ET PAS MULTIPLIER??? PCQ COURANT?
    Gain_dB_b = 20*np.log10(Gain_b)
    Gain_dB_b=Gain_dB_b/alpha

    num_b = num/(B1*gm16) #AGAIN PQ?

    # Open-loop
    H_b = num_b / den

    # Closed-loop
    G_b = H_b / (1 + H_b)

    mag_H_b = 20 * np.log10(abs(H_b))
    mag_G_b = 20 * np.log10(abs(G_b))
    angle_H_b = get_tf_angle(H_b)
    angle_G_b = get_tf_angle(G_b)



    #for the write in file
    gmid10=B1*gmid17

    #apparemment besoin de ca ( pour dans txt poru saturation ou non )
    vdsat1 = float(scint.interp1d(M1['GMID'], M1['VDSAT'])(gmid1))
    vdsat2 = float(scint.interp1d(M2['GMID'], M2['VDSAT'])(gmid2))
    vdsat3 = float(scint.interp1d(M3['GMID'], M3['VDSAT'])(gmid3))
    vdsat4 = float(scint.interp1d(M4['GMID'], M4['VDSAT'])(gmid4))
    vdsat6 = float(scint.interp1d(M6['GMID'], M6['VDSAT'])(gmid6))
    vdsat5 = float(scint.interp1d(M5['GMID'], M5['VDSAT'])(gmid5))
    vdsat7 = float(scint.interp1d(M7['GMID'], M7['VDSAT'])(gmid7))
    vdsat8 = float(scint.interp1d(M8['GMID'], M8['VDSAT'])(gmid8))
    
    
    
    vdsat16 = float(scint.interp1d(M16['GMID'], M16['VDSAT'])(gmid16))
    vdsat17 = float(scint.interp1d(M17['GMID'], M17['VDSAT'])(gmid17))
    

    
    VIN = (1.2 - np.abs(vdsat5) - np.abs(vsg1) ) /2.0

    vdsat_list = [vdsat1, vdsat2, vdsat3, vdsat4, vdsat5, vdsat6, vdsat7, vdsat8, vdsat16, vdsat17]
    Flag_Correct_Saturation = transitor_are_satured(vdsat_list, 1.2)

    # # ---------------------------------------------------------------------------
    # # -------------------- Design algorithm - FF -----------------------------
    # # ---------------------------------------------------------------------------

    #DES AUTRES ETAGES
    #B1017=B1
    #id17=ibias
    #id8=ibias
    #dim M8, M3
    #gm1 = g1
    #g2 = gm7
    #gm16
    #Cl
    #id3

    #gain separes
    g1=gm1
    g2 = gm7
    B1017=B1
    gboost=B1017*gm16
    gff = 0.5*g1*B314*B1213 #sinon PM<45

    #M12
    idn12 = float(scint.interp1d(M12['GMID'], M12['IN'])(gmid12))
    id12= 2**id17*gff/g1
    W_over_L12 = id12/idn12
    W12 = W_over_L12*L12

    #M13
    W_over_L13 = B1213*W_over_L12
    W13 = W_over_L13 *L13

    #M14
    W_over_L14 = W_over_L3 *B314
    W14 = W_over_L14*L14

    #M15
    id14 = id3*B314
    id13 = id12*B1213
    id15 = id13+id14
    B158 = id15/id8
    W_over_L15 = B158*W8_over_L8 
    W15 = W_over_L15*L15

    #M11
    id10 = B1017*id17
    id11 = id12-id10
    B811 = id11/id8
    W_over_L11 = B811*W8_over_L8 
    W11 = W_over_L11 *L11


    #les capas
    g3=gboost
    Cc1 = (CL*g1)/(gff-g1) #c'est censé etre égal à Cf
    Cc2=(CL*g2*g3)/((gff-g1)**2)


    omegau = 2*np.pi*fT
    a = (Cc2 *( g2*gff - g1*g3 )) / (g1*g2*g3)
    b = Cc1*Cc2*( (gff-g1) / (g1*g2*g3) )
    c = Cc2 * ( ( Cc1*g2 - Cc1*g3 + CL*g2) / (Cc1*g2*g3) )
    d = (Cc2*CL) / (g2*g3)
    omega0_0 = (-np.sqrt(a**2 - 4*b) -a) / (2*b)
    omega0_1 = (np.sqrt(a**2 - 4*b) -a) / (2*b)
    omega0_2 = (-np.sqrt(c**2 - 4*d) -c) / (2*d)
    omega0_3 = (np.sqrt(c**2 - 4*d) -c) / (2*d)
    array = [omega0_0,omega0_1,omega0_2,omega0_3]
    
    if a**2 - 4*b < 0:
        Flag_Mathematical_integrity = 0
    if c**2 - 4*d < 0:
        Flag_Mathematical_integrity = 0

    if omega0_0 == omega0_1 == omega0_2 == omega0_3:
        omega0 = omega0_0
        #print("good")
    else:
        omega0 = np.mean(array)
        #print("mean")            
    #
    #Verif par rapport à l'article
    #

    verif=np.zeros(4)
    #condition sur Cc1 = Cf? 
    if np.abs((Cf -Cc1) / Cf) < 0.01:
        verif[0]=1.0
    #condition sur Cc2 = Cc1 = Cl/4
    if np.abs((Cc2 -Cf) / Cf) < 0.01:
        verif[1]=10.0
    #verif omu=om0
    if np.abs((omegau-omega0)/omega0) < 0.01:
        verif[2]=100.0
    #verif valeur om0
    if np.abs((omega0-(gff-g1)/CL) / ((gff-g1)/CL) ) < 0.01 :
        verif[3]=1000.0
    
    Flag_Correct_Design = verif.sum()
    
        
    # # ---------------------------------------------------------------------------
    # # --------------------------- Bode plots ------------------------------------
    # # ---------------------------------------------------------------------------
    
    #Miller    
    """   
    if plot:
        
        fig,ax = plt.subplots(2, 1, figsize=(20,10),sharex=True)
        ax[0].semilogx(f, mag_H,label='Open loop',linewidth=3)
        ax[1].semilogx(f, angle_H,linewidth=3)
        ax[0].semilogx(f, mag_G,label='Closed loop',linewidth=3)
        ax[1].semilogx(f, angle_G,linewidth=3)
        ax[0].legend(fontsize=16,framealpha=1,fancybox=False,edgecolor='k')
        ax[1].grid(True, which='both')
        ax[0].grid(True, which='both')
        ax[1].set_xlim((fmin,fmax))
        ax[1].set_xlabel("Frequency [Hz]",fontsize=20)
        ax[0].set_ylabel("Magnitude [dB]",fontsize=20)
        ax[1].set_ylabel("Phase [°]",fontsize=20)
        ax[1].set_ylim((min_angle,max_angle))
        ax[0].set_ylim((min_dB,max_dB))
        plt.title("Bode Miller")
        plt.show()

    arg_Pm = np.argmin(abs(mag_H - 0))
    Pm = angle_H[arg_Pm] + 180
    f_Pm = float(scint.interp1d(mag_H, f)(0))

    arg_Gm = np.argmin(abs(angle_H + 180))
    Gm = - mag_H[arg_Gm]
    #f_Gm = float(scint.interp1d(angle_H, f)(-120))

    if plot:

        fig,ax = plt.subplots(2, 1, figsize=(20,10),sharex=True)
        ax[0].semilogx(f, mag_H,label='Open loop',linewidth=3)
        ax[1].semilogx(f, angle_H,linewidth=3)
        ax[0].legend(fontsize=16,framealpha=1,fancybox=False,edgecolor='k')
        ax[0].grid(True, which='both')
        ax[1].grid(True, which='both')
        ax[1].set_xlim((fmin,fmax))
        ax[1].set_xlabel("Frequency [Hz]",fontsize=20)
        ax[0].set_ylabel("Magnitude [dB]",fontsize=20)
        ax[1].set_ylabel("Phase [°]",fontsize=20)
        ax[1].set_ylim((min_angle,max_angle))
        ax[0].set_ylim((min_dB,max_dB))
        
        ax[0].hlines(0, fmin, fmax, color='grey', linestyles='dashed')
        ax[1].hlines(-180, fmin, fmax, color='grey',linestyles='dashed')
        ax[0].vlines(f_Pm,min_dB,0,color='grey',linestyles='dashed')
        ax[1].vlines(f_Pm,Pm-180,max_angle,color='grey',linestyles='dashed')
        ax[1].vlines(f_Pm,-180,Pm-180,color='k')
    #   ax[1].vlines(f_Gm,-180,max_angle,color='grey',linestyles='dashed')
    #  ax[0].vlines(f_Gm,min_dB,-Gm,color='grey',linestyles='dashed')
    # ax[0].vlines(f_Gm,-Gm,0,color='k')
        #ax[0].set_title(r"$G_m =${:2.2f} dB at {:3.1f} MHz     $\Phi_m =${:2.2f}° at {:3.1f} kHz".format(Gm,f_Gm/settings.mega,Pm,f_Pm/settings.kilo),fontsize=20)
        plt.title("Miller open")
        plt.show()

    #Boost
    
    min_dB = -100; max_dB = 150
    min_angle = -360; max_angle = 180
    plot=True
    if plot:
        
        fig,ax = plt.subplots(2, 1, figsize=(20,10),sharex=True)
        ax[0].semilogx(f, mag_H_b,label='Open loop',linewidth=3)
        ax[1].semilogx(f, angle_H_b,linewidth=3)
        ax[0].semilogx(f, mag_G_b,label='Closed loop',linewidth=3)
        ax[1].semilogx(f, angle_G_b,linewidth=3)
        ax[0].legend(fontsize=16,framealpha=1,fancybox=False,edgecolor='k')
        ax[1].grid(True, which='both')
        ax[0].grid(True, which='both')
        ax[1].set_xlim((fmin,fmax))
        ax[1].set_xlabel("Frequency [Hz]",fontsize=20)
        ax[0].set_ylabel("Magnitude [dB]",fontsize=20)
        ax[1].set_ylabel("Phase [°]",fontsize=20)
        ax[1].set_ylim((min_angle,max_angle))
        ax[0].set_ylim((min_dB,max_dB))
        plt.title("Bode Miller + ampli")
        plt.show()   

    """
    W_array = [W1, W2, W3, W4, W5, W6, W7, W8, W10, W11, W12, W13, W14, W15, W16, W17]
    L_array = [L1, L2, L3, L4, L5, L6, L7, L8, L10, L11, L12, L13, L14, L15, L16, L17]
    
    Flag_result_is_physical = Result_is_physical(Cc1, Cc2, CL, ibias, VIN, W_array)
    Flags = [Flag_Correct_Saturation, Flag_Correct_Design, Flag_result_is_physical, Flag_Mathematical_integrity]
    
    
    return  Cf, Cc1, Cc2, CL, ibias, VIN, W_array, L_array, Flags, omega0, omegau


if __name__ == "__main__":
    
    fT_design = 13.7e06
    CL_design = 10e-12 
    #phase margin>45
    Avout_design=60 

    # # ---------------------------------------------------------------------------
    # # ----------------------------- Design choices ------------------------------
    # # ---------------------------------------------------------------------------


    #MILLER
    B=1 #gmid5,8 sont les memes

    gmid1_list = [13,14,15,16]
    gmid3_list =  [11,12,13,14]
    gmid5_list = [8,9,10,11]
    gmid6_list = [8,9,10,11] 

    #BOOST
    gmid16_list= [5,7,8,10,12,13,14,15,16,17,18,19,20,21,22]
    gmid17_list= [5,7,8,10,12,13,14,15,16,17,18,19,20,21,22]

    #FF
    gmid12_list= [5,7,8,10,12,13,14,15,16,17,18,19, 20,21,22]
    B1213_list = [10]
    
    B314=1
    #M13 et M14 mm 
    
    Cf_design = 1e-9

    file_result = open("result.txt", "w")
    file_result.write("gmid16 gmid17 gmid12 B12-13 Flag_Correct_Saturation Flag_Correct_Design Flag_result_is_physical Flag_result_is_mathematical omega0 omegau Cf Cc1 Cc2\n")
    
    import time
    t0 = time.time()
    
    iter_1 = 0.0
    N1 = len(gmid16_list) 
    N2 = len(gmid17_list)
    for gmid16 in (gmid16_list) : 
        iter_2 = 0.0
        for gmid17 in (gmid17_list) :
            iter_2 += 1.0
            print(f'progression : { (iter_1*N2 + iter_2) /(N1*N2) *100}%')
            for gmid12 in (gmid12_list) : 
                for B1213 in (B1213_list) : 
                    for gmid1 in (gmid1_list) :
                        for gmid3 in (gmid3_list) :
                            for gmid5 in (gmid5_list) :
                                for gmid6 in (gmid6_list) :
                                    Cf, Cc1, Cc2, CL, ibias, VIN, W_array, L_array, Flags, omega0, omegau = Miller_Sizing (fT_design, CL_design, Cf_design, Avout_design, gmid1, gmid3, gmid5, gmid6, gmid12, gmid16, gmid17, B, B1213, B314)

                                    file_result.write("{} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(gmid16, gmid17, gmid12, B1213, Flags[0] , Flags[1], Flags[2], Flags[3], omega0, omegau, Cf, Cc1, Cc2 ))
        iter_1 += 1.0
        
        
    t1 = time.time()
    print("Time elapsed: ", t1-t0)
    
    print("Done")
            