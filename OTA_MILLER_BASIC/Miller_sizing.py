
"""
Miller 1st order Python sizing

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

import codecs



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


M1 = pch_lvt; M2 = pch_lvt
M3 = nch_lvt; M4 = nch_lvt
M5 = pch_lvt; M6 = pch_lvt
M7 = nch_lvt; M8 = pch_lvt



global gm1,gm2,gm3,gm4,gm5,gm6,gm7,gm8
global g1,g2,g3,g4,g5,g6,g7,g8
global id1,id2,id3,id4,id5,id6,id7,id8
global idn1,idn2,idn3,idn4,idn5,idn6,idn7,idn8
global vsg1,vsg2,vgs3,vgs4,vsg5,vsg8,vsg6,vgs7

global VIN

def Miller_sizing(Cf, fT, CL, L, gmid1, gmid3, gmid5, gmid6, gmid7, gmid8):
    
    L1 = L; L2 = L1;
    L3 = L; L4 = L3; 
    L5 = L; L6 = L5;
    L7 = L; 
    L8 = L;
    B=1;
    
    gmid2 = gmid1; gmid4 = gmid3;
    
    
    Cfavant = Cf /1000.0
    i = 0
    
    Cfchoix = Cf
    Cfavantchoix = Cfavant

    while np.abs((Cf-Cfavant)/Cf) > 0.000001:
        i+=1
        Cfavant=Cf
    # # ---------------------------------------------------------------------------
    # # ----------------------------- Design algorithm ----------------------------
    # # ---------------------------------------------------------------------------

    # M1, M2, + cascode M3 M4
        gm1 = 2*np.pi * fT * Cf
        id1 = gm1 /gmid1
        idn1 = float(scint.interp1d(M1['GMID'], M1['IN'])(gmid1))
        vsg1 = float(scint.interp1d(M1['GMID'], M1['VGS'])(gmid1))
        W_over_L1 = id1 / idn1
        W1 = W_over_L1 * L1
        W2 = W1
        id2 = id1
        vsg2 = vsg1


        id3= id1
        id4= id1
        idn3 = float(scint.interp1d(M3['GMID'], M3['IN'])(gmid3))
        vgs3 = float(scint.interp1d(M3['GMID'], M3['VGS'])(gmid3))
        vgs4=vgs3
        W_over_L3 = id3 / idn3
        W3 = W_over_L3 * L3
        W4 = W3


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


    #transi de Miller M7
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





    

    #apparemment besoin de ca
        vdsat1 = float(scint.interp1d(M1['GMID'], M1['VDSAT'])(gmid1))
        vdsat2 = float(scint.interp1d(M2['GMID'], M2['VDSAT'])(gmid2))
        vdsat3 = float(scint.interp1d(M3['GMID'], M3['VDSAT'])(gmid3))
        vdsat4 = float(scint.interp1d(M4['GMID'], M4['VDSAT'])(gmid4))
        vdsat6 = float(scint.interp1d(M6['GMID'], M6['VDSAT'])(gmid6))
        vdsat5 = float(scint.interp1d(M5['GMID'], M5['VDSAT'])(gmid5))
        vdsat7 = float(scint.interp1d(M7['GMID'], M7['VDSAT'])(gmid7))
        vdsat8 = float(scint.interp1d(M8['GMID'], M8['VDSAT'])(gmid8))
    #    vdsat9 = float(scint.interp1d(M9['GMID'], M9['VDSAT'])(gmid9))

        vea1 = float(scint.interp1d(M1['GMID'], M1['VEA'])(gmid1))
        vea7 = float(scint.interp1d(M7['GMID'], M7['VEA'])(gmid7)) 

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
        
        VIN = (1.2 - np.abs(vdsat5) - np.abs(vsg1) ) /2.0

    #A REDEFINIR!!!
    Gain = (gm1*gm7)/(g1*g7); Gain_dB = 20*np.log10(Gain) #trouver formule du gain
    #Poled = g1/C1
    #Polend = g7/CL
    Poled = (g1*g7)/(gm7*Cf)
    Polend = gm7 *(Cf/(C1*CL + (C1+CL)*Cf))
    Zeron = gm7/Cf 

    # # ---------------------------------------------------------------------------
    # # --------------------------- Bode plots -------------------------------------
    # # ---------------------------------------------------------------------------

    plot = False

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
        plt.show()


    # # ---------------------------------------------------------------------------
    # # --------------------- Write results to file -------------------------------
    # # ---------------------------------------------------------------------------

    filename = "Miller_design.txt"
    encoding = "utf-8"
    with open(filename, "w", encoding=encoding) as file:
        file.write("Choix:\n")
        file.write("\tFt: {:4.1f} [MHz]\n".format(fT/settings.mega))
        file.write("\tCL: {:4.1f} [pF]\n".format(CL/settings.pico))
        file.write("\tgmid1: {:4.1f} \n".format(gmid1))
        file.write("\tgmid2: {:4.1f} \n".format(gmid2))
        file.write("\tgmid3: {:4.1f} \n".format(gmid3))
        file.write("\tgmid4: {:4.1f} \n".format(gmid4))
        file.write("\tgmid5: {:4.1f} \n".format(gmid5))
        file.write("\tgmid6: {:4.1f} \n".format(gmid6))
        file.write("\tgmid7: {:4.1f} \n".format(gmid7))
        file.write("\tgmid8: {:4.1f} \n".format(gmid8))
        file.write("\tB: {:4.1f} \n".format(B))
        file.write("\tCf {:3.2f} [pF]\n".format(Cfchoix/settings.pico))
        file.write("\tCfavant {:3.2f} [pF]\n".format(Cfavantchoix/settings.pico))
        file.write("Longueurs voir dimensions et tous les transis sont des n.pch_lvt.\n \n")
        
        file.write("Performance:\n")
        file.write("\tGain: {:3.2f} [dB]\n".format(Gain_dB))
        file.write("\tTransition frequency: {:4.1f} [MHz]\n".format(f_Pm/settings.mega))
        file.write("\tDominant pole: {:3.2f} [kHz]\n".format(Poled/2/np.pi/settings.kilo))

        #file.write("\tNon-dominant pole 1: {:3.2f} [MHz]\n".format(Polep/2/np.pi/settings.mega))
        file.write("\tNon-dominant pole 2: {:3.2f} [MHz]\n".format(Polend/2/np.pi/settings.mega))
        file.write("\tPhase margin: {:2.2f} [°]\n".format(Pm))
        file.write("\tBias current: {:3.2f} [µA]\n".format(ibias/settings.micro))
        file.write("\n")
        file.write("DC point:\n")
        file.write("M1:\t\tVSG: {:1.2f} [V]\tVDSsat: {:1.2f} [V]\n".format(vsg1,abs(vdsat1)))
        file.write("M2:\t\tVSG: {:1.2f} [V]\tVDSsat: {:1.2f} [V]\n".format(vsg2,abs(vdsat2)))
        file.write("M3:\t\tVGS: {:1.2f} [V]\tVSDsat: {:1.2f} [V]\n".format(vgs3,vdsat3))
        file.write("M4:\t\tVGS: {:1.2f} [V]\tVSDsat: {:1.2f} [V]\n".format(vgs4,vdsat4))
        file.write("M5:\t\tVSG: {:1.2f} [V]\tVSDsat: {:1.2f} [V]\n".format(vsg5,abs(vdsat5)))
        file.write("M6:\t\tVSG: {:1.2f} [V]\tVSDsat: {:1.2f} [V]\n".format(vsg6,abs(vdsat6)))
        file.write("M7:\t\tVGS: {:1.2f} [V]\tVDSsat: {:1.2f} [V]\n".format(vgs7,vdsat7))
        file.write("M8:\t\tVSG: {:1.2f} [V]\tVDSsat: {:1.2f} [V]\n".format(vsg8,abs(vdsat8)))
        #file.write("M9:\tVGS: {:1.2f} [V]\tVDSsat: {:1.2f} [V]\n".format(vgs9,vdsat9))
        file.write("\n")
        file.write("Dimensions:\n")
        file.write("\tW1: {:3.2f} [µm]\tL1: {:1.2f} [µm]\n".format(W1 / settings.micro,L1 / settings.micro))
        file.write("\tW2: {:3.2f} [µm]\tL2: {:1.2f} [µm]\n".format(W2 / settings.micro,L2 / settings.micro))
        file.write("\tW3: {:3.2f} [µm]\tL3: {:1.2f} [µm]\n".format(W3 / settings.micro,L3 / settings.micro))
        file.write("\tW4: {:3.2f} [µm]\tL4: {:1.2f} [µm]\n".format(W4 / settings.micro,L4 / settings.micro))
        file.write("\tW5: {:3.2f} [µm]\tL5: {:1.2f} [µm]\n".format(W5 / settings.micro,L5 / settings.micro))
        file.write("\tW6: {:3.2f} [µm]\tL6: {:1.2f} [µm]\n".format(W6 / settings.micro,L6 / settings.micro))
        file.write("\tW7: {:3.2f} [µm]\tL7: {:1.2f} [µm]\n".format(W7 / settings.micro,L7 / settings.micro))
        file.write("\tW8: {:3.2f} [µm]\tL8: {:1.2f} [µm]\n".format(W8 / settings.micro,L8 / settings.micro))
        file.write("\n")
        file.write("Capacités :\n")
        file.write("\tC1: {:3.2f} [pF]\n".format(C1 / settings.pico))
        file.write("\tCL: {:3.2f} [pF]\n".format(CL / settings.pico))
        file.write("\tCf: {:3.2f} [pF]\n".format(Cf / settings.pico))
        #file.write("\tW9: {:3.2f} [µm]\tL9: {:1.2f} [µm]\n".format(W9 / settings.micro,L9 / settings.micro))
        #file.write("\tW10: {:3.2f} [µm]\tL10: {:1.2f} [µm]\n".format(W10 / settings.micro,L10 / settings.micro))

        W_list =  np.array([W1, W2, W3, W4, W5, W6, W7, W8])/settings.micro
        L_list = np.array([L1, L2, L3, L4, L5, L6, L7, L8])/settings.micro
        vdsat_list = np.abs(np.array([vdsat1, vdsat2, vdsat3, vdsat4, vdsat5, vdsat6, vdsat7, vdsat8]))
        return Cf/settings.pico, CL/settings.pico, W_list, L_list, ibias/settings.micro, VIN, vdsat_list

if __name__ == '__main__':
        
    # # ---------------------------------------------------------------------------
    # # ----------------------------- Specifications ------------------------------
    # # ---------------------------------------------------------------------------

    fT = 13.7e06
    CL = 10e-12 

    # # ---------------------------------------------------------------------------
    # # ----------------------------- Design choices ------------------------------
    # # ---------------------------------------------------------------------------

    Cf = 3e-9
    Cfavant = 3e-12
    i = 0
    L = 1e-6



    gmid1= 15.0 
    gmid2= 15.0 
    gmid3= 13.0 
    gmid4= 13.0 
    gmid5= 10.0 
    gmid6= 10.0 
    gmid7= 13.0 
    gmid8= 10.0 
    
    
    Cf, CL, W_list, L_list, ibias, VIN, vdsat_list = Miller_sizing(Cf, fT, CL, L, gmid1, gmid3, gmid5, gmid6, gmid7, gmid8)


    #print W_list
    print("W_list = ", W_list)
    
    print(VIN)
    
    