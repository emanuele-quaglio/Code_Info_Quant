# -*- coding: utf-8 -*-
"""
porter thomas plots
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.optimize import curve_fit
import scipy.stats as stats
from math import exp, log, floor, ceil, sqrt
#_____________________________________________
import matplotlib.font_manager as font_manager
mpl.rcParams['font.family']='serif'
cmfont = font_manager.FontProperties(fname=mpl.get_data_path() + '/fonts/ttf/cmr10.ttf')
mpl.rcParams['font.serif']=cmfont.get_name()
mpl.rcParams['mathtext.fontset']='cm'
mpl.rcParams['axes.unicode_minus']=False
#___________________________________________

def compat(x, sx,y, sy):
    return abs(x-y)/sqrt(sx*sx+sy*sy)
def vectorizedint(x):
    for i in range(len(x)):
        x[i]=int(x[i])
    return x
def vectorizedfloat(x):
    for i in range(len(x)):
        x[i]=float(x[i])
    return x
def PORTERTHOMAS(filename, freqs="teo"):
    data=[]
    with open(filename, 'r') as file:
        file.readline()
        while True:
            line = file.readline().split()
            if not line:
                break
            Q=line[0]
            M=line[1]
            C=line[2]
            I=line[3]
            if((M,C)!=('500','60') or I not in ['3', '6', '10'] ): #
                continue
            if freqs=="teo":
                line=line[8:-(len(line)-8)//2]#TEORICHE! riga che seleziona le freqs teoriche da un file le cui righe sono 'par1 par2 ...parK freqs_teo.... freqs_emp (stesso numero di elementi di freqs teo)'
            elif freqs=="emp":
                line=line[-(len(line)-8)//2:]#EMPIRICHE! riga che seleziona le freqs empiriche da un file le cui righe sono 'par1 par2 ...parK freqs_teo.... freqs_emp (stesso numero di elementi di freqs teo)'

            line=vectorizedfloat(line)
            line=sorted(line, reverse=True)
            
            '''
            line=[freq for freq in line if freq >0]
            my_norm=sum(np.array(line))
            line=[freq/my_norm for freq in line]
            '''
            for mybins in [11]: #range(2,30,1): 
                plt.yscale('log')
                n, bins, patchess=plt.hist(line, bins=mybins, density=True, edgecolor='gray', color='lightgray')
                plt.plot()
                x=np.linspace(0, bins[-1], 100)
                #plt.plot(x, PT_distr(x, len(line),  bins, normalized=False), label='Porter-Thomas')
                plt.plot(x, PT_distr(x, len(line),  bins, normalized=True), label='Porter-Thomas', color='blue', linestyle='--')
                #popt, pcov=curve_fit(decay_distr, bins[:len(bins)-1] , n, bounds=([0, 0],[np.inf, np.inf]))
                
                y_fit=[log(yi) for yi in n if yi>0]
                x_fit=[bins[i] for i in range(len(n)) if n[i]>0] 
                popt, pcov=curve_fit(lin_distr, x_fit , y_fit)
                
                a,b=exp(popt[1]), -popt[0]
                sa,sb= exp(popt[1])*sqrt(pcov[1][1]), sqrt(pcov[0][0])
                #fitlabel='fit: '+str(round(a))+'*e^(-'+str(round(b))+'x)'
                #fitlabel='fit: $C_1=$'+str(round(a))+'$\pm$'+str(round(sqrt(pcov[0][0])))+', $C_2=$'+str(round(10*b)/10)+'$\pm$'+str(round(10*sqrt(pcov[1][1]))/10)
                fitlabel='fit: $h=C_1 e^{-C_2 q}$'
                #plt.plot(x, decay_distr(x, *popt),  label=fitlabel)
                plt.plot(x, np.exp(lin_distr(x, *popt)),  label=fitlabel, color='black')
                #print(popt)
                
                print('a=',round(a),round(sa), 'b=', round(b),round(sb))
                print('compat C1='+str(compat(a,sa,512,0)), 'compat C2='+str(compat(b,sb,512,0)))
                print('\('+str(round(a))+' \pm '+str(round(sa))+' \) & '+str(compat(a,sa,512,0))+' & \('+str(round(b))+' \pm '+str(round(sb))+'\) & '+str(compat(b,sb,512,0))+' \\')
                            #plt.title('3-3-2, sample '+str(sample_number))
                #plt.title('$n=8, c=$'+str(C)+'$, I=$'+str(I)+', cube')
                plt.xlabel(('p','q')[freqs=='teo'], fontsize=12) #+', '+str(len(bins)-1)+' bins'
                plt.ylabel('h', fontsize=12)
                plt.legend(loc="upper right")
                pdftitle='grafici_finali/PTq8safety_I'+I+'.pdf'
                #plt.savefig(pdftitle)
                plt.show()
def lin_distr(x, a, b):
    x=np.array(x)
    return a*x+b
def exp_distr(x, a, b, c):
    x=np.array(x)
    return a*np.exp(b*x)+c
def decay_distr(x, a, b):
    x=np.array(x)
    return a*np.exp(-b*x)
def PT_distr(x, N, bins, normalized=False):
    histo=decay_distr(x, N, N)
    if normalized==True:
        temphisto=decay_distr(bins, N, N)
        sum=0
        for i in range(len(bins)-1):
           sum+=(bins[i+1]-bins[i])*temphisto[i] 
        histo/=sum
    return histo

PORTERTHOMAS("AT_q3x3_vs_l.txt", freqs="teo")