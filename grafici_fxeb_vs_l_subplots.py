# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 22:11:25 2023

@author: Emanuele
"""

# -*- coding: utf-8 -*-
"""
@author: Emanuele
"""
# -*- coding: utf-8 -*-
"""
@author: Emanuele
"""
# -*- coding: utf-8 -*-
#SCRIPT GRAFICI TESI TRIENNALE INFO QUANT
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.optimize import curve_fit
import scipy.stats as stats
from math import exp, log, floor, ceil, sqrt
from matplotlib.ticker import ScalarFormatter
import warnings
from matplotlib.mathtext import MathTextWarning
from matplotlib.ticker import FuncFormatter
#plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
#mpl.rcParams['font.family']='serif'
#mpl.rcParams['mathtext.fontset']='cm'
#__________________________________________________
import matplotlib.font_manager as font_manager
mpl.rcParams['font.family']='serif'
cmfont = font_manager.FontProperties(fname=mpl.get_data_path() + '/fonts/ttf/cmr10.ttf')
mpl.rcParams['font.serif']=cmfont.get_name()
mpl.rcParams['mathtext.fontset']='cm'
mpl.rcParams['axes.unicode_minus']=False

def dummy_symbol_fix(ax=None):
    if ax is None:
        ax = plt.gca()
    fig = ax.get_figure()
    # Force the figure to be drawn
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=MathTextWarning)
        fig.canvas.draw()
    # Remove '\mathdefault' from all minor tick labels
    labels = [label.get_text().replace('\mathdefault', '')
              for label in ax.get_xminorticklabels()]
    ax.set_xticklabels(labels, minor=True)
    labels = [label.get_text().replace('\mathdefault', '')
              for label in ax.get_yminorticklabels()]
    ax.set_yticklabels(labels, minor=True)
#_____________________________________________________


#funzione per analisi distribuzioni frequenze teoriche output circuito
#prende in input    -file con righe: 'q': line[0],'m':line[1],  'l': line[2], 'i':line[3], 'flin':line[4], 'sflin':line[5], 'flog':line[6], 'sflog':line[7], 'freqs': frequenze riordinate in ordine decrescente
#                   -funzione fittante distr(x, a,b,c,...)
#                   -f_x ed f_y funzioni da applicare ai dati x e y prima del fit (eventualmente identità)
def ramp(x, offset):
    if x>offset:
        return x*0.82
    else:
        return x
def vectorizedint(x):
    for i in range(len(x)):
        x[i]=int(x[i])
    return x
def vectorizedfloat(x):
    for i in range(len(x)):
        x[i]=float(x[i])
    return x
def data_organizer(filename, ftype1='flin', ftype2='flintt', ftype3='flog', ftype4='flogtt',  general_label=''):
    data=[]
    with open(filename, 'r') as file:
        file.readline()
        while True:
            line = file.readline().split()
            if not line: 
                break
            line[:4], line[4:]=vectorizedint(line[:4]), vectorizedfloat(line[4:])
            data.append({'q': line[0],'m':line[1],  'l':line[2], 'i':line[3], ftype1:line[4], ftype2:line[5], ftype3:line[6], ftype4:line[7], 'fs_t': line[8:-(len(line)-8)//2], 'fs_e': line[-(len(line)-8)//2:]})
    return data, (ftype1, ftype2,ftype3,ftype4), general_label
def data_averager(data, ftypes,  general_label=''):
    ftype1,ftype2,ftype3,ftype4=ftypes[0],ftypes[1],ftypes[2],ftypes[3]
    pset=[]
    for row in data:
        if list(row.items())[:3] not in pset:
            pset.append(list(row.items())[:3])
    av_data=[]
    for p in pset:
        analogs_ftype1=[row[ftype1] for row in data if list(row.items())[:3]==p]
        analogs_ftype2=[row[ftype2] for row in data if list(row.items())[:3]==p]
        analogs_ftype3=[row[ftype3] for row in data if list(row.items())[:3]==p]
        analogs_ftype4=[row[ftype4] for row in data if list(row.items())[:3]==p]
        av_ftype1, av_ftype2, av_ftype3, av_ftype4=sum(analogs_ftype1)/len(analogs_ftype1), sum(analogs_ftype2)/len(analogs_ftype2), sum(analogs_ftype3)/len(analogs_ftype3),sum(analogs_ftype4)/len(analogs_ftype4)
        s_ftype1, s_ftype2, s_ftype3, s_ftype4=sqrt(sum([pow(obj-av_ftype1, 2) for obj in analogs_ftype1])/len(analogs_ftype1)),sqrt(sum([pow(obj-av_ftype2, 2) for obj in analogs_ftype2])/len(analogs_ftype2)), sqrt(sum([pow(obj-av_ftype3, 2) for obj in analogs_ftype3])/len(analogs_ftype3)), sqrt(sum([pow(obj-av_ftype4, 2) for obj in analogs_ftype4])/len(analogs_ftype4))
        newstuffdict={'av'+ftype1 : av_ftype1, 's'+ftype1 : s_ftype1, 'av'+ftype2 : av_ftype2, 's'+ftype2 : s_ftype2, 'av'+ftype3 : av_ftype3, 's'+ftype3 : s_ftype3,'av'+ftype4 : av_ftype4, 's'+ftype4 : s_ftype4 }
        betterdict=dict(p)
        betterdict.update(newstuffdict)
        av_data.append(betterdict)
    return(av_data)
def data_displ_averager(data, ftypes,  general_label=''):
    ftype1,ftype2,ftype3,ftype4=ftypes[0],ftypes[1],ftypes[2],ftypes[3]
    pset=[]
    for row in data:
        if list(row.items())[:3] not in pset:
            pset.append(list(row.items())[:3])
    av_data=[]
    for p in pset:
        analogs_ftype1=[row[ftype1] for row in data if list(row.items())[:3]==p]
        analogs_ftype2=[row[ftype2] for row in data if list(row.items())[:3]==p]
        analogs_ftype3=[row[ftype3] for row in data if list(row.items())[:3]==p]
        analogs_ftype4=[row[ftype4] for row in data if list(row.items())[:3]==p]
        #print(analogs_ftype1, analogs_ftype2, analogs_ftype3, analogs_ftype4)
        #sostituisco flog con abs(flog-flogtt)/flogtt
        analogs_ftype3=[abs(analogs_ftype3[i]-analogs_ftype4[i])/analogs_ftype4[i] for i in range(len(analogs_ftype3)) ]
        #print(analogs_ftype1, analogs_ftype2, analogs_ftype3, analogs_ftype4)
        av_ftype1, av_ftype2, av_ftype3, av_ftype4=sum(analogs_ftype1)/len(analogs_ftype1), sum(analogs_ftype2)/len(analogs_ftype2), sum(analogs_ftype3)/len(analogs_ftype3),sum(analogs_ftype4)/len(analogs_ftype4)
        #print('averages: ',av_ftype1, av_ftype2, av_ftype3, av_ftype4)
        s_ftype1, s_ftype2, s_ftype3, s_ftype4=sqrt(sum([pow(obj-av_ftype1, 2) for obj in analogs_ftype1])/len(analogs_ftype1)),sqrt(sum([pow(obj-av_ftype2, 2) for obj in analogs_ftype2])/len(analogs_ftype2)), sqrt(sum([pow(obj-av_ftype3, 2) for obj in analogs_ftype3])/len(analogs_ftype3)), sqrt(sum([pow(obj-av_ftype4, 2) for obj in analogs_ftype4])/len(analogs_ftype4))
        newstuffdict={'av'+ftype1 : av_ftype1, 's'+ftype1 : s_ftype1, 'av'+ftype2 : av_ftype2, 's'+ftype2 : s_ftype2, 'av'+ftype3 : av_ftype3, 's'+ftype3 : s_ftype3,'av'+ftype4 : av_ftype4, 's'+ftype4 : s_ftype4 }
        print(newstuffdict)
        betterdict=dict(p)
        betterdict.update(newstuffdict)
        av_data.append(betterdict)
    return(av_data)
def lin_distr(x, a, b): 
    x=np.array(x)
    return a*x+b
def exp_distr(x, a, b, c):
    x=np.array(x)
    return a*np.exp(b*x)+c
def trim_edges(x, da=1/5, a=3/5):
    return x[floor(len(x)*da):ceil(len(x)*a)]
def customlog(x):
    x=np.array(x)
    nonzerox=[]
    for xi in x:
        if xi!=0:
            nonzerox.append(xi)        
    custom0=min(nonzerox)*1E-6
    for i in range(len(x)):
        if x[i]==0:
            x[i]=custom0
    return np.log(x)
def X2vsQLscatter(data, general_label=''):
    plt.scatter([row['q'] for row in data], [row['l'] for row in data], c=[row['X2'] for row in data])
    plt.colorbar()
    plt.title('colormap del χ2'+' '+general_label)
    plt.xlabel('n qubits')
    plt.ylabel('n cycles')
    #plt.show()
def X2vsMplot(data, general_label=''):
    fig, ax = plt.subplots(1,1)
    ax.set_xscale("log")
    ax.step([row['m'] for row in data if 'X2' in row], [row['X2'] for row in data if 'X2' in row])
    ax.set(title='χ2 vs # meas, q=8, c=27, i=10'+general_label)
    ax.set_xlabel('# meas')
    ax.set_ylabel('χ2')
    #plt.show()
def FXEBvsMplot(data, ftypes, general_label='', myax=plt.gca()):
    #myax.set_xscale('log')
    #myax.set_yscale('log')
    myax.loglog()
    print([data[i]['avflog'] for i in range(len(data))])
    if 'sflog' in data[0]:
        x=[row['m'] for row in data]
        y=[row['avflog'] for row in data]
        sigmay=[row['sflog'] for row in data]
        myax.fill_between(x, [abs(y[i]-sigmay[i]) for i in range(len(y))], [y[i]+sigmay[i] for i in range(len(y))], alpha=0.15, linewidth=1, edgecolor='black')
        myax.plot(x,y)
        #plt.errorbar(x,y, yerr=sigmay, fmt='.')
    else:
       iset=set()
       for row in data:
           iset.add(row['i'])
       for i in iset:
           plt.plot([row['m'] for row in data if row['i']==i], [row[ftype]-row['flin(t,t)'] for row in data if row['i']==i])
   
    #plt.title('$n=8, c=20, N_I=10$'+general_label, y=0.9)
    myax.axvline(x=500, linestyle='dashed', color='gray', label='$N_S=500$')
    
    #plt.show()
def FXEBvsLplot(data, ftypes, general_label='', ax=None):
    ax.set_xscale('log')
    ax.set_yscale('log')
    #plt.ylim(ymin=0.2)
    if 'sflin' in data[0]:
        #ax.errorbar([row['l'] for row in data], [row['avflin'] for row in data], yerr=[row['sflin'] for row in data], fmt='.')
        #plt.title('$n=8, N_S=500, N_I=10$'+general_label)
        
        x=[row['l'] for row in data]
        y=[row['avflin'] for row in data]
        sigmay=[row['sflin'] for row in data]
        ax.fill_between(x, [y[i]-sigmay[i] for i in range(len(y))], [y[i]+sigmay[i] for i in range(len(y))], alpha=0.5, linewidth=1, edgecolor='black', facecolor='orange')
        ax.plot(x,y, color='orangered')
    else:
        for i in range(1,11):
            plt.plot([row['l'] for row in data if row['i']==i], [row['flin'] for row in data if row['i']==i] ) #, marker='o'
            #plt.title('$f_{XEB, lin}$ vs # cycles, q=8 (332), m=500, i=1'+general_label)
        #plt.title('$n=8, N_S=500, N_I=10$'+general_label)
    
    
    ax.minorticks_on()
    #for axis in [ax.xaxis, ax.yaxis]:
        #axis.set_major_formatter(ScalarFormatter())
    #ax.set_scientific(False)
    ax.tick_params(which='minor') #, labelcolor='white', labelsize=0
    #ax.grid(which='minor', linestyle='--')
    #ax.grid(which='major', linestyle='--')
    ax.axhline(y=1, linestyle='dashed', color='black')
    #dummy_symbol_fix(ax)
    
    #plt.savefig('grafici_finali/safety.pdf')
    #plt.show()
    
def FXEBvsQplot(data, general_label=''):
    plt.errorbar([row['q'] for row in data if 'X2' in row], [row['flog'] for row in data if 'X2' in row], yerr=[row['sflog'] for row in data if 'X2' in row], fmt='.')
    plt.title('fXEB-log vs # qubits, c=27, m=10^4, i=10'+general_label)
    plt.xlabel('# qubits')
    plt.ylabel('fXEB-log')
    #plt.show()        

fig,axes=plt.subplots(2,2, sharey=True, sharex=True) #,figsize=(6.7,5.06)
types=['8','4x2', '332s', '3x3']
labels=['a','b','c','d']
for i in [0,1]:
    for j in [0,1]:
        mydata, myftypes, gl=data_organizer("AT_q"+types[2*i+j]+"_vs_l.txt", 'flin', 'flintt', 'flog', 'flogtt',
                                    ', linear chain')
        FXEBvsLplot(data_averager(mydata,myftypes ),myftypes, gl, ax=axes[i][j])
        axes[i][j].set_title('('+labels[2*i+j]+')', loc='right', y=0.83)
for i in [0,1]:
    axes[i][0].set_ylabel('$f_{XEB, lin}$', fontsize=14)
    axes[i][0].set_yticks([1, 10])
    axes[i][0].get_yaxis().set_major_formatter(ScalarFormatter())
for j in [0,1]:
    axes[1][j].set_xlabel('$c$', fontsize=15)
    axes[1][j].set_xticks([2, 10, 60])
    axes[1][j].get_xaxis().set_major_formatter(ScalarFormatter())

fig.tight_layout()
fig.savefig('grafici_finali/flin_vs_l_multiplot_logx_cleaner.pdf')
plt.show()
#FXEBvsLplot(mydata,myftypes, gl)
"""
md=data_averager(mydata,myftypes )
print(sum([row['avflin'] for row in md])/len(md))
print(sum([row['sflin'] for row in md])/len(md))
print(max([row['avflin'] for row in md][7:]))
print(max([row['sflin'] for row in md][7:]))
"""
"""
#roba per FXEBvsMplot
fig, ((ax1,ax2),(ax3,ax4))=plt.subplots(2,2,sharex=True, sharey=True, )
mydata, myftypes, gl=data_organizer("AT_q8_vs_m.txt", 'flin', 'flintt', 'flog', 'flogtt',
                                    ', linear chain')
FXEBvsMplot(data_displ_averager(mydata,myftypes),myftypes, gl, ax1)
mydata, myftypes, gl=data_organizer("AT_q4x2_vs_m.txt", 'flin', 'flintt', 'flog', 'flogtt',
                                    ', double chain')
FXEBvsMplot(data_displ_averager(mydata,myftypes),myftypes, gl, ax2)
mydata, myftypes, gl=data_organizer("AT_q332s_vs_m.txt", 'flin', 'flintt', 'flog', 'flogtt',
                                    ', L-tromino')
FXEBvsMplot(data_displ_averager(mydata,myftypes),myftypes, gl, ax3)
mydata, myftypes, gl=data_organizer("AT_q3x3_vs_m.txt", 'flin', 'flintt', 'flog', 'flogtt',
                                    ' ')
FXEBvsMplot(data_displ_averager(mydata,myftypes),myftypes, gl, ax4)

#plt.savefig('grafici_finali/q8_vs_m.pdf')
ax3.set_xlabel('$N_S$', fontsize=12)
ax4.set_xlabel('$N_S$', fontsize=12)
ax1.set_ylabel('$\epsilon$', fontsize=15)
ax3.set_ylabel('$\epsilon$', fontsize=15)

'''
ax1.set_title('(a)', x=0.1,y=0.04)
ax2.set_title('(b)', x=0.1,y=0.04)
ax3.set_title('(c)', x=0.1,y=0.04)
ax4.set_title('(d)', x=0.1,y=0.04)
ax2.legend(loc="upper right")
'''
ax1.set_title('(a)', loc='right', y=0.8)
ax2.set_title('(b)',  loc='right', y=0.8)
ax3.set_title('(c)',  loc='right', y=0.8)
ax4.set_title('(d)',  loc='right', y=0.8)
ax2.legend(loc="lower left")
fig.tight_layout()

fig.savefig('grafici_finali/safety.pdf')
plt.show()
"""
