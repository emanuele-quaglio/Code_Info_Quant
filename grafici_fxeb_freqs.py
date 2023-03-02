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
from math import exp, log, floor, ceil


#funzione per analisi distribuzioni frequenze teoriche output circuito
#prende in input    -file con righe: 'q': line[0],'m':line[1],  'l': line[2], 'i':line[3], 'flin':line[4], 'sflin':line[5], 'flog':line[6], 'sflog':line[7], 'freqs': frequenze riordinate in ordine decrescente
#                   -funzione fittante distr(x, a,b,c,...)
#                   -f_x ed f_y funzioni da applicare ai dati x e y prima del fit (eventualmente identità)
def vectorizedint(x):
    for i in range(len(x)):
        x[i]=int(x[i])
    return x
def vectorizedfloat(x):
    for i in range(len(x)):
        x[i]=float(x[i])
    return x
def analisi_distr_freqs_theo(filename, distr, f_x=lambda x:x, f_y=lambda y:y, f_xy=lambda xy:xy, f_xy_domain=False, my_ylabel='frequenze teoriche', title_add='', my_bounds=(-np.inf, np.inf)):
    data=[]
    badline=set()
    with open(filename, 'r') as file:
        min_chi=0
        count=0
        index_min_chi=0
        file.readline()
        while True:
            line = file.readline().split()
            if not line:
                break
            line[:4], line[4:]=vectorizedint(line[:4]), vectorizedfloat(line[4:])
            x, y=np.array(range(len(line)-8), dtype=np.int64), np.array(sorted(line[8:])[::-1], dtype=np.float64)
            x_fit, y_fit=f_x(f_xy(x)), f_y(f_xy(y))
            
            #plt.plot(f_x(range(len(y))), f_y(y), color='black')
            #plt.show()
            
            popt, pcov=curve_fit(distr, x_fit , y_fit, bounds=my_bounds)
            try:
                chi_square, p_value = stats.chisquare(y_fit, distr(x_fit, *tuple(popt)), ddof=len(popt))    
                data.append({'q': line[0],'m':line[1],  'l': line[2], 'i':line[3], 'flin':line[4], 'sflin':line[5], 'flog':line[6], 'sflog':line[7], 'freqs': y, 'param':tuple(popt), 'pcov': pcov, 'X2': chi_square, 'pval': p_value})
                if (count==0 or chi_square<min_chi):
                    min_chi=chi_square
                    index_min_chi=count
            except ValueError:
                badline.add(count)
                data.append({'q': line[0],'m':line[1],  'l': line[2], 'i':line[3], 'flin':line[4], 'sflin':line[5], 'flog':line[6], 'sflog':line[7], 'freqs': y, 'param':tuple(popt), 'pcov': pcov})
            count+=1
    #stampa un grafico della funzione fittata
    def fitted_distr(index):
        row=data[index]
        plt.title(' q: '+str(row['q'])+', l: '+str(row['l'])+', m: '+str(row['m'])+title_add)
        plt.xlabel('sorted bitstrings')
        plt.ylabel(my_ylabel)
        #condizione if che in base a un opzione di analisis distr applica f_xy anche a i dati normali nel plot, ad es per mostrarne un sottinsieme
        def cond_f_xy(x):
            if f_xy_domain==True:
                return f_xy(x)
            else:
                return x
        N=len(row['freqs'])
        plt.plot(f_x(cond_f_xy(range(N))), f_y(cond_f_xy(row['freqs'])), color='black')
        #print ("FREQS: ", row['freqs'])
        x_fit=f_x(cond_f_xy(np.linspace(0, N, 100)))
        plt.plot(x_fit, distr(x_fit, *row['param']))
        plt.show()
        print(index, 'q: ', data[index]['q'], 'm: ', data[index]['m'], ' l: ', data[index]['l'] )
        if 'X2' in data[index]:
            print(' X^2: ', data[index]['X2'])
        
        #print([row['X2'] for row in data])
    #stampo uno scatter a colori del chi^2 in funzione di n_qubits "q" e n_layers "l"
    #funzione che ritorna bool con condizioni generali per sottoselezione dei dati nei grafici
    def Condition(row):
        return row['q']==8
    
    def X2vsQLscatter():
        plt.scatter([row['q'] for row in data], [row['l'] for row in data], c=[row['X2'] for row in data])
        plt.colorbar()
        plt.title('colormap del χ2'+' '+title_add)
        plt.xlabel('n qubits')
        plt.ylabel('n cycles')
        plt.show()
    def X2vsMplot():
        fig, ax = plt.subplots(1,1)
        ax.set_xscale("log")
        ax.step([row['m'] for row in data if 'X2' in row], [row['X2'] for row in data if 'X2' in row])
        ax.set(title='χ2 vs # meas, q=8, c=27, i=10'+title_add)
        ax.set_xlabel('# meas')
        ax.set_ylabel('χ2')
        plt.show()
    def FXEBvsMplot():
        fig, ax = plt.subplots(1,1)
        ax.set_xscale("log")
        ax.errorbar([row['m'] for row in data if 'X2' in row], [row['flog'] for row in data if 'X2' in row], yerr=[row['sflog'] for row in data if 'X2' in row], fmt='.')
        ax.set(title='fXEB-log vs # meas, q=8, c=27, i=10'+title_add)
        ax.set_xlabel('# meas')
        ax.set_ylabel('fXEB-log')
        plt.show()
    def FXEBvsLplot():
        plt.errorbar([row['l'] for row in data if( 'X2' in row and Condition(row)) ], [(row['flin']+1)/103*256-1 for row in data if( 'X2' in row and Condition(row))  ], yerr=[row['sflin']*256/103 for row in data if 'X2' in row and Condition(row)], fmt='.')
        plt.title('fXEB-lin vs # cycles, q=8, m=10^4, i=30'+title_add)
        plt.xlabel('# cycles')
        plt.ylabel('fXEB-lin')
        plt.show()
    def FXEBvsQplot():
        plt.errorbar([row['q'] for row in data if 'X2' in row], [row['flog'] for row in data if 'X2' in row], yerr=[row['sflog'] for row in data if 'X2' in row], fmt='.')
        plt.title('fXEB-log vs # qubits, c=27, m=10^4, i=10'+title_add)
        plt.xlabel('# qubits')
        plt.ylabel('fXEB-log')
        plt.show()        
    
    #STAMPA DI DIVERSI GRAFICI, CHIAMANDO LE FUNZIONI QUI SOPRA DEFINITE:
    '''
    print(badline)
    print(len(data))
    '''
    
    for t in range(len(data)):
        fitted_distr(t)
    
    #X2vsMplot()
    #FXEBvsLplot();
    #FXEBvsMplot()
    """
    plt.plot([row['l'] for row in data if row['q']==4], [row['X2'] for row in data if row['q']==4] )
    plt.show()
    plt.plot([row['q'] for row in data if row['l']==12], [row['X2'] for row in data if row['l']==12] )
    plt.show()
    """
    #print(data)

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
        
#MAIN
#exp_param_bounds=([0, -np.inf, -np.inf],[np.inf, 0, np.inf])
#analisi_distr_freqs_theo('analisi_fxeb_q8_l1_50_trim.txt', exp_distr, my_bounds=exp_param_bounds, my_ylabel='measured frequencies', title_add=', subset: 0.2-0.6')
#analisi_distr_freqs_theo('analisi_fxeb_q8_l10204050.txt', lin_distr, f_y=lambda x:customlog(x), my_ylabel='frequenze misurate (log)')      
#analisi_distr_freqs_theo('analisi_fxeb_q8_l10204050.txt', lin_distr, f_y=lambda x:customlog(x), f_xy=lambda xy:trim_edges(xy, da=0.2, a=0.7), my_ylabel='frequenze misurate (log)', f_xy_domain=True, title_add='subset 0.2-0.7, zoom')      

