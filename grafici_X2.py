# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
@author: Emanuele
"""
#SCRIPT GRAFICI TESI TRIENNALE INFO QUANT
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.optimize import curve_fit
import scipy.stats as stats
from math import exp, log, floor, ceil


#funzione per analisi distribuzioni frequenze teoriche output circuito
#prende in input    -file con righe: n_qubit, n_layers, frequenze teoriche riordinate in ordine decrescente
#                   -funzione fittante distr(x, a,b,c,...)
#                   -f_x ed f_y funzioni da applicare ai dati x e y prima del fit (eventualmente identità)
def analisi_distr_freqs_theo(filename, distr, f_x=lambda x:x, f_y=lambda y:y, my_ylabel='frequenze teoriche'):
    data=[]
    with open(filename, 'r') as file:
        min_chi=0
        count=0
        index_min_chi=0
        file.readline()
        while True:
            line = file.readline().split()
            if not line:
                break
            line[0], line[1], line[2:]=int(line[0]), int(line[1]), [float(word) for word in line[2:]]
            x, y=np.array(range(len(line)-2), dtype=np.int64), np.array(line[2:], dtype=np.float64)
            x_fit, y_fit=f_x(x), f_y(y)
            popt, pcov=curve_fit(distr, x_fit , y_fit)
            chi_square, p_value = stats.chisquare(5/min(y_fit)*y_fit, 5/min(y_fit)*distr(x_fit, *tuple(popt)), ddof=len(popt))
            data.append({'q': line[0], 'l': line[1], 'freqs': y, 'param':tuple(popt), 'pcov': pcov, 'X2': chi_square, 'pval': p_value})
            if (count==0 or chi_square<min_chi):
                min_chi=chi_square
                index_min_chi=count
            count+=1
            """
            if count==4:
                print(np.ndarray.tolist(y_fit))
            """
    #stampa un grafico della funzione fittata
    def fitted_distr(index):
        row=data[index]
        plt.title(' q: '+str(row['q'])+', l: '+str(row['l']))
        plt.xlabel('elementi di base riordinati')
        plt.ylabel(my_ylabel)
        N=len(row['freqs'])
        plt.plot(f_x(range(N)), f_y(row['freqs']), color='black')
        #print ("FREQS: ", row['freqs'])
        x_fit=f_x(np.linspace(0, N, 100))
        plt.plot(x_fit, distr(x_fit, *row['param']))
        plt.show()
        print(index, 'q: ', data[index]['q'], ' l: ', data[index]['l'], ' X^2: ', data[index]['X2'] )
        #print([row['X2'] for row in data])
    #stampo uno scatter a colori del chi^2 in funzione di n_qubits "q" e n_layers "l" 
    def X2vsQLscatter():
        plt.scatter([row['q'] for row in data], [row['l'] for row in data], c=[row['X2'] for row in data])
        plt.colorbar()
        plt.title('colormap del X^2')
        plt.xlabel('n qubits')
        plt.ylabel('n layers')
        plt.show()
    #stampo un grafico della distribuzione fittata
    X2vsQLscatter()
    for t in range(24):
        fitted_distr(t) 
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
def exp_distr(x, a, b):
    x=np.array(x)
    return a*np.exp(b*x)
def trim_edges(x):
    return x[floor(len(x)*1/5):ceil(len(x)*3/5)]
    
#MAIN
analisi_distr_freqs_theo('analisi_distr_freqs_theo_X2.txt', lin_distr, f_x=trim_edges, f_y=lambda x:trim_edges(np.log(x)))
#analisi_distr_freqs_theo('analisi_distr_freqs_theo_X2.txt', lin_distr, f_y=lambda x:np.log(x))      

