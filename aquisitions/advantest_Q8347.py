# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 21:56:05 2019

Script para obtenção dos espectros com o Advantest Q8347 via GPIB.
O programa salva arquivos na pasta onde está script, um .txt e um .png.
O endereço do OSA pode ser alterado na linha 32 (rm.open_resource)

@author: Thiago
"""
from pathlib import Path
import pyvisa
from numpy import savetxt, asfarray, column_stack
import matplotlib
import matplotlib.pyplot as plt

def checkSTB(instrument, t):   #FUNÇÃO BASEADA NO "wait of spectrometer" do Gabriel
    from time import sleep  #função para verificar o status byte do aparelho
    sleep(1)                #Como o Advantest nao tem a função do GPIB (SRQ),
    i = True                #é preciso usar uma funçõa de mais baixo nivel
    while i:
        stb = instrument.read_stb()
        if stb == 1:        #Quando o STB é igual a 1, o Advantest terminou
            i = False       #de fazer a ação que havia sido solicitada
        else:
            sleep(t)
    return

rm = pyvisa.ResourceManager()
rm.list_resources()

osa = rm.open_resource('GPIB0::8::INSTR') #osa livre
osa.chunk_size = 65535 #configuraçoes da comunicaçao
osa.timeout = 20000 #configuraçoes da comunicaçao
osa.read_termination = '\n' #configuraçoes da comunicaçao


j=1
cont = True
while cont:
    print('Fazendo aquisição...')
    osa.write('MEA1') #faz medição SINGLE
    checkSTB(osa, 0.5) #espera final da medição
    x = osa.query('OSD1') #coleta dados do eixo
    y = osa.query('OSD0') #coleta dados do eixo y

    x = x[5:] #corta cabeçalho
    y = y[5:] #corta cabeçalho
    x_array = asfarray(x.split(',')) #converte para array de floats do numpy
    y_array = asfarray(y.split(',')) #converte para array de floats do numpy
    data = column_stack((x_array, y_array)) #cria array com dados
    folder = 'data/19052023/'
    test = 'y_plus_ref_'
    #handling files
    path_file = (Path().absolute()/(folder+test+str(j)+'.txt'))
    if (path_file.exists()):
        j+=1
    filename = folder+test+ str(j) + '.txt'
    savetxt(filename, data) #salva arquivo txt
    print('Arquivo salvo: '+ filename)

    fig, ax = plt.subplots()
    ax.plot(x_array, y_array)
    plt.ticklabel_format(axis='x', style='sci',scilimits=(-9,-9),useMathText=True)
    ax.grid()
    ax.minorticks_on()
    ax.set(xlabel='Wavelength')
    fig.savefig(filename[:len(filename)-4]+'.png', dpi=100)
    plt.show()
    print('Imagem salva: '+filename[:len(filename)-4]+'.png' )

    new = input('\nPress ENTER to new aquisition or N to exit\n')
    # new = "n"
    # if new.lower() == "n":
    cont = False