#Correlacion cruzada en el tiempo y en frecuencia, y estimacion del retardo

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "Experimentos-fft"))
from dft_fft import fft, ifft, siguiente_potencia_de_dos


def correlacion_directa(rx, tx):
    #Definicion: para cada retardo, la suma de productos con la replica emitida
    salida = np.zeros(rx.size - tx.size + 1)
    for n in range(salida.size):
        salida[n] = np.dot(rx[n:n + tx.size], tx)
    return salida


def correlacion_fft(rx, tx):
    #La correlacion es una convolucion con la replica invertida, y convolucionar
    #en el tiempo equivale a multiplicar en frecuencia: R = RX * conj(TX)
    N = siguiente_potencia_de_dos(rx.size + tx.size)
    RX = fft(np.concatenate([rx, np.zeros(N - rx.size)]))
    TX = fft(np.concatenate([tx, np.zeros(N - tx.size)]))
    return np.real(ifft(RX*np.conj(TX)))[:rx.size - tx.size + 1]


def correlacion_fft_referencia(rx, tx):
    #Mismo algoritmo pero con la FFT de numpy, para separar el costo del
    #algoritmo del costo de nuestra implementacion
    N = siguiente_potencia_de_dos(rx.size + tx.size)
    RX = np.fft.fft(rx, N)
    TX = np.fft.fft(tx, N)
    return np.real(np.fft.ifft(RX*np.conj(TX)))[:rx.size - tx.size + 1]


def envolvente(rx, tx):
    #Anular las frecuencias negativas devuelve la senal analitica: el pico deja
    #de oscilar a la frecuencia portadora y no salta de un ciclo al vecino
    N = siguiente_potencia_de_dos(rx.size + tx.size)
    R = fft(np.concatenate([rx, np.zeros(N - rx.size)])) * \
        np.conj(fft(np.concatenate([tx, np.zeros(N - tx.size)])))
    R[N//2 + 1:] = 0
    R[1:N//2] *= 2
    return np.abs(ifft(R))[:rx.size - tx.size + 1]


def interpolar(curva, k):
    #Interpolacion parabolica sobre los tres puntos alrededor del maximo
    if 0 < k < curva.size - 1:
        y0, y1, y2 = curva[k-1:k+2]
        den = y0 - 2*y1 + y2
        if den != 0:
            return k + 0.5*(y0 - y2)/den
    return float(k)


def picos(curva, cantidad, separacion):
    #Los maximos mas altos, exigiendo una separacion minima entre ellos
    encontrados = []
    quedan = curva.copy()
    for _ in range(cantidad):
        k = int(np.argmax(quedan))
        encontrados.append(interpolar(curva, k))
        quedan[max(0, k - separacion):k + separacion] = -np.inf
    return sorted(encontrados)
