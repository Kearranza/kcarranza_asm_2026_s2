#Simulacion de la senal transmitida y sus ecos con retardos conocidos

import numpy as np

FS = 20000.0
C_SONIDO = 343.0
F0 = 2000.0
F1 = 5000.0
DURACION_CHIRP = 0.004
SEPARACION = 0.07


def chirp(fs=FS, duracion=DURACION_CHIRP, f0=F0, f1=F1):
    #Mismo chirp que emite el microcontrolador, con ventana Tukey del 20 %
    n = int(round(duracion*fs))
    t = np.arange(n)/fs
    k = (f1-f0)/duracion
    m = n//10
    w = np.ones(n)
    w[:m] = 0.5-0.5*np.cos(np.pi*np.arange(m)/m)
    w[-m:] = w[:m][::-1]
    return w*np.sin(2*np.pi*(f0*t + 0.5*k*t**2))


def retardar(x, retardo_muestras, largo, fs=FS):
    #Retardo fraccionario exacto: en frecuencia es multiplicar por exp(-j*2*pi*f*tau)
    N = 1 << int(np.ceil(np.log2(largo + x.size + 64)))
    X = np.fft.rfft(np.concatenate([x, np.zeros(N-x.size)]))
    f = np.fft.rfftfreq(N)
    return np.fft.irfft(X*np.exp(-2j*np.pi*f*retardo_muestras), N)[:largo]


def retardo_de_distancia(d, fs=FS):
    #El eco recorre 2d; el camino directo recorre la separacion entre transductores
    return (2*d - SEPARACION)/C_SONIDO*fs


def escena(distancias, reflectividades, largo, snr_db=20.0, fs=FS, semilla=1):
    #Devuelve la senal recibida: camino directo, un eco por objeto y ruido blanco
    tx = chirp(fs)
    rx = 3.0*retardar(tx, 0.0, largo, fs)
    for d, rho in zip(distancias, reflectividades):
        rx += rho/(2*d)*retardar(tx, retardo_de_distancia(d, fs), largo, fs)

    rng = np.random.default_rng(semilla)
    potencia_ruido = np.mean(rx**2)/(10**(snr_db/10))
    return tx, rx + rng.normal(0.0, np.sqrt(potencia_ruido), largo)
