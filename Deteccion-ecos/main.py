#Etapa 3: deteccion de ecos por correlacion. Corre los experimentos y deja todo en resultados/

import csv
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "Experimentos-fft"))
from graficas import plt, AZUL, NARANJA, VERDE, GRIS

import correlacion as co
import senales_eco as se

RAIZ = os.path.dirname(os.path.abspath(__file__))
SALIDA = os.path.join(RAIZ, "resultados")
FIGS = os.path.join(SALIDA, "figs")

DISTANCIAS = [0.30, 0.80, 1.50]
REFLECTIVIDADES = [0.6, 0.6, 0.6]
LARGO = 800
SNR_DB = 20.0
DISTANCIA_MINIMA = 0.25


def distancia_de_retardo(lag):
    return (lag*se.C_SONIDO/se.FS + se.SEPARACION)/2


def guarda():
    return int(np.ceil(se.retardo_de_distancia(DISTANCIA_MINIMA)))


def estimar(rx, tx, cuantos):
    env = co.envolvente(rx, tx)
    buscable = env.copy()
    buscable[:guarda()] = -np.inf
    separacion = int(0.10*2*se.FS/se.C_SONIDO)
    return env, co.picos(buscable, cuantos, separacion)


def escenario(semilla_objeto=1, semilla_fondo=2):
    #Dos capturas del mismo escenario: con los objetos y sin ellos
    tx, rx = se.escena(DISTANCIAS, REFLECTIVIDADES, LARGO, SNR_DB, semilla=semilla_objeto)
    _, fondo = se.escena([], [], LARGO, SNR_DB, semilla=semilla_fondo)
    return tx, rx, fondo


def validar(tx, rx, fondo):
    print("== Deteccion de ecos ==")
    filas = []
    for etiqueta, senal in [("sin restar el fondo", rx), ("restando el fondo", rx - fondo)]:
        env, estimados = estimar(senal, tx, len(DISTANCIAS))
        print(f"\n  {etiqueta}:")
        print("    real [m]  estimado [m]   error [cm]")
        for real, lag in zip(DISTANCIAS, estimados):
            d = distancia_de_retardo(lag)
            print(f"    {real:8.2f} {d:13.3f} {100*(d-real):12.1f}")
            filas.append({"metodo": etiqueta, "real_m": real,
                          "estimado_m": round(d, 4), "error_cm": round(100*(d-real), 2)})
    with open(os.path.join(SALIDA, "deteccion.csv"), "w", newline="", encoding="utf-8") as f:
        escritor = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
        escritor.writeheader()
        escritor.writerows(filas)


def comparar(tx, rx):
    directa = co.correlacion_directa(rx, tx)
    propia = co.correlacion_fft(rx, tx)
    referencia = co.correlacion_fft_referencia(rx, tx)
    e1 = np.max(np.abs(directa - propia))/np.max(np.abs(directa))
    e2 = np.max(np.abs(directa - referencia))/np.max(np.abs(directa))
    print(f"\n== Directa contra FFT ==")
    print(f"  diferencia relativa con nuestra FFT: {e1:.2e}")
    print(f"  diferencia relativa con numpy.fft:   {e2:.2e}")
    return directa, propia


def cronometrar(funcion, *args, repeticiones=3):
    mejor = float("inf")
    for _ in range(repeticiones):
        inicio = time.perf_counter()
        funcion(*args)
        mejor = min(mejor, time.perf_counter() - inicio)
    return mejor


def tiempos():
    tx = se.chirp()
    filas = []
    print("\n== Tiempos de ejecucion ==")
    print("      N   directa [ms]   FFT propia [ms]   FFT numpy [ms]")
    for p in range(8, 15):
        N = 2**p
        _, rx = se.escena(DISTANCIAS, REFLECTIVIDADES, N, SNR_DB)
        t_dir = cronometrar(co.correlacion_directa, rx, tx)
        t_propia = cronometrar(co.correlacion_fft, rx, tx)
        t_numpy = cronometrar(co.correlacion_fft_referencia, rx, tx)
        filas.append({"N": N, "directa_ms": round(t_dir*1e3, 3),
                      "fft_propia_ms": round(t_propia*1e3, 3),
                      "fft_numpy_ms": round(t_numpy*1e3, 4)})
        print(f"  {N:5d} {t_dir*1e3:13.3f} {t_propia*1e3:17.3f} {t_numpy*1e3:16.4f}")
    with open(os.path.join(SALIDA, "tiempos.csv"), "w", newline="", encoding="utf-8") as f:
        escritor = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
        escritor.writeheader()
        escritor.writerows(filas)
    return filas


def barrido_de_ruido():
    #Hasta que SNR sigue siendo fiable la deteccion del eco mas lejano
    tx = se.chirp()
    filas = []
    print("\n== Efecto del ruido sobre el eco mas lejano ==")
    for snr in range(30, -11, -5):
        errores = []
        for semilla in range(20):
            _, rx = se.escena(DISTANCIAS, REFLECTIVIDADES, LARGO, snr, semilla=semilla)
            _, fondo = se.escena([], [], LARGO, snr, semilla=semilla + 100)
            _, estimados = estimar(rx - fondo, tx, len(DISTANCIAS))
            errores.append(abs(distancia_de_retardo(estimados[-1]) - DISTANCIAS[-1]))
        aciertos = 100*np.mean(np.array(errores) < 0.05)
        filas.append({"snr_db": snr, "aciertos_pct": aciertos})
        print(f"  SNR {snr:4d} dB -> {aciertos:5.1f} % de aciertos")
    return filas


def figuras(tx, rx, fondo, directa, propia, filas_tiempo, filas_ruido):
    fig, ejes = plt.subplots(2, 1, figsize=(9, 5))
    ejes[0].plot(np.arange(tx.size)/se.FS*1e3, tx, color=AZUL)
    ejes[0].set_ylabel("Amplitud"); ejes[0].set_xlabel("Tiempo [ms]")
    ejes[0].set_title("Chirp transmitido (2 a 5 kHz, 4 ms)", loc="left", fontsize=9)
    ejes[1].plot(np.arange(rx.size)/se.FS*1e3, rx, color=NARANJA, linewidth=0.8)
    ejes[1].set_ylabel("Amplitud"); ejes[1].set_xlabel("Tiempo [ms]")
    ejes[1].set_title(f"Recibido: camino directo, tres ecos y ruido (SNR {SNR_DB:.0f} dB)",
                      loc="left", fontsize=9)
    for eje in ejes: eje.grid(alpha=0.25, linewidth=0.5)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "senales.png")); plt.close(fig)

    env_sin, _ = estimar(rx, tx, len(DISTANCIAS))
    env_con, estimados = estimar(rx - fondo, tx, len(DISTANCIAS))
    d = distancia_de_retardo(np.arange(env_con.size))*100
    fig, ejes = plt.subplots(2, 1, figsize=(9, 5.5))
    ejes[0].plot(d, directa, color=AZUL, linewidth=1.4, label="directa")
    ejes[0].plot(d, propia, color=NARANJA, linewidth=0.9, linestyle="--", label="por FFT")
    ejes[0].set_ylabel("Correlacion"); ejes[0].legend(fontsize=8)
    ejes[0].set_title("Las dos implementaciones dan el mismo resultado", loc="left", fontsize=9)
    ejes[1].plot(d, env_sin, color=GRIS, linewidth=0.9, label="sin restar el fondo")
    ejes[1].plot(d, env_con, color=VERDE, linewidth=1.2, label="restando el fondo")
    for real in DISTANCIAS:
        ejes[1].axvline(100*real, color="#c00000", linestyle=":", linewidth=0.9)
    ejes[1].set_xlabel("Distancia [cm]"); ejes[1].set_ylabel("Envolvente")
    ejes[1].set_title("Lineas rojas: posicion real de los tres objetos", loc="left", fontsize=9)
    ejes[1].legend(fontsize=8)
    for eje in ejes:
        eje.set_xlim(0, 200)
        eje.grid(alpha=0.25, linewidth=0.5)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "correlacion.png")); plt.close(fig)

    n = [f["N"] for f in filas_tiempo]
    fig, eje = plt.subplots(figsize=(7, 4))
    eje.loglog(n, [f["directa_ms"] for f in filas_tiempo], "o-", color=NARANJA, ms=4, label="directa")
    eje.loglog(n, [f["fft_propia_ms"] for f in filas_tiempo], "s-", color=AZUL, ms=4, label="FFT propia")
    eje.loglog(n, [f["fft_numpy_ms"] for f in filas_tiempo], "d-", color=VERDE, ms=4, label="FFT de numpy")
    eje.set_xlabel("Muestras de la senal recibida $N$"); eje.set_ylabel("Tiempo [ms]")
    eje.grid(True, which="both", alpha=0.25, linewidth=0.5); eje.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "tiempos.png")); plt.close(fig)

    fig, eje = plt.subplots(figsize=(7, 3.4))
    eje.plot([f["snr_db"] for f in filas_ruido], [f["aciertos_pct"] for f in filas_ruido],
             "o-", color=VERDE, ms=4)
    eje.set_xlabel("SNR [dB]"); eje.set_ylabel("Aciertos [%]")
    eje.set_title("Deteccion del eco a 1,5 m contra el nivel de ruido", loc="left", fontsize=9)
    eje.grid(alpha=0.25, linewidth=0.5)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "ruido.png")); plt.close(fig)


def main():
    os.makedirs(FIGS, exist_ok=True)
    tx, rx, fondo = escenario()
    validar(tx, rx, fondo)
    directa, propia = comparar(tx, rx)
    filas_tiempo = tiempos()
    filas_ruido = barrido_de_ruido()
    figuras(tx, rx, fondo, directa, propia, filas_tiempo, filas_ruido)
    print(f"\nResultados en {SALIDA}")


if __name__ == "__main__":
    main()
