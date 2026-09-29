# Grafica una captura del prototipo y estima la distancia al objeto
# Uso: python3 ver_captura.py captura.txt [fondo.txt] [chirp.txt]

import sys
import numpy as np
import matplotlib.pyplot as plt

C = 343.0
FS_POR_DEFECTO = 20000.0
GUARDA_CM = 8.0
# Separacion entre el centro del parlante y el del microfono. El sistema mide
# (2d - s)/2, asi que hay que sumarle s/2 para obtener la distancia al objeto.
SEPARACION_CM = 7.0


def leer(ruta):
    fs, valores = FS_POR_DEFECTO, []
    for linea in open(ruta, encoding="utf-8", errors="ignore"):
        linea = linea.strip()
        if linea.startswith("#"):
            for campo in linea.split():
                if campo.startswith("fs="):
                    fs = float(campo[3:])
        else:
            for campo in linea.replace(",", " ").split():
                if campo.isdigit():
                    valores.append(int(campo))
    return fs, np.array(valores, dtype=float)


def chirp_de_referencia(fs, n=80, f0=2000.0, f1=5000.0):
    # Misma formula que el firmware, por si no se volco la tabla real
    t = np.arange(n)/fs
    k = (f1-f0)/(n/fs)
    m = n//10
    w = np.ones(n)
    w[:m] = 0.5-0.5*np.cos(np.pi*np.arange(m)/m)
    w[-m:] = w[:m][::-1]
    return w*np.sin(2*np.pi*(f0*t + 0.5*k*t**2))


def envolvente(rx, tx):
    # Correlacion analitica: se anulan las frecuencias negativas antes de la IFFT
    N = 1 << int(np.ceil(np.log2(len(rx)+len(tx))))
    R = np.fft.fft(rx, N)*np.conj(np.fft.fft(tx, N))
    R[N//2+1:] = 0
    R[1:N//2] *= 2
    return np.abs(np.fft.ifft(R))[:len(rx)-len(tx)+1]


def interpolar(env, k):
    if 0 < k < len(env)-1:
        y0, y1, y2 = env[k-1:k+2]
        den = y0-2*y1+y2
        if den != 0:
            return k + 0.5*(y0-y2)/den
    return float(k)


archivos = [a for a in sys.argv[1:] if not a.startswith("-")]
ruta = archivos[0] if archivos else "captura.txt"
fs, cuentas = leer(ruta)
if cuentas.size == 0:
    sys.exit(f"No se encontraron muestras en {ruta}")

x = cuentas - cuentas.mean()
tx = chirp_de_referencia(fs)
fondo = None
if len(archivos) > 1:
    _, crudo = leer(archivos[1])
    if crudo.size == x.size:
        fondo = crudo - crudo.mean()
    else:
        print(f"AVISO: {archivos[1]} tiene {crudo.size} muestras y no {x.size}, se ignora")
if len(archivos) > 2:
    tx = leer(archivos[2])[1] - leer(archivos[2])[1].mean()

mv = cuentas*3300.0/4095.0
print(f"muestras {cuentas.size} | fs {fs:.0f} Hz | continua {mv.mean():.0f} mV | "
      f"pico a pico {mv.max()-mv.min():.0f} mV")
if cuentas.max() >= 4090 or cuentas.min() <= 5:
    print("ATENCION: la senal se recorta, baje la amplitud")

env_cruda = envolvente(x, tx)
k_directo = int(np.argmax(env_cruda))

# El camino directo es la referencia de tiempo: cancela los retardos del parlante,
# el amplificador y el ADC, que son fijos y desconocidos.
env = envolvente(x-fondo, tx) if fondo is not None else env_cruda
distancia = (np.arange(env.size)-k_directo)*C/(2*fs)*100 + SEPARACION_CM/2

guarda = k_directo + int(np.ceil(GUARDA_CM/100*2*fs/C))
busqueda = env.copy()
busqueda[:guarda] = 0
k = int(np.argmax(busqueda))
d = (interpolar(env, k)-k_directo)*C/(2*fs)*100 + SEPARACION_CM/2
piso = np.median(busqueda[guarda:])

print(f"camino directo en la muestra {k_directo}")
print(f"objeto detectado a {d:.1f} cm   (pico {env[k]:.0f}, piso {piso:.0f}, "
      f"relacion {20*np.log10(env[k]/max(piso, 1e-9)):.1f} dB)")
if fondo is None:
    print("sin resta de fondo: el pico puede ser un reflejo del entorno y no el objeto")

fig, ejes = plt.subplots(3, 1, figsize=(9, 8))
ejes[0].plot(np.arange(cuentas.size)/fs*1e3, mv, color="#1f4e79", linewidth=0.9)
ejes[0].set_xlabel("Tiempo [ms]"); ejes[0].set_ylabel("Entrada [mV]")
ejes[0].set_title("Captura", loc="left")

f = np.fft.rfftfreq(x.size, 1/fs)
esp = np.abs(np.fft.rfft(x*np.hanning(x.size)))
ejes[1].plot(f, 20*np.log10(esp/max(esp.max(), 1e-12)), color="#c55a11", linewidth=0.9)
ejes[1].set_ylim(-80, 5); ejes[1].set_xlabel("Frecuencia [Hz]"); ejes[1].set_ylabel("[dB]")
ejes[1].set_title("Espectro", loc="left")

if fondo is not None:
    ejes[2].plot(distancia, env_cruda, color="#bbbbbb", linewidth=0.9, label="sin restar")
ejes[2].plot(distancia, env, color="#2e6f40", linewidth=1.1, label="correlacion")
ejes[2].axvline(d, color="#c00000", linestyle=":", linewidth=1.2, label=f"objeto: {d:.0f} cm")
ejes[2].axvline(0, color="k", linestyle=":", linewidth=0.8)
ejes[2].set_xlim(-15, distancia.max())
ejes[2].set_xlabel("Distancia respecto al camino directo [cm]")
ejes[2].set_ylabel("Correlacion")
ejes[2].set_title("Envolvente de la correlacion", loc="left")
ejes[2].legend(fontsize=8)

for eje in ejes:
    eje.grid(alpha=0.25, linewidth=0.5)
fig.tight_layout()
salida = ruta.rsplit(".", 1)[0] + ".png"
fig.savefig(salida, dpi=140)
print(f"grafica guardada en {salida}")
