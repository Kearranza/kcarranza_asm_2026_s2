# Pide una captura al ESP32 por el puerto serial y la guarda en un archivo
# Uso: python3 capturar_serial.py COM5 [captura.txt] [--fondo]
# Necesita pyserial:  pip install pyserial

import subprocess
import sys
import time

import serial

puerto = sys.argv[1] if len(sys.argv) > 1 else "COM5"
salida = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].startswith("-") else "captura.txt"
comando = b"r" if "--fondo" in sys.argv else b"p"

with serial.Serial(puerto, 115200, timeout=2) as puerto_serie:
    time.sleep(2)
    puerto_serie.reset_input_buffer()

    puerto_serie.write(comando)
    time.sleep(1.5)
    print(puerto_serie.read_all().decode(errors="ignore").strip())

    puerto_serie.write(b"d")
    lineas = []
    while True:
        linea = puerto_serie.readline().decode(errors="ignore").strip()
        if not linea:
            break
        lineas.append(linea)
        if linea.startswith("# fin"):
            break

if len(lineas) < 2:
    sys.exit("No llego el volcado. Revisá el puerto y que no esté abierto el monitor serial del IDE.")

with open(salida, "w", encoding="utf-8") as f:
    f.write("\n".join(lineas) + "\n")
print(f"guardado en {salida}")

subprocess.run([sys.executable, "ver_captura.py", salida])
