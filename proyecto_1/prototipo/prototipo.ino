// Prototipo de puesta en marcha del radar acustico - ESP32 + KY-037
// Genera un chirp por el DAC y captura por el ADC. El analisis se hace en la PC.

const int PIN_DAC = 25;
const int PIN_MIC = 34;

const float FS = 20000.0;
const int PERIODO_US = 50;
const int N_TX = 80;
const int N_REC = 400;
const float F0 = 2000.0;
const float F1 = 5000.0;

uint8_t chirp[N_TX];
uint16_t captura[N_REC];
bool hubo_emision = false;
int32_t jitter_max = 0;
int atrasadas = 0;
int peor_indice = -1;
int amplitud = 100;

void generarChirp() {
  float T = N_TX / FS;
  float k = (F1 - F0) / T;
  int m = N_TX / 10;
  for (int n = 0; n < N_TX; n++) {
    float t = n / FS;
    // Ventana Tukey del 20 %: evita el clic del parlante al arrancar y al parar
    float w = 1.0;
    if (n < m) w = 0.5 - 0.5 * cos(PI * n / m);
    else if (n >= N_TX - m) w = 0.5 - 0.5 * cos(PI * (N_TX - 1 - n) / m);
    float s = sin(2 * PI * (F0 * t + 0.5 * k * t * t));
    chirp[n] = (uint8_t)(127.5 + amplitud * w * s);
  }
}

void capturar(bool emitir) {
  jitter_max = 0;
  atrasadas = 0;
  peor_indice = -1;
  // La primera lectura del ADC calibra el conversor y tarda mucho mas que las siguientes
  for (int i = 0; i < 8; i++) { analogRead(PIN_MIC); dacWrite(PIN_DAC, 128); }
  int64_t t0 = esp_timer_get_time();
  for (int n = 0; n < N_REC; n++) {
    int64_t objetivo = t0 + (int64_t)n * PERIODO_US;
    while (esp_timer_get_time() < objetivo) {}
    int32_t desvio = (int32_t)(esp_timer_get_time() - objetivo);
    if (desvio > jitter_max) { jitter_max = desvio; peor_indice = n; }
    if (desvio > 10) atrasadas++;
    if (emitir) dacWrite(PIN_DAC, n < N_TX ? chirp[n] : 128);
    captura[n] = analogRead(PIN_MIC);
  }
  dacWrite(PIN_DAC, 128);
  hubo_emision = emitir;
}

void capturarTono(float f) {
  for (int i = 0; i < 8; i++) { analogRead(PIN_MIC); dacWrite(PIN_DAC, 128); }
  int64_t t0 = esp_timer_get_time();
  for (int n = 0; n < N_REC; n++) {
    while (esp_timer_get_time() < t0 + (int64_t)n * PERIODO_US) {}
    dacWrite(PIN_DAC, (uint8_t)(127.5 + amplitud * sin(2 * PI * f * n / FS)));
    captura[n] = analogRead(PIN_MIC);
  }
  dacWrite(PIN_DAC, 128);
}

float amplitudEnFrecuencia(float f) {
  // DFT de un solo bin: rechaza el ruido fuera de esa frecuencia
  double media = 0;
  for (int n = 0; n < N_REC; n++) media += captura[n];
  media /= N_REC;
  double I = 0, Q = 0;
  for (int n = 0; n < N_REC; n++) {
    double ang = 2 * PI * f * n / FS;
    I += (captura[n] - media) * cos(ang);
    Q += (captura[n] - media) * sin(ang);
  }
  return 2.0 * sqrt(I * I + Q * Q) / N_REC;
}

void barrido() {
  const float frecuencias[] = {1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000};
  Serial.printf("Respuesta de la cadena completa (amplitud de salida %d/127):\n", amplitud);
  Serial.println("  Hz     amplitud [cuentas]");
  for (unsigned i = 0; i < sizeof(frecuencias) / sizeof(float); i++) {
    capturarTono(frecuencias[i]);
    uint16_t minimo = 4095, maximo = 0;
    for (int n = 0; n < N_REC; n++) {
      if (captura[n] < minimo) minimo = captura[n];
      if (captura[n] > maximo) maximo = captura[n];
    }
    bool recorta = (maximo >= 4090 || minimo <= 5);
    Serial.printf("  %5.0f   %7.1f %s\n", frecuencias[i],
                  amplitudEnFrecuencia(frecuencias[i]), recorta ? " RECORTA" : "");
    delay(80);
  }
}

void estadisticas() {
  long suma = 0;
  uint16_t minimo = 4095, maximo = 0;
  for (int n = 0; n < N_REC; n++) {
    suma += captura[n];
    if (captura[n] < minimo) minimo = captura[n];
    if (captura[n] > maximo) maximo = captura[n];
  }
  float media = (float)suma / N_REC;
  double acum = 0;
  for (int n = 0; n < N_REC; n++) acum += pow(captura[n] - media, 2);
  float rms = sqrt(acum / N_REC);

  Serial.printf("  continua %.1f cuentas (%.0f mV)\n", media, media * 3300.0 / 4095.0);
  Serial.printf("  minimo %u, maximo %u, pico a pico %u cuentas (%.0f mV)\n",
                minimo, maximo, maximo - minimo, (maximo - minimo) * 3300.0 / 4095.0);
  Serial.printf("  rms %.1f cuentas (%.1f mV)\n", rms, rms * 3300.0 / 4095.0);
  Serial.printf("  jitter maximo %ld us en la muestra %d (periodo %d us)\n",
                (long)jitter_max, peor_indice, PERIODO_US);
  Serial.printf("  muestras atrasadas mas de 10 us: %d de %d\n", atrasadas, N_REC);
  if (maximo >= 4090 || minimo <= 5) Serial.println("  ATENCION: la senal se esta recortando, baje la ganancia");
}

void volcarCaptura() {
  Serial.printf("# captura fs=%.0f n=%d emision=%d\n", FS, N_REC, hubo_emision ? 1 : 0);
  for (int n = 0; n < N_REC; n++) {
    Serial.print(captura[n]);
    Serial.print((n % 20 == 19 || n == N_REC - 1) ? '\n' : ',');
  }
  Serial.println("# fin");
}

void volcarChirp() {
  Serial.printf("# chirp fs=%.0f n=%d f0=%.0f f1=%.0f\n", FS, N_TX, F0, F1);
  for (int n = 0; n < N_TX; n++) {
    Serial.print(chirp[n]);
    Serial.print((n % 20 == 19 || n == N_TX - 1) ? '\n' : ',');
  }
  Serial.println("# fin");
}

void tonoDePrueba() {
  Serial.println("Tono de 3 kHz durante 1 s...");
  int64_t t0 = esp_timer_get_time();
  for (long n = 0; n < (long)FS; n++) {
    while (esp_timer_get_time() < t0 + (int64_t)n * PERIODO_US) {}
    dacWrite(PIN_DAC, (uint8_t)(127.5 + amplitud * sin(2 * PI * 3000.0 * n / FS)));
  }
  dacWrite(PIN_DAC, 128);
  Serial.println("Listo. Si no se escucho nada, revise amplificador y parlante.");
}

void pruebaDeContinua() {
  const uint8_t valores[] = {0, 64, 128, 192, 255};
  Serial.println("Medi D25 contra GND con el multimetro, antes del capacitor:");
  for (unsigned i = 0; i < sizeof(valores); i++) {
    dacWrite(PIN_DAC, valores[i]);
    Serial.printf("  DAC = %3u  ->  esperado %.2f V   (5 s para medir)\n",
                  valores[i], valores[i] * 3.3 / 255.0);
    delay(5000);
  }
  dacWrite(PIN_DAC, 128);
  Serial.println("Listo, el DAC vuelve a 128 (1,65 V).");
}

void tonoContinuo() {
  Serial.println("Tono continuo de 3 kHz. Mandá cualquier tecla para parar.");
  while (Serial.available()) Serial.read();
  long n = 0;
  int64_t t0 = esp_timer_get_time();
  while (!Serial.available()) {
    while (esp_timer_get_time() < t0 + (int64_t)n * PERIODO_US) {}
    dacWrite(PIN_DAC, (uint8_t)(127.5 + amplitud * sin(2 * PI * 3000.0 * n / FS)));
    n++;
    if (n > 1000000) { n = 0; t0 = esp_timer_get_time(); }
  }
  while (Serial.available()) Serial.read();
  dacWrite(PIN_DAC, 128);
  Serial.println("Tono detenido.");
}

void menu() {
  Serial.println("\n--- Prototipo radar acustico ---");
  Serial.println("  t  tono de prueba de 1 s");
  Serial.println("  r  captura sin emitir (ruido de fondo)");
  Serial.println("  p  ping: emite el chirp y captura");
  Serial.println("  d  volcar la ultima captura");
  Serial.println("  x  volcar la tabla del chirp");
  Serial.println("  b  barrido: respuesta de la cadena entre 1 y 7 kHz");
  Serial.println("  v  prueba de continua del DAC (para medir con multimetro)");
  Serial.println("  c  tono continuo hasta que mandes otra tecla");
  Serial.printf("  + -  subir o bajar la amplitud de salida (ahora %d de 127)\n", amplitud);
}

void setup() {
  Serial.begin(115200);
  delay(300);
  analogReadResolution(12);
  analogSetPinAttenuation(PIN_MIC, ADC_11db);
  generarChirp();
  dacWrite(PIN_DAC, 128);
  menu();
}

void loop() {
  if (!Serial.available()) return;
  char c = Serial.read();
  switch (c) {
    case 't': tonoDePrueba(); break;
    case 'r': capturar(false); Serial.println("Ruido de fondo:"); estadisticas(); break;
    case 'p': capturar(true);  Serial.println("Ping:");           estadisticas(); break;
    case 'd': volcarCaptura(); break;
    case 'x': volcarChirp();   break;
    case 'b': barrido(); break;
    case 'v': pruebaDeContinua(); break;
    case 'c': tonoContinuo(); break;
    case '+': case '-':
      amplitud = constrain(amplitud + (c == '+' ? 15 : -15), 10, 127);
      generarChirp();
      Serial.printf("Amplitud de salida: %d de 127\n", amplitud);
      break;
    case '\n': case '\r': break;
    default: menu();
  }
}
