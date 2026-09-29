// Radar acustico en ESP32: genera el chirp, captura el eco, correlaciona por FFT y estima la distancia

const int PIN_DAC = 25;
const int PIN_MIC = 34;

const float FS = 20000.0;
const int PERIODO_US = 50;
const int N_TX = 80;
const int N_REC = 400;
const int N_FFT = 512;
const int N_ENV = N_REC - N_TX + 1;
const float F0 = 2000.0;
const float F1 = 5000.0;

const float C_SONIDO = 343.0;
const float SEPARACION_CM = 7.0;
const float GUARDA_CM = 8.0;
const float UMBRAL_DB = 10.0;
const int PROMEDIO = 8;

int amplitud = 100;
uint8_t chirpDac[N_TX];
float chirpF[N_TX];
float re[N_FFT], im[N_FFT];
float txRe[N_FFT], txIm[N_FFT];
float fondo[N_REC], rx[N_REC], env[N_ENV];
bool hayFondo = false;
int kDirecto = 0;


void fft(float *pre, float *pim, int n, bool inversa) {
  for (int i = 1, j = 0; i < n; i++) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      float t = pre[i]; pre[i] = pre[j]; pre[j] = t;
      t = pim[i]; pim[i] = pim[j]; pim[j] = t;
    }
  }
  for (int tam = 2; tam <= n; tam <<= 1) {
    int mitad = tam >> 1;
    float ang = (inversa ? 2.0f : -2.0f) * PI / tam;
    for (int i = 0; i < mitad; i++) {
      float wr = cosf(ang * i), wi = sinf(ang * i);
      for (int inicio = 0; inicio < n; inicio += tam) {
        int a = inicio + i, b = a + mitad;
        float br = pre[b] * wr - pim[b] * wi;
        float bi = pre[b] * wi + pim[b] * wr;
        pre[b] = pre[a] - br; pim[b] = pim[a] - bi;
        pre[a] += br;         pim[a] += bi;
      }
    }
  }
  if (inversa) {
    for (int i = 0; i < n; i++) { pre[i] /= n; pim[i] /= n; }
  }
}

void prepararChirp() {
  float T = N_TX / FS;
  float k = (F1 - F0) / T;
  int m = N_TX / 10;
  for (int n = 0; n < N_TX; n++) {
    float t = n / FS;
    // Ventana Tukey del 20 %: evita el clic del parlante al arrancar y al parar
    float w = 1.0;
    if (n < m) w = 0.5 - 0.5 * cos(PI * n / m);
    else if (n >= N_TX - m) w = 0.5 - 0.5 * cos(PI * (N_TX - 1 - n) / m);
    chirpF[n] = w * sin(2 * PI * (F0 * t + 0.5 * k * t * t));
    chirpDac[n] = (uint8_t)(127.5 + amplitud * chirpF[n]);
  }
  // Se guarda el espectro conjugado del chirp: la correlacion es un producto en frecuencia
  for (int i = 0; i < N_FFT; i++) { re[i] = (i < N_TX) ? chirpF[i] : 0.0f; im[i] = 0.0f; }
  fft(re, im, N_FFT, false);
  for (int i = 0; i < N_FFT; i++) { txRe[i] = re[i]; txIm[i] = -im[i]; }
}

void capturar(float *destino, int repeticiones) {
  static int32_t acumulador[N_REC];
  for (int n = 0; n < N_REC; n++) acumulador[n] = 0;

  for (int r = 0; r < repeticiones; r++) {
    // La primera lectura del ADC calibra el conversor y tarda mucho mas que las siguientes
    for (int i = 0; i < 8; i++) { analogRead(PIN_MIC); dacWrite(PIN_DAC, 128); }
    int64_t t0 = esp_timer_get_time();
    for (int n = 0; n < N_REC; n++) {
      while (esp_timer_get_time() < t0 + (int64_t)n * PERIODO_US) {}
      dacWrite(PIN_DAC, n < N_TX ? chirpDac[n] : 128);
      acumulador[n] += analogRead(PIN_MIC);
    }
    dacWrite(PIN_DAC, 128);
    delay(30);
  }

  float media = 0;
  for (int n = 0; n < N_REC; n++) media += acumulador[n];
  media /= (float)N_REC * repeticiones;
  for (int n = 0; n < N_REC; n++) destino[n] = acumulador[n] / (float)repeticiones - media;
}

void correlacionar(const float *senal) {
  for (int i = 0; i < N_FFT; i++) { re[i] = (i < N_REC) ? senal[i] : 0.0f; im[i] = 0.0f; }
  fft(re, im, N_FFT, false);
  for (int i = 0; i < N_FFT; i++) {
    float pr = re[i] * txRe[i] - im[i] * txIm[i];
    float pi = re[i] * txIm[i] + im[i] * txRe[i];
    re[i] = pr; im[i] = pi;
  }
  // Anular las frecuencias negativas hace que la IFFT devuelva la senal analitica,
  // asi el pico se busca sobre la envolvente y no sobre una correlacion que oscila
  for (int i = N_FFT / 2 + 1; i < N_FFT; i++) { re[i] = 0.0f; im[i] = 0.0f; }
  for (int i = 1; i < N_FFT / 2; i++) { re[i] *= 2.0f; im[i] *= 2.0f; }
  fft(re, im, N_FFT, true);
  for (int n = 0; n < N_ENV; n++) env[n] = sqrtf(re[n] * re[n] + im[n] * im[n]);
}

int indiceDelMaximo(int desde, int hasta) {
  int k = desde;
  for (int n = desde; n < hasta; n++) if (env[n] > env[k]) k = n;
  return k;
}

float interpolar(int k) {
  if (k <= 0 || k >= N_ENV - 1) return k;
  float den = env[k - 1] - 2 * env[k] + env[k + 1];
  return (den == 0) ? k : k + 0.5f * (env[k - 1] - env[k + 1]) / den;
}

void calibrar() {
  Serial.println("Calibrando: quite el objeto y no se mueva...");
  capturar(fondo, PROMEDIO);
  correlacionar(fondo);
  // El camino directo es la referencia de tiempo: cancela los retardos fijos
  // del parlante, el amplificador y el ADC sin tener que medirlos
  kDirecto = indiceDelMaximo(0, N_ENV);
  hayFondo = true;
  Serial.printf("Listo. Camino directo en la muestra %d.\n", kDirecto);
}

void medir() {
  if (!hayFondo) { Serial.println("Falta calibrar: mande k sin el objeto."); return; }

  capturar(rx, PROMEDIO);
  for (int n = 0; n < N_REC; n++) rx[n] -= fondo[n];
  correlacionar(rx);

  int guarda = kDirecto + (int)ceilf(GUARDA_CM / 100.0f * 2 * FS / C_SONIDO);
  int k = indiceDelMaximo(guarda, N_ENV);

  float piso = 0;
  for (int n = guarda; n < N_ENV; n++) piso += env[n];
  piso /= (N_ENV - guarda);
  float calidad = 20 * log10f(env[k] / (piso > 0 ? piso : 1e-9f));

  float d = (interpolar(k) - kDirecto) * C_SONIDO / (2 * FS) * 100 + SEPARACION_CM / 2;

  if (calidad < UMBRAL_DB) Serial.printf("Sin deteccion  (mejor pico %.1f dB)\n", calidad);
  else Serial.printf("Distancia: %6.1f cm   calidad %4.1f dB\n", d, calidad);
}

void tonoDePrueba() {
  int64_t t0 = esp_timer_get_time();
  for (long n = 0; n < (long)FS; n++) {
    while (esp_timer_get_time() < t0 + (int64_t)n * PERIODO_US) {}
    dacWrite(PIN_DAC, (uint8_t)(127.5 + amplitud * sin(2 * PI * 3000.0 * n / FS)));
  }
  dacWrite(PIN_DAC, 128);
}

void menu() {
  Serial.println("\n--- Radar acustico ---");
  Serial.println("  k  calibrar el fondo (sin el objeto)");
  Serial.println("  m  medir una vez");
  Serial.println("  l  medir en bucle hasta que mandes otra tecla");
  Serial.println("  t  tono de prueba de 1 s");
  Serial.printf("  + -  amplitud de salida (ahora %d de 127)\n", amplitud);
}

void setup() {
  Serial.begin(115200);
  delay(300);
  analogReadResolution(12);
  analogSetPinAttenuation(PIN_MIC, ADC_11db);
  prepararChirp();
  dacWrite(PIN_DAC, 128);
  menu();
}

void loop() {
  if (!Serial.available()) return;
  char c = Serial.read();
  switch (c) {
    case 'k': calibrar(); break;
    case 'm': medir(); break;
    case 'l':
      while (Serial.available()) Serial.read();
      while (!Serial.available()) { medir(); delay(300); }
      while (Serial.available()) Serial.read();
      break;
    case 't': tonoDePrueba(); break;
    case '+': case '-':
      amplitud = constrain(amplitud + (c == '+' ? 15 : -15), 10, 127);
      prepararChirp();
      Serial.printf("Amplitud: %d de 127\n", amplitud);
      break;
    case '\n': case '\r': break;
    default: menu();
  }
}
