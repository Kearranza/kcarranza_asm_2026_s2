// Radar acustico en ESP32: chirp, FFT radix-2 y espectro de referencia

const int PIN_DAC = 25;
const int PIN_MIC = 34;

const float FS = 20000.0;
const int PERIODO_US = 50;
const int N_TX = 80;
const int N_FFT = 512;
const float F0 = 2000.0;
const float F1 = 5000.0;

int amplitud = 100;
uint8_t chirpDac[N_TX];
float chirpF[N_TX];
float re[N_FFT], im[N_FFT];
float txRe[N_FFT], txIm[N_FFT];


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

void emitirChirp() {
  int64_t t0 = esp_timer_get_time();
  for (int n = 0; n < N_TX; n++) {
    while (esp_timer_get_time() < t0 + (int64_t)n * PERIODO_US) {}
    dacWrite(PIN_DAC, chirpDac[n]);
  }
  dacWrite(PIN_DAC, 128);
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
  Serial.println("  e  emitir el chirp una vez");
  Serial.println("  t  tono de prueba de 1 s");
  Serial.printf("  + -  amplitud de salida (ahora %d de 127)\n", amplitud);
}

void setup() {
  Serial.begin(115200);
  delay(300);
  prepararChirp();
  dacWrite(PIN_DAC, 128);
  menu();
}

void loop() {
  if (!Serial.available()) return;
  char c = Serial.read();
  switch (c) {
    case 'e': emitirChirp(); break;
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
