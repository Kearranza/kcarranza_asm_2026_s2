// Radar acustico en ESP32: generacion del chirp y salida por el DAC

const int PIN_DAC = 25;
const int PIN_MIC = 34;

const float FS = 20000.0;
const int PERIODO_US = 50;
const int N_TX = 80;
const float F0 = 2000.0;
const float F1 = 5000.0;

int amplitud = 100;
uint8_t chirpDac[N_TX];
float chirpF[N_TX];


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
