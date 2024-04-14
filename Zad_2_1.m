% Parametry
% PIERWSZA KROPKA
sigma = sqrt(0.64); % Odchylenie standardowe, przykładowa wartość dla wariancji σ^2 = 0.64
Tp = 0.001; % Okres próbkowania
N = 2000; % Liczba próbek
n = 0:N-1; % Wektor indeksów próbek

% Definicja sygnałów
x = sin(2*pi*5*n*Tp) + 0.5*sin(2*pi*10*n*Tp) + 0.25*sin(2*pi*30*n*Tp); % Sygnał x(nTp)
e = sigma * randn(1, N); % Sygnał e(nTp) - szum biały

% Definicja i zastosowanie filtru
H = tf(0.1, [1 -0.9], Tp); % Dyskretny filtr dolnoprzepustowy H(q^-1)
tn = n*Tp; % Wektor czasu
v = lsim(H, e, tn); % Sygnał v(nTp) po filtracji
% Pierwsza kropka
% Wyświetlenie wykresów
figure;
subplot(3,1,1);
plot(tn, x);
title('Sygnał x(nTp)');
xlabel('Czas [s]');
ylabel('Amplituda');

subplot(3,1,2);
plot(tn, e);
title('Szum biały e(nTp)');
xlabel('Czas [s]');
ylabel('Amplituda');

subplot(3,1,3);
plot(tn, v);
title('Sygnał v(nTp) po filtracji');
xlabel('Czas [s]');
ylabel('Amplituda');
% Druga kropka
% Obliczenie FFT i przeskalowanie
X = fft(x) *2/N;

% Obliczenie modułu i normalizacja widma do wysokości prążków
X_magnitude = abs(X);

% Skala częstotliwości
f = (0:N-1)*(1/(Tp*N));

% Wykreślenie widma amplitudowego
figure;
stem(f, X_magnitude);
title('Widmo amplitudowe |X_N(j\omega_k)| sygnału x(nT_p)');
xlabel('Częstotliwość [Hz]');
ylabel('|X_N(j\omega_k)|');
xlim([0 50]); % Ograniczenie zakresu dla lepszej czytelności
% TRZECIA KROPKA
Tp = 0.001; % Okres próbkowania

% Definicja sygnału x(nTp) dla różnych wartości N
N_values = [1000, 200, 100];

for N = N_values
    n = 0:N-1; % Wektor próbek
    x = sin(2*pi*5*n*Tp) + 0.5*sin(2*pi*10*n*Tp) + 0.25*sin(2*pi*30*n*Tp);
    X = fft(x) *2 /N;
    X_magnitude = abs(X);
    f = (0:N-1)*(1/(Tp*N));

    % Wykreślenie widma amplitudowego
    figure;
    stem(f, X_magnitude);
    title(['Widmo amplitudowe |X_N(j\omega_k)|, N = ' num2str(N)]);
    xlabel('Częstotliwość [Hz]');
    ylabel('|X_N(j\omega_k)|');
    xlim([0 50]);
end
% CZWARA KROPKA
% Parametry
Tp = 0.001; % Okres próbkowania
N = 2000; % Liczba próbek
n = 0:N-1; % Wektor próbek

% Generacja sygnału
x = sin(2*pi*5*n*Tp) + 0.5*sin(2*pi*10*n*Tp) + 0.25*sin(2*pi*30*n*Tp);

% Energia sygnału w czasie
E_time = Tp * sum(x.^2);

% Obliczenie DFT sygnału x
X = fft(x);

% Energia sygnału w dziedzinie częstotliwości (z uwzględnieniem skalowania)
E_freq = (1/(N*Tp)) * sum(abs(X).^2);
% PERSEVALA TWIERDZENIE
% Porównanie energii sygnału w czasie i częstotliwości
fprintf('Energia sygnału w czasie: %f\n', E_time);
fprintf('Energia sygnału w częstotliwości: %f\n', E_freq/1000000);
fprintf('Różnica: %e\n', abs(E_time - E_freq/1000000));
% PIĄTA KROPKA
sigma_2=[0.64, 1.0];
sigma=sqrt(sigma_2);
e1=sigma(1)*randn(1,N);
e2=sigma(2)*randn(1,N);
N=2000;
n=0:N-1;
est_14_e1=(Tp/N)*power(abs(sum(e1.*exp(-complex(0,1)*n.*omega_k))),2);
est_14_e2=(Tp/N)*power(abs(sum(e2.*exp(-complex(0,1)*n.*omega_k))),2);

r_xx = xcorr(e1, 'biased');
okna = [N, round(N/2), round(N/4)];
figure;

for i = 1:length(okna)
    Mw = okna(i);
    w = hamming(Mw)';
    padding = (length(r_xx) - Mw) / 2;
    w_padded = [zeros(1, floor(padding)), w, zeros(1, ceil(padding))];
    r_xx_win = r_xx .* w_padded;
    phi_xx = fft(r_xx_win);
    f = (0:length(phi_xx)-1) / (length(phi_xx)*Tp);
    plot(f, abs(phi_xx));
    hold on;
end
%xlim([0 50])
legend('N', 'N/2', 'N/4');
title('Estymata gęstości widmowej mocy e1(nTp)');

N=2000;
n=0:N-1;
est_14_e1=(Tp/N)*power(abs(sum(e1.*exp(-complex(0,1)*n.*omega_k))),2);
est_14_e2=(Tp/N)*power(abs(sum(e2.*exp(-complex(0,1)*n.*omega_k))),2);

r_xx = xcorr(e2, 'biased');
okna = [N, round(N/2), round(N/4)];
figure;

for i = 1:length(okna)
    Mw = okna(i);
    w = hamming(Mw)';
    padding = (length(r_xx) - Mw) / 2;
    w_padded = [zeros(1, floor(padding)), w, zeros(1, ceil(padding))];
    r_xx_win = r_xx .* w_padded;
    phi_xx = fft(r_xx_win);
    f = (0:length(phi_xx)-1) / (length(phi_xx)*Tp);
    plot(f, abs(phi_xx));
    hold on;
end

legend('N', 'N/2', 'N/4');
title('Estymata gęstości widmowej mocy e2(nTp)');

% Siódma kropka
N=2000;
n=0:N-1;
est_14_v1=(Tp/N)*power(abs(sum(v1.*exp(-complex(0,1)*n.*omega_k))),2);
est_14_v2=(Tp/N)*power(abs(sum(v2.*exp(-complex(0,1)*n.*omega_k))),2);

r_xx2 = xcorr(v1, 'biased')';
r_xx2_u=xcorr(v1, 'unbiased')';
okna = [N];
figure;

for i = 1:length(okna)
    Mw2 = okna(i);
    w2 = hamming(Mw2)';
    padding2 = (length(r_xx) - Mw2) / 2;
    w_padded2 = [zeros(1, floor(padding2)), w2, zeros(1, ceil(padding2))];
    r_xx_win2 = r_xx2 .* w_padded2;
    phi_xx2 = fft(r_xx_win2);
    f2 = (0:length(phi_xx2)-1) / (length(phi_xx2)*Tp);
    plot(f2, abs(phi_xx2));
    hold on;
end
legend('N');
title('Estymata gęstości widmowej mocy v1(nTp)');

figure;

for i = 1:length(okna)
    Mw2 = okna(i);
    w2 = hamming(Mw2)';
    padding2 = (length(r_xx) - Mw2) / 2;
    w_padded2 = [zeros(1, floor(padding2)), w2, zeros(1, ceil(padding2))];
    r_xx_win2 = r_xx2_u .* w_padded2;
    phi_xx2 = fft(r_xx_win2);
    f2 = (0:length(phi_xx2)-1) / (length(phi_xx2)*Tp);
    plot(f2, abs(phi_xx2));
    hold on;
end

legend('N');
title('Estymata gęstości widmowej mocy v2(nTp)');

