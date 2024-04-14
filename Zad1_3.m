% Parametry
% PIERWSZA KROPKA
sigma = sqrt(0.64); % Odchylenie standardowe
Tp = 0.001; % Okres próbkowania
N = 2000; % Liczba próbek
n = 0:N-1; % Wektor chwil czasu dyskretnego
tn = n * Tp; % Wektor chwil czasu rzeczywistego

% Definicja sygnałów
e = sigma * randn(1, N); % Sygnał biały szum
x = sin(2 * pi * 5 * tn); % Sygnał x(nTp)
y = x + e; % Sygnał y(nTp)

% Definicja filtru H(q^-1)
H = tf([0.1], [1 -0.9], Tp);

% Odpowiedź filtru na wymuszenie e
v = lsim(H, e, tn);

% Wyświetlanie wykresów
figure;
subplot(4, 1, 1);
plot(tn, e);
title('Sygnał e(nTp)');
xlabel('Czas [s]');
ylabel('Amplituda');

subplot(4, 1, 2);
plot(tn, x);
title('Sygnał x(nTp)');
xlabel('Czas [s]');
ylabel('Amplituda');

subplot(4, 1, 3);
plot(tn, y);
title('Sygnał y(nTp) = x(nTp) + e(nTp)');
xlabel('Czas [s]');
ylabel('Amplituda');

subplot(4, 1, 4);
plot(tn, v);
title('Sygnał v(nTp) = H(q^-1) * e(nTp)');
xlabel('Czas [s]');
ylabel('Amplituda');


% DRUGA KROPKA

% Zakładając, że mamy już zdefiniowane sygnały e, x, y, v oraz wektor czasu tn

% Autokorelacja dla e
[Re_biased, lags] = xcorr(e, 'biased');
Re_unbiased = xcorr(e, 'unbiased');

% Podobnie dla x
[Rx_biased, ~] = xcorr(x, 'biased');
Rx_unbiased = xcorr(x, 'unbiased');

% Dla y
[Ry_biased, ~] = xcorr(y, 'biased');
Ry_unbiased = xcorr(y, 'unbiased');

% Dla v
[Rv_biased, ~] = xcorr(v, 'biased');
Rv_unbiased = xcorr(v, 'unbiased');


% Wykreślenie autokorelacji dla sygnału e
figure;
subplot(2,1,1);
plot(lags, Re_biased);
title('Autokorelacja sygnału e - estymator obciążony');
xlabel('Przesunięcie i');
ylabel('Autokorelacja');

subplot(2,1,2);
plot(lags, Re_unbiased);
title('Autokorelacja sygnału e - estymator nieobciążony');
xlabel('Przesunięcie i');
ylabel('Autokorelacja');

figure;
subplot(2,1,1);
plot(lags, Rx_biased);
title('Autokorelacja sygnału x - estymator obciążony');
xlabel('Przesunięcie i');
ylabel('Autokorelacja');

subplot(2,1,2);
plot(lags, Rx_unbiased);
title('Autokorelacja sygnału x - estymator nieobciążony');
xlabel('Przesunięcie i');
ylabel('Autokorelacja');

figure;
subplot(2,1,1);
plot(lags, Ry_biased);
title('Autokorelacja sygnału y - estymator obciążony');
xlabel('Przesunięcie i');
ylabel('Autokorelacja');

subplot(2,1,2);
plot(lags, Ry_unbiased);
title('Autokorelacja sygnału y - estymator nieobciążony');
xlabel('Przesunięcie i');
ylabel('Autokorelacja');

figure;
subplot(2,1,1);
plot(lags, Rv_biased);
title('Autokorelacja sygnału v - estymator obciążony');
xlabel('Przesunięcie i');
ylabel('Autokorelacja');

subplot(2,1,2);
plot(lags, Rv_unbiased);
title('Autokorelacja sygnału v - estymator nieobciążony');
xlabel('Przesunięcie i');
ylabel('Autokorelacja');
% Obliczanie autokorelacji
lags = [-n:n]
[R_ee, lags] = xcorr(e, 'unbiased');
[R_xx, ~] = xcorr(x, 'unbiased');
[R_yy, ~] = xcorr(y, 'unbiased');
[R_vv, ~] = xcorr(v, 'unbiased');

% Wybór zakresu [0; N-1] z wyników autokorelacji
zeroLagIndex = find(lags == 0); % Znajdź indeks odpowiadający zerowemu opóźnieniu
R_ee_positive = R_ee(zeroLagIndex:end); % Autokorelacja e dla opóźnień [0; N-1]
R_xx_positive = R_xx(zeroLagIndex:end); % Autokorelacja x dla opóźnień [0; N-1]
R_yy_positive = R_yy(zeroLagIndex:end); % Autokorelacja y dla opóźnień [0; N-1]
R_vv_positive = R_vv(zeroLagIndex:end); % Autokorelacja v dla opóźnień [0; N-1]
%{
% Wyświetlanie wyników (opcjonalnie)
figure;
subplot(4, 1, 1);
plot(lags(zeroLagIndex:end), R_ee_positive);
title('Autokorelacja e(nTp)');

subplot(4, 1, 2);
plot(lags(zeroLagIndex:end), R_xx_positive);
title('Autokorelacja x(nTp)');

subplot(4, 1, 3);
plot(lags(zeroLagIndex:end), R_yy_positive);
title('Autokorelacja y(nTp)');

subplot(4, 1, 4);
plot(lags(zeroLagIndex:end), R_vv_positive);
title('Autokorelacja v(nTp)');
%}
% Obliczanie i rysowanie funkcji autokorelacji
signals = {e, x, y, v};
signalNames = {'e(nTp)', 'x(nTp)', 'y(nTp)', 'v(nTp)'};

% czwarta kropka
figure;
for i = 1:length(signals)
    [acor, lags] = xcorr(signals{i}, 'biased'); % Obliczenie autokorelacji
    subplot(4, 1, i);
    plot(lags, acor);
    title(['Autokorelacja sygnału ' signalNames{i}]);
    xlabel('Przesunięcie');
    ylabel('Autokorelacja');
end

% PIATA KROPKA
% Obliczenie autokorelacji dla e i v
[Re_biased, lags] = xcorr(e, 'biased');
Rv_biased = xcorr(v, 'biased');

% Wykreślenie estymaty funkcji autokorelacji dla e
figure;
subplot(2,1,1);
plot(lags, Re_biased);
title('Autokorelacja sygnału e(nTp) - estymator obciążony');
xlabel('Przesunięcie i');
ylabel('Autokorelacja');
xlim([-N+1, N-1]); % Ograniczenie osi X dla lepszej czytelności

% Wykreślenie estymaty funkcji autokorelacji dla v
subplot(2,1,2);
plot(lags, Rv_biased);
title('Autokorelacja sygnału v(nTp) - estymator obciążony');
xlabel('Przesunięcie i');
ylabel('Autokorelacja');
xlim([-N+1, N-1]); % Ograniczenie osi X dla lepszej czytelności


% SZÓSTA KROPKA
% Obliczenie funkcji korelacji wzajemnej między y i x
[Ryx, lags] = xcorr(y, x, 'biased'); % Używam estymatora obciążonego dla spójności

% Wykreślenie funkcji korelacji wzajemnej
figure;
plot(lags, Ryx);
title('Korelacja wzajemna pomiędzy y(nTp) i x(nTp)');
xlabel('Przesunięcie i');
ylabel('Korelacja wzajemna');
grid on; % Dodaje siatkę dla lepszej czytelności

