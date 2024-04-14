clear all
clc
% Wczytanie danych
load('StochasticProcess.mat'); % Zakładamy, że StochasticProcess jest macierzą danych

% Wybór kilku realizacji do wykreślenia
time = StochasticProcess(1, :); % Czas jest pierwszym wierszem
realizations = StochasticProcess(2:end, :); % Realizacje są w pozostałych wierszach

% Tworzenie subplotów dla pierwszych pięciu realizacji
figure;
for i = 1:5 % Pętla dla pierwszych pięciu realizacji
    subplot(5, 1, i); % Tworzy subplot w pięciu wierszach i jednej kolumnie
    plot(time, realizations(i, :));
    title(sprintf('Realizacja %d', i));
    xlabel('Czas');
    ylabel('Amplituda');
    if i == 1
        % legend('Realizacja 1', 'Location', 'Best');
    elseif i == 5
        xlabel('Czas');
    end
end

% Obliczenie estymat dla każdej realizacji
m_i = mean(realizations,2); % Średnie dla każdej realizacji
sigma2_i = var(realizations, 0, 2); % Wariancje dla każdej realizacji

% Obliczenie estymat po czasie
m = mean(realizations); % Średnie dla każdego momentu czasu
sigma2 = var(realizations); % Wariancje dla każdego momentu czasu

% Średnie wartości estymat
mean_m_i = mean(m_i);
mean_sigma2_i = mean(sigma2_i);
mean_m = mean(m);
mean_sigma2 = mean(sigma2);

figure;
subplot(2,1,1)
plot(1:500,m_i)
title('Estymaty średniej \mu_i i \mu');
xlabel('Realizacja');
ylabel('Średnia');
legend('\mu_i', '\mu', 'Location', 'Best');
subplot(2,1,2)
plot(time,m)
title('Estymaty średniej \mu_i i \mu');
xlabel('Czas');
ylabel('Średnia');
legend('\mu_i', '\mu', 'Location', 'Best');

figure;
subplot(2,1,1)
plot(1:500,sigma2_i)
title('Estymaty wariancji \sigma^2_i i \sigma^2');
xlabel('Realizacja');
ylabel('Wariancja');
legend('\sigma^2_i', '\sigma^2', 'Location', 'Best');
subplot(2,1,2)
plot(time,sigma2)
title('Estymaty wariancji \sigma^2_i i \sigma^2');
xlabel('Czas');
ylabel('Wariancja');
legend('\sigma^2_i', '\sigma^2', 'Location', 'Best');



