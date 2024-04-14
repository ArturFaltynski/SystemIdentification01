clc
close all

Tp=0.001;
N=2000;
n=[0:N-1];
sigma_2=[0.64, 1.0];
sigma=sqrt(sigma_2);
tn=n*Tp;

x=sin(2*pi()*5*n*Tp)+0.5*sin(2*pi()*10*n*Tp)+0.25*sin(2*pi()*30*n*Tp);
e1=sigma(1)*randn(1,N);
e2=sigma(2)*randn(1,N);
H=tf([0.1],[1 -0.9],Tp);
v1=lsim(H,e1,tn);
v2=lsim(H,e2,tn);

%pkt1
figure
subplot(3,2,1)
plot(tn,x);
title('x(nTp)')
subplot(3,2,3)
plot(tn,e1);
title('e1(nTp) dla sigma=0.64')
subplot(3,2,4)
plot(tn,e2);
title('e2(nTp) dla sigma=1')
subplot(3,2,5)
plot(tn,v1);
title('v1(nTp) dla sigma=0.64')
subplot(3,2,6)
plot(tn,v2);
title('v2(nTp) dla sigma=1')
%suptitle('Wykresy wszystkich sygna³ów w dziedzinie czasu dyskretnego')

%pkt2
figure
xn=2/(N*Tp).*abs(Tp*fft(x));
stem(tn,xn)
hold on
plot(tn,x)
legend('|X_N (jw_k)|','x(nTp)');
title('dyskretne widmo amplitudowe |X_N (jw_k)|')

%pkt3
N=[1000,200,100];
n=[0:N(1)-1];
x1=sin(2*pi()*5*n*Tp)+0.5*sin(2*pi()*10*n*Tp)+0.25*sin(2*pi()*30*n*Tp);
xn1=2/(N(1))*abs(Tp*fft(x1));
tn1=n*Tp;
n=[0:N(2)-1];
x2=sin(2*pi()*5*n*Tp)+0.5*sin(2*pi()*10*n*Tp)+0.25*sin(2*pi()*30*n*Tp);
xn2=2/(N(2))*abs(Tp*fft(x2));
tn2=n*Tp;
n=[0:N(3)-1];
x3=sin(2*pi()*5*n*Tp)+0.5*sin(2*pi()*10*n*Tp)+0.25*sin(2*pi()*30*n*Tp);
xn3=2/(N(3))*abs(Tp*fft(x3));
tn3=n*Tp;
figure
subplot(3,1,1)
stem(tn1,xn1);
title('|X_N (jw_k)| dla N=1000');
subplot(3,1,2)
stem(tn2,xn2);
title('|X_N (jw_k)| dla N=200');
subplot(3,1,3)
stem(tn3,xn3);
title('|X_N (jw_k)| dla N=100');

%pkt4
N=2000;

eps_n1=Tp*sum(power(x,2))

k=[0:N-1];
omega_k=(2*pi()*k)/N;
xn_jomegak=abs(Tp*fft(x));
eps_n2=sum(1/(N*Tp)*power(abs(xn_jomegak),2))

%phi_xx=1/(N*Tp)*power(abs(xn_jomegak),2);
%eps_n3=sum(phi_xx)

%pkt5/pkt6
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
title('Estymata gêstoci widmowej mocy e1(nTp)');

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
%xlim([0 50])
legend('N', 'N/2', 'N/4');
title('Estymata gêstoci widmowej mocy e2(nTp)');

%pkt7
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
title('Estymata gêstoci widmowej mocy v1(nTp)');

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
title('Estymata gêstoci widmowej mocy v2(nTp)');
%xlim([0 50])