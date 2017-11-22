clear
clc
close all

disp('Лабораторна робота #5')
disp('КОРЕЛЯЦІЙНИЙ АНАЛІЗ БІОСИГНАЛІВ')
disp('Виконав:  , група ... ННІІДС')

%=== Завдання #1.1 ===
% Моделювання стаціонарного випадкового процесу
fs = 20; 
f0 = 5 ;
A = 0.5;
D = 0.25;
n=512;
p=3.1416;
t=(0:(n-1))/fs; % вектор часу
h=t*2*p*f0;
s =A*sin(h);
v = sqrt(D)*randn(1,length(t));  
x=s+v;  
% Графіки отриманих сигналів
figure(1)
subplot (3, 1, 1); plot(t, v), grid on;
title ('Сигнал x(t)');
ylabel('Амплітуда'); xlabel('Відліки');
subplot (3, 1, 2); plot (t, s); grid on;
title ('Сигнал s(t)');
ylabel('Амплітуда'); xlabel('Відліки');
subplot (3, 1, 3); plot (t, x); grid on;
title ('Сигнал y(t)');
ylabel('Амплітуда'); xlabel('Відліки');
%=== Завдання #1.2 ===
% Обчислення оцінки дисперсії шуму, дисперсії сигналу та відношення шум/сигнал
fprintf('Дисперсія шуму = %4.3g\n',var(v))
fprintf('Дисперсія шуму = %4.3g\n',var(v))
fprintf('Дисперсія шуму = %4.3g\n',snr(s,v))
%=== Завдання #1.3 ===
% Обчислення незміщеної та не зміщеної оцінки АКФ змодельованого процесу
maxlag_02=fix(0.2*length(x)); %максимальна затримка
maxlag_09 = fix(0.9*length(x));%максимальна затримка
akf_biased_02 = xcorr(x,x,maxlag_02,'unbiased'); %оцінки:
akf_biased_09 = xcorr(x,x,maxlag_09,'unbiased');
akf_unbiased_02 = xcorr(x,x,maxlag_02, 'biased');
akf_unbiased_09 = xcorr(x,x,maxlag_09, 'biased');
figure(2) %графіки отриманих оцінок
subplot (4, 1, 1); plot (akf_biased_02); grid on;
title('Графік зміщеної оцінки сигналу(maxlag = 0.2)');;
subplot (4, 1, 2); plot (akf_unbiased_02); grid on;
title('Графік не зміщеної оцінки сигналу(maxlag = 0.9)');
subplot (4, 1, 3); plot (akf_biased_09); grid on
title('Графік зміщеної оцінки сигналу(maxlag = 0.2)');
subplot (4, 1, 4); plot (akf_unbiased_09); grid on
title('Графік не зміщеної оцінки сигналу(maxlag = 0.9)');
%=== Завдання #1.5 ===
% Обчислення оцінки АКФ змодельованого процесу при збільшенні тривалості процесу
t1 = 100;% векторb часу
t2 = 1000;
h1=t1*2*p*f0;%параметри сигналу
s1 =sin(h1);
v1 = sqrt(D)*randn(1,length(t1));  %параметри шуму
x1=s1+v1; %модельований сигнал з тривалістю 100  
h2=t2*2*p*f0;
s2 =sin(h2);
v2 = sqrt(D)*randn(1,length(t2));  
x2=s2+v2; %модельований сигнал з тривалістю 1000  
maxlag1=10;
maxlag2=30;
akf_un=xcorr(x1,x1,maxlag1, 'unbiased');
akf_bias=xcorr(x2,x2,maxlag2, 'unbiased');
figure(3) %Графіки 
subplot (2, 1, 1); plot (akf_un); grid on;
subplot (2, 1, 2); plot (akf_bias); grid on;
%=== Завдання #2.1 ===
% Завантаження сигналу ЕЕГ файл (eeg1-p4.dat) 
eeg_p4 =load('eeg1-p4.dat'); % сигнал EEU
eeg1_p4= detrend(eeg_p4);
figure(5)
title ('Сигнал ЕЕГ')
plot(eeg1_p4);
%=== Завдання #2.2 ===
% Виділення епохи сигналу ЕЕГ від t1 = 4,7 с до t2 = 5,8 с
fs=100;
t1 =4.7; 
t2 =5.8;
t = (0:length(eeg1_p4)-1)/fs;
n1 = (t1*fs) + 1; 
n2 = (t2*fs) + 1;
figure(6)
title ('Епоха сигналу від t1 = 4,7 с до t2 = 5,8 с')
plot(t(n1:n2),eeg1_p4(n1:n2));
%== Завдання #2.3 ===
% Обчислення незміщеної оцінки АКФ сигналу ЕЕГ
maxlag3 =fix(0.9*length(eeg1_p4))
akf = xcorr(eeg1_p4,eeg1_p4, maxlag3, 'unbiased');
figure(7)
plot(akf);
max(akf);

%=== Завдання #2.4 ===
% Обчислення спектральної щільності сигналу
Sxx = abs(fft(eeg1_p4)/length(eeg1_p4));
A=Sxx';
B=[A(1),2.*A(2:end-1),A(end)];
Y=fftshift(B);
N = length(eeg1_p4);
f = (0:N-1)/N*fs;
figure(8)
plot(f,Y), grid on;

%=== Завдання #2.5 ===
% Завантаження сигналу ЕЕГ файл (eeg1-f3.dat)
fs = 100;
eeg1_f3 = load('eeg1-f3.dat');
eeg1_f3 = detrend(eeg1_f3);
t = (0:length(eeg1_f3)-1)/fs;
figure(9)
plot(t,eeg1_f3)
 
% Виділення епохи сигналу ЕЕГ від t1 = 4,2 с до t2 = 4,96 с
fs=100;
t1 =4.2; 
t2 =4.96;
t = (0:length(eeg1_f3)-1)/fs;
n1 = (t1*fs) + 1; 
n2 = (t2*fs) + 1;
figure(10)
plot(t(n1:n2),eeg1_f3(n1:n2));

% Обчислення незміщеної оцінки АКФ сигналу ЕЕГ
maxlag4 =20;
akf2 = xcorr(eeg1_f3,eeg1_f3, maxlag4, 'unbiased');
figure(11)
plot(akf2);
max(akf2);

% Обчислення спектральної щільності сигналу
Sxx1 = abs(fft(eeg1_f3))/length(eeg1_f3);
C = Sxx1';
Z = fftshift(C);
N2 = length(eeg1_f3);
f = (0:N2-1)/N2*fs;
figure(12)
plot(f,Z), grid on;

%=== Завдання #3.1 ===
% Завантаження сигналу ЕЕГ файл (eeg1-o2.dat)
fs = 100;
eeg1_o2 = load('eeg1-o2.dat');
eeg1_o2 = detrend(eeg1_o2);
t = (0:length(eeg1_o2)-1)/fs;
figure(13)
plot(t,eeg1_o2);

% Завантаження сигналу ЕЕГ файл (eeg1-c3.dat)
fs = 100;
eeg1_c3 = load('eeg1-c3.dat');
eeg1_c3 = detrend(eeg1_c3);
t = (0:length(eeg1_c3)-1)/fs;
figure(14)
plot(t,eeg1_c3);

%=== Завдання #3.2 ===
% Виділення епохи сигналів ЕЕГ від t1 = 5,71 с до t2 = 6,78 с

fs=100;
t1 =5.71; 
t2 =6.78;
t11 = (0:length(eeg1_o2)-1)/fs;
t22 = (0:length(eeg1_c3)-1)/fs;
n1 = (t1*fs) + 1; 
n2 = (t2*fs) + 1;
figure(15)
subplot (2, 1, 1); plot(t11(n1:n2),eeg1_o2(n1:n2));
subplot (2, 1, 2); plot(t22(n1:n2),eeg1_c3(n1:n2));

%=== Завдання #3.3 ===
% Обчислення ВКФ сигналів ЕЕГ
maxlag1=10;
vkf1 = xcorr(eeg1_o2,eeg1_c3, maxlag1, 'biased');
figure(16)
plot(vkf1);

%=== Завдання #3.4 ===
% Обчислення взаємної спектральної щільності сигналів ЕЕГ
Sxy = abs(fft(vkf1)/length(vkf1));
Axy = Sxy';
Bxy = fftshift(Axy);
N = length(vkf1);
f = (0:N-1)/N*fs;
figure(17)
plot(f,Bxy), grid on;

%=== Завдання #3.5 ===
% Завантаження сигналу ЕЕГ файл (eeg1-p3.dat) та сигналу (eeg1-p4.dat)
fs = 100;
eeg1_p3 = load('eeg1-p3.dat');
eeg1_p3 = detrend(eeg1_p3);
eeg1_p4 = load('eeg1-p4.dat');
eeg1_p4 = detrend(eeg1_p4);
t_p3 = (0:length(eeg1_p3)-1)/fs;
t_p4 = (0:length(eeg1_p4)-1)/fs;
figure(18)
subplot (2, 1, 1); plot(t_p3, eeg1_p3);
subplot (2, 1, 2); plot(t_p4, eeg1_p4);

% Виділення епохи сигналів ЕЕГ від t1 = 4,7 с до t2 = 5,8 с
fs=100;
t1= 4,7; 
t2= 5,8;
t_p3 = (0:length(eeg1_p3)-1)/fs;
t_p4 = (0:length(eeg1_p4)-1)/fs;
n1 = (t1*fs) + 1; 
n2 = (t2*fs) + 1;
figure(19)
subplot (2, 1, 1); plot(t11(n1:n2),eeg1_p3(n1:n2));
subplot (2, 1, 2); plot(t22(n1:n2),eeg1_p4(n1:n2));

% Обчислення ВКФ сигналів ЕЕГ
maxlag1=100;
vkf2 = xcorr(eeg1_o2,eeg1_c3, maxlag1, 'biased');
figure(20)
plot(vkf2);

% Обчислення взаємної спектральної щільності сигналів ЕЕГ
Svkf2 = abs(fft(vkf2)/length(vkf2));
Axy2 = Svkf2';
Bxy2 = fftshift(Axy2);
N = length(vkf2);
f = (0:N-1)/N*fs;
figure(21)
plot(f,Bxy2), grid on;

%=== Завдання #3.6 ===
% Завантаження сигналу ЕЕГ файл (eeg1-f3.dat)
% Завантаження сигналу ЕЕГ файл (eeg1-f4.dat)
fs = 100;
eeg1_f3 = load('eeg1-f3.dat');
eeg1_f3 = detrend(eeg1_f3);
eeg1_f4 = load('eeg1-f4.dat');
eeg1_f4 = detrend(eeg1_f4);
t_f3 = (0:length(eeg1_f3)-1)/fs;
t_f4 = (0:length(eeg1_f4)-1)/fs;
figure(22)
subplot (2, 1, 1); plot(t_f3, eeg1_f3);
subplot (2, 1, 2); plot(t_f4, eeg1_f4);

% Виділення епохи сигналів ЕЕГ від t1 = 4,13 с до t2 = 4,96 с t3 = 1,43 с до t4 = 2,26 с
fs=100;
t1 =4.13; 
t2 =4.96;
t3 =1.43;
t4 =2.26;
n1 = (t1*fs) + 1; 
n2 = (t2*fs) + 1;
n3 = (t3*fs) + 1;
n4 = (t4*fs) + 1;
figure(23)
subplot (2, 1, 1); plot(t11(n1:n2),eeg1_f3(n1:n2));
subplot (2, 1, 2); plot(t22(n3:n4),eeg1_f4(n3:n4));

% Обчислення ВКФ сигналів ЕЕГ
maxlag=100;
vkf3 = xcorr(eeg1_f3, eeg1_f4, maxlag);
figure(24)
plot(vkf3);

% Обчислення взаємної спектральної щільності сигналів ЕЕГ
Svkf3 = abs(fft(vkf3)/length(vkf3));
Axy3 = Svkf3';
Bxy3 = fftshift(Axy3);
N = length(vkf3);
f = (0:N-1)/N*fs;
figure(25)
plot(f,Bxy3), grid on;





