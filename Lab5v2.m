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
n=512; %відліки
p=3.1416;
t=(0:(n-1))/fs; % вектор часу
h=t*2*p*f0; %параметри сигналу
s =A*sin(h);
v = sqrt(D)*randn(1,length(t));% параметри шуму  
x=s+v;  %сигнал+шум

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
fprintf('Відношення сигнал/шум = %4.3g\n',snr(s,v))

%=== Завдання #1.3-1.4 ===
% Обчислення незміщеної та не зміщеної оцінки АКФ змодельованого процесу
akf_biased_02 = xcorr(x,x,fix(0.2*length(x)),'unbiased');%оцінки не зміщені
akf_biased_09 = xcorr(x,x,fix(0.9*length(x)),'unbiased');%оцінки не зміщені
akf_unbiased_02 = xcorr(x,x,fix(0.2*length(x)),'biased');%оцінки зміщені
akf_unbiased_09 = xcorr(x,x,fix(0.9*length(x)),'biased');%оцінки зміщені
figure(2) %графіки отриманих оцінок
subplot (4, 1, 1); plot (akf_biased_02); grid on;
title('Графік зміщеної оцінки сигналу(maxlag = 0.2)');
ylabel('Значення'); xlabel('Відліки');
subplot (4, 1, 2); plot (akf_unbiased_02); grid on;
title('Графік не зміщеної оцінки сигналу(maxlag = 0.9)');
ylabel('Значення'); xlabel('Відліки');
subplot (4, 1, 3); plot (akf_biased_09); grid on
title('Графік зміщеної оцінки сигналу(maxlag = 0.2)');
ylabel('Значення'); xlabel('Відліки');
subplot (4, 1, 4); plot (akf_unbiased_09); grid on
title('Графік не зміщеної оцінки сигналу(maxlag = 0.9)');
ylabel('Значення'); xlabel('Відліки');

%=== Завдання #1.5 ===
% Обчислення оцінки АКФ змодельованого процесу при збільшенні тривалості процесу
t1 = 0:100;% вектор часу 1
t2 = 0:1000; % вектор часу 2
h1=t1*2*p*f0;%параметри сигналу
s1 =sin(h1);
v1 = sqrt(D)*randn(1,length(t1));  %параметри шуму
x1=s1+v1; %модельований сигнал з тривалістю 100  
h2=t2*2*p*f0;
s2 =sin(h2);
v2 = sqrt(D)*randn(1,length(t2));  
x2=s2+v2; %модельований сигнал з тривалістю 1000  
akf_un=xcorr(x1,x1,10, 'unbiased'); %оцінки
akf_bias=xcorr(x2,x2,100, 'unbiased');%оцінки
figure(3) %Графіки отриманих оцінок
subplot (2, 1, 1); plot (akf_un); grid on;
title('Графік оцінки сигналу(тривалість 100)');
ylabel('Значення'); xlabel('Відліки');
subplot (2, 1, 2); plot (akf_bias); grid on;
title('Графік оцінки сигналу(тривалість 1000)');
ylabel('Значення'); xlabel('Відліки');

%=== Завдання #2.1 ===
% Завантаження сигналу ЕЕГ файл (eeg1-p4.dat) 
eeg_p4 =load('eeg1-p4.dat'); % сигнал EEГ
eeg1_p4= detrend(eeg_p4); %центрування сигналу
figure(4)
plot(eeg1_p4);
title ('Сигнал ЕЕГ')
ylabel('Амплітуда'); xlabel('Відліки');

%=== Завдання #2.2 ===
% Виділення епохи сигналу ЕЕГ від t1 = 4,7 с до t2 = 5,8 с
fs=100;% xfcnjnf lbcrhtnbpfws]
t = (0:length(eeg1_p4)-1)/fs; % вектор часу
n1 = (4.7*fs) + 1; 
n2 = (5.8*fs) + 1;
figure(5)
plot(t(n1:n2),eeg1_p4(n1:n2));
title ('Епоха сигналу від t1 = 4,7 с до t2 = 5,8 с')
ylabel('Амплітуда'); xlabel('Відліки');

%== Завдання #2.3 ===
% Обчислення незміщеної оцінки АКФ сигналу ЕЕГ
akf = xcorr(eeg1_p4,eeg1_p4, fix(0.9*length(eeg1_p4)), 'unbiased'); %обчислення оцынки
figure(6)
plot(akf);
title('Графік не зміщеної оцінки епохи сигналу eeg1-p4.dat');
ylabel('Значення'); xlabel('Відліки');

%=== Завдання #2.4 ===
% Обчислення спектральної щільності сигналу
Sx = abs(fft(eeg1_p4)/length(eeg1_p4)); %обчислення модуля швидкого перетворення Фурьє
B=Sx(1:length(eeg1_p4)/2+1);%Отримання односторонньої спектральної щільності 
B(2:end-1)=2*B(2:end-1);
f =fs*(0:(length(eeg1_p4)/2))/length(eeg1_p4);
figure(7)
plot(f,B), grid on;
title('Спектральна щільність');
ylabel('Sxx'); xlabel('Відліки');

%=== Завдання #2.5 ===
% Завантаження сигналу ЕЕГ файл (eeg1-f3.dat)
eeg_f3 = load('eeg1-f3.dat'); %завантаження сигналу
eeg1_f3 = detrend(eeg_f3); %центрування 
t = (0:length(eeg1_f3)-1)/fs; % вектор часу
figure(8)
plot(t,eeg1_f3)%побудова графіку
title ('Сигнал ЕЕГ eeg1-f3.dat')
ylabel('Амплітуда'); xlabel('Відліки');

% Виділення епохи сигналу ЕЕГ від t1 = 4,2 с до t2 = 4,96 с
t_f3 = (0:length(eeg1_f3)-1)/fs;
n1 = (4.2*fs) + 1; %формування відліків
n2 = (4.96*fs) + 1; %формування відліків
figure(8)
plot(t_f3(n1:n2),eeg1_f3(n1:n2));
title ('Епоха сигналу від t1 = 4,2 с до t2 = 4,96 с')
ylabel('Амплітуда'); xlabel('Відліки');

% Обчислення незміщеної оцінки АКФ сигналу ЕЕГ
akf2 = xcorr(eeg1_f3,eeg1_f3,fix(0.2*length(eeg1_f3)) , 'unbiased'); %обчислення оцінки
figure(10)
plot(akf2);
title('Графік не зміщеної оцінки сигналу ЕЕГ сигналу');
ylabel('Значення'); xlabel('Відліки');

% Обчислення спектральної щільності сигналу
Sxx1 = abs(fft(eeg1_f3))/length(eeg1_f3);% обчислення модуля спектральної щільності
Bx1=Sxx1(1:length(eeg1_f3)/2+1);%Отримання односторонньої спектральної щільності 
Bx1(2:end-1)=2*Bx1(2:end-1);
f =fs*(0:(length(eeg1_f3)/2))/length(eeg1_f3);
figure(11)
plot(f,Bx1), grid on;
title('Спектральна щільність');
ylabel('Sxx1'); xlabel('Відліки');

% Завантаження сигналу  ЕЕГ файл (eeg1-o2.dat)
% Завантаження сигналу ЕЕГ файл (eeg1-c3.dat)
eeg_o2 = load('eeg1-o2.dat');  % Завантаження сигналу  ЕЕГ файл (eeg1-o2.dat)
eeg1_o2 = detrend(eeg_o2);  %центрування 
eeg_c3 = load('eeg1-c3.dat'); % Завантаження сигналу ЕЕГ файл (eeg1-c3.dat)
eeg1_c3 = detrend(eeg_c3); %центрування 
t_o2 = (0:length(eeg1_o2)-1)/fs; % вектор часу
t_c3 = (0:length(eeg1_c3)-1)/fs; % вектор часу
figure(12)
subplot(2,1,1); plot(t_o2,eeg1_o2);
title ('Сигнал ЕЕГ eeg1-o2.dat')
ylabel('Амплітуда'); xlabel('Відліки');
subplot(2,1,2); plot(t_c3,eeg1_c3);
title ('Сигнал ЕЕГ eeg1_c3')
ylabel('Амплітуда'); xlabel('Відліки');

%=== Завдання #3.2 ===
% Виділення епохи сигналів ЕЕГ від t1 = 5,71 с до t2 = 6,78 с
t_o2 = (0:length(eeg1_o2)-1)/fs; % вектор часу
t_c3 = (0:length(eeg1_c3)-1)/fs; % вектор часу
n1 = (5.71*fs) + 1;   %формування відліків
n2 = (6.78*fs) + 1;   %формування відліків
figure(13)
subplot (2, 1, 1); plot(t_o2(n1:n2),eeg1_o2(n1:n2));
title ('Епоха сигналу eeg1_o2 від t1 = 4,2 с до t2 = 4,96 с')
ylabel('Амплітуда'); xlabel('Відліки');
subplot (2, 1, 2); plot(t_c3(n1:n2),eeg1_c3(n1:n2));
title ('Епоха сигналу eeg1_c3 від t1 = 4,2 с до t2 = 4,96 с')
ylabel('Амплітуда'); xlabel('Відліки');

%=== Завдання #3.3 ===
% Обчислення ВКФ сигналів ЕЕГ
vkf1 = xcorr(eeg1_o2,eeg1_c3,'biased'); %обчислення ВКФ
figure(14)
plot(vkf1);
title ('ВКФ сигналів ЕЕГ eeg1_o2 та eeg1_c3')
ylabel('Значення'); xlabel('Відліки');

%=== Завдання #3.4 ===
% Обчислення взаємної спектральної щільності сигналів ЕЕГ
Sxy = abs(fft(vkf1)/length(vkf1));  % обчислення модуля спектральної щільності
Bz=Sxy(1:length(vkf1)/2+1);%Отримання односторонньої спектральної щільності 
Bz(2:end-1)=2*Bz(2:end-1);
f =fs*(0:(length(vkf1)/2))/length(vkf1);
figure(15)
plot(f,Bz), grid on;
title ('Взаємна спектральна щільность заданих сигналів ЕЕГ')
ylabel('Bz'); xlabel('Частота');

%=== Завдання #3.5 ===
% Завантаження сигналу ЕЕГ файл (eeg1-p3.dat) та сигналу (eeg1-p4.dat)
eeg1_p3 = load('eeg1-p3.dat'); % завантаження
eeg1_p3 = detrend(eeg1_p3);
eeg1_p4 = load('eeg1-p4.dat');  % завантаження
eeg1_p4 = detrend(eeg1_p4);
t_p3 = (0:length(eeg1_p3)-1)/fs; % вектор часу
t_p4 = (0:length(eeg1_p4)-1)/fs; % вектор часу
figure(16)
subplot (2, 1, 1); plot(t_p3, eeg1_p3);
title ('Графік сигналу eeg1-p3.dat')
subplot (2, 1, 2); plot(t_p4, eeg1_p4);
title ('Графік сигналу eeg1-p4.dat')

% Виділення епохи сигналів ЕЕГ від t1 = 4,7 с до t2 = 5,8 с
t_p3 = (0:length(eeg1_p3)-1)/fs;
t_p4 = (0:length(eeg1_p4)-1)/fs;
n1 = (4.7*fs) + 1; 
n2 = (5.8*fs) + 1;
figure(17)
subplot (2, 1, 1); plot(t_p3(n1:n2),eeg1_p3(n1:n2));
title ('Епоха сигналу eeg1_p3 від 4,7с до 5,8с')
ylabel('Амплітуда'); xlabel('Відліки');
subplot (2, 1, 2); plot(t_p4(n1:n2),eeg1_p4(n1:n2));
title ('Епоха сигналу eeg1-p4 від 4,7с до 5,8с')
ylabel('Амплітуда'); xlabel('Відліки');
% Обчислення ВКФ сигналів ЕЕГ
vkf2 = xcorr(eeg1_p3,eeg1_p4); 
figure(18)
plot(vkf2);
title ('% Обчислення ВКФ сигналів ЕЕГ eeg1_p3 та eeg1_p4')
ylabel('Значення'); xlabel('Відліки');

% Обчислення взаємної спектральної щільності сигналів ЕЕГ
Svkf2 = abs(fft(vkf2)/length(vkf2)); %обчислення модуля швидкого перетворення Фурьє
Bn=Svkf2(1:length(vkf2)/2+1);%Отримання односторонньої спектральної щільності 
Bn(2:end-1)=2*Bn(2:end-1);
f =fs*(0:(length(vkf2)/2))/length(vkf2);
figure(19)
plot(f,Bn), grid on;
title ('Взаємна спектральна щільность сигналів eeg1_p3 та eeg1-p4')
ylabel('Значення'); xlabel('Частота');

%=== Завдання #3.6 ===
% Завантаження сигналу ЕЕГ файл (eeg1-f3.dat)
% Завантаження сигналу ЕЕГ файл (eeg1-f4.dat)
eeg1_f3 = load('eeg1-f3.dat');
eeg1_f3 = detrend(eeg1_f3);
eeg1_f4 = load('eeg1-f4.dat');
eeg1_f4 = detrend(eeg1_f4);
t_f3 = (0:length(eeg1_f3)-1)/fs;
t_f4 = (0:length(eeg1_f4)-1)/fs;
figure(20)
subplot (2, 1, 1); plot(t_f3, eeg1_f3);
title ('Графік сигналу eeg1-f3.dat')
ylabel('Амплітуда'); xlabel('Відліки');
subplot (2, 1, 2); plot(t_f4, eeg1_f4);
title ('Графік сигналу eeg1-f4.dat')
ylabel('Амплітуда'); xlabel('Відліки');

% Виділення епохи сигналів ЕЕГ від t1 = 4,13 с до t2 = 4,96 с t3 = 1,43 с до t4 = 2,26 с
n1=(4.13*fs) + 1;%формування відліків
n2=(4.96*fs) + 1;%формування відліків
n3=(1.43*fs) + 1;%формування відліків
n4=(2.26*fs) + 1;%формування відліків
figure(21)
subplot (2, 1, 1); plot(t_f3(n1:n2),eeg1_f3(n1:n2));
title ('Епоха сигналу eeg1_f3 від 4,13с до 4.96с')
ylabel('Амплітуда'); xlabel('Відліки');
subplot (2, 1, 2); plot(t_f4(n3:n4),eeg1_f4(n3:n4));
title ('Епоха сигналу eeg1_f4 від 1.43с до 2.26с')
ylabel('Амплітуда'); xlabel('Відліки');

% Обчислення ВКФ сигналів ЕЕГ
vkf3 = xcorr(eeg1_f3, eeg1_f4, 250); % обчислення ВКФ
figure(22)
plot(vkf3);
title ('ВКФ епох сигналів eeg1_f3 та eeg1_f4')
ylabel('Значення'); xlabel('Відліки');

% Обчислення взаємної спектральної щільності сигналів ЕЕГ
Sx = abs(fft(vkf3)/length(vkf3)); %обчислення модуля швидкого перетворення Фурьє
Bq=Sx(1:length(vkf3)/2+1);%Отримання односторонньої спектральної щільності 
Bq(2:end-1)=2*Bq(2:end-1);
f =fs*(0:(length(vkf3)/2))/length(vkf3);
figure(23)
plot(f,Bq), grid on;
title ('Взаємна спектральна щільность епох сигналів eeg1_p3 та eeg1-p4')
ylabel('Значення'); xlabel('частота');



