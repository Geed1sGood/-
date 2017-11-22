clear
clc
close all

disp('����������� ������ #5')
disp('������ֲ���� ���˲� ��������˲�')
disp('�������:  , ����� ... �Ͳ���')

%=== �������� #1.1 ===
% ����������� ������������� ����������� �������

fs = 20;
f0 = 5 ;
A = 0.5;
D = 0.25; 
n = 512; 
p=3.1416;
t = [0:0.02:5];		% ������ ����
h=t*2*p*f0;
s =sin(h);
v = D*randn(1,length(t));  
x=s+v;  

% ������� ��������� �������

figure(1)
subplot (3, 1, 1); plot(t, v), grid on;
title ('������ x(t)');
ylabel('��������');
xlabel('³����');
subplot (3, 1, 2); plot (t, s); grid on;
title ('������ s(t)');
ylabel('��������');
xlabel('³����');
subplot (3, 1, 3); plot (t, x); grid on;
title ('������ y(t)');
ylabel('��������');
xlabel('³����');

%=== �������� #1.2 ===
% ���������� ������ ������� ����, ������� ������� �� ��������� ���/������

disp(var(v))
disp(var(x))
l=snr(s,v)

%=== �������� #1.3 ===
% ���������� �������� ������ ��� �������������� �������

t1=-5:0.02:5;
maxlag1=1;
maxlag2=50;
maxlag3=100;
ran = xcorr(x,x,'unbiased');
rlag12=xcorr(x,x,maxlag1, 'unbiased');
rlag13=xcorr(x,x,maxlag2, 'unbiased');
rlag14=xcorr(x,x,maxlag3, 'unbiased');
figure(2)
subplot (3, 1, 1); plot (rlag12); grid on;
subplot (3, 1, 2); plot (rlag13); grid on;
subplot (3, 1, 3); plot (rlag14); grid on;
figure(3)
rlag15=xcorr(x,x,maxlag1+6, 'biased');
rlag16=xcorr(x,x,maxlag2+6, 'biased');
rlag17=xcorr(x,x,maxlag3+6, 'biased');
subplot (3, 1, 1); plot (rlag15); grid on;
subplot (3, 1, 2); plot (rlag16); grid on;
subplot (3, 1, 3); plot (rlag17); grid on;

fs = 202;
f0 = 5 ;
A = 0.5;
D = 0.25; 
n = 512; 
p=3.1416
t = [0:0.02:5];		% ������ ����
h=t*2*p*f0;
s =sin(h);
v = D*randn(1,length(t));  
x=s+v;  

% ������� ��������� �������

figure(4)
subplot (3, 1, 1); plot(t, v), grid on;
title ('������ x(t)');
ylabel('��������');
xlabel('³����');
subplot (3, 1, 2); plot (t, s); grid on;
title ('������ s(t)');
ylabel('��������');
xlabel('³����');
subplot (3, 1, 3); plot (t, x); grid on;
title ('������ y(t)');
ylabel('��������');
xlabel('³����');

%=== �������� #1.2 ===

maxlag1=fix(0.2*length(x));
maxlag2=100;
maxlag3=3;
ran = xcorr(x,x,'unbiased');
rlag12=xcorr(x,x,maxlag1, 'unbiased');
rlag13=xcorr(x,x,maxlag2, 'unbiased');
rlag14=xcorr(x,x,maxlag3, 'unbiased');
figure(5)

subplot (3, 1, 1); plot (rlag12); grid on;
subplot (3, 1, 2); plot (rlag13); grid on;
subplot (3, 1, 3); plot (rlag14); grid on;
max(rlag15);
figure(6)
rlag15=xcorr(x,x,maxlag1+6, 'biased');
rlag16=xcorr(x,x,maxlag2+6, 'biased');
rlag17=xcorr(x,x,maxlag3+6, 'biased');
subplot (3, 1, 1); plot (rlag15); grid on;
subplot (3, 1, 2); plot (rlag16); grid on;
subplot (3, 1, 3); plot (rlag17); grid on;

%=== �������� #2.1 ===
% ������������ ������� ��� ���� (eeg1-p4.dat)
fs =100; 
eeg1_p4 =load('eeg1-p4.dat'); % ������ EEU
eeg1 = detrend(eeg1_p4);
figure(7)   
plot(eeg1);

%=== �������� #2.2 ===
% �������� ����� ������� ��� �� t1 = 4,7 � �� t2 = 5,8 �

fs=100;
t1 =4.7; 
t2 =5.8;
t = (0:length(eeg1_p4)-1)/fs;
n1 = (t1*fs) + 1; 
n2 = (t2*fs) + 1;
figure(8)
plot(t(n1:n2),eeg1_p4(n1:n2));

%=== �������� #2.3 ===
% ���������� �������� ������ ��� ������� ���
maxlag3 =fix(0.9*length(eeg1_p4))
akf = xcorr(eeg1_p4,eeg1_p4, maxlag3, 'unbiased');
figure(8)
plot(akf);
max(akf);

%=== �������� #2.4 ===
% ���������� ����������� �������� �������
Sxx = abs(fft(eeg1_p4))/length(eeg1_p4);
A=Sxx';
B=[A(1),2.*A(2:end-1),A(end)];
Y=fftshift(B);
N = length(eeg1_p4);
f = (0:N-1)/N*fs;
figure(9)
plot(f,Y), grid on;

%=== �������� #2.5 ===
% ������������ ������� ��� ���� (eeg1-f3.dat)
fs = 100;
eeg1_f3 = load('eeg1-f3.dat');
eeg1_f3 = detrend(eeg1_f3);
t = (0:length(eeg1_f3)-1)/fs;
figure(10)
plot(t,eeg1_f3)
 
% �������� ����� ������� ��� �� t1 = 4,2 � �� t2 = 4,96 �
fs=100;
t1 =4.2; 
t2 =4.96;
t = (0:length(eeg1_f3)-1)/fs;
n1 = (t1*fs) + 1; 
n2 = (t2*fs) + 1;
figure(11)
plot(t(n1:n2),eeg1_f3(n1:n2));

% ���������� �������� ������ ��� ������� ���
maxlag4 =20;
akf2 = xcorr(eeg1_f3,eeg1_f3, maxlag4, 'unbiased');
figure(12)
plot(akf2);
max(akf2);

% ���������� ����������� �������� �������
Sxx1 = abs(fft(eeg1_f3))/length(eeg1_f3);
C = Sxx1';
Z = fftshift(C);
N2 = length(eeg1_f3);
f = (0:N2-1)/N2*fs;
figure(13)
plot(f,Z), grid on;

%=== �������� #3.1 ===
% ������������ ������� ��� ���� (eeg1-o2.dat)
fs = 100;
eeg1_o2 = load('eeg1-o2.dat');
eeg1_o2 = detrend(eeg1_o2);
t = (0:length(eeg1_o2)-1)/fs;
figure(14)
plot(t,eeg1_o2);

% ������������ ������� ��� ���� (eeg1-c3.dat)
fs = 100;
eeg1_c3 = load('eeg1-c3.dat');
eeg1_c3 = detrend(eeg1_c3);
t = (0:length(eeg1_c3)-1)/fs;
figure(15)
plot(t,eeg1_c3);

%=== �������� #3.2 ===
% �������� ����� ������� ��� �� t1 = 5,71 � �� t2 = 6,78 �

fs=100;
t1 =5.71; 
t2 =6.78;
t11 = (0:length(eeg1_o2)-1)/fs;
t22 = (0:length(eeg1_c3)-1)/fs;
n1 = (t1*fs) + 1; 
n2 = (t2*fs) + 1;
figure(16)
subplot (2, 1, 1); plot(t11(n1:n2),eeg1_o2(n1:n2));
subplot (2, 1, 2); plot(t22(n1:n2),eeg1_c3(n1:n2));

%=== �������� #3.3 ===
% ���������� ��� ������� ���
maxlag1=1;
vkf1 = xcorr(eeg1_o2,eeg1_c3,'biased');
figure(17)
plot(vkf1);

%=== �������� #3.4 ===
% ���������� ������ ����������� �������� ������� ���
Sxy = abs(fft(vkf1)/length(vkf1));
Axy = Sxy';
Bxy = fftshift(Axy);
N = length(vkf1);
f = (0:N-1)/N*fs;
figure(18)
plot(f,Bxy), grid on;

%=== �������� #3.5 ===
% ������������ ������� ��� ���� (eeg1-p3.dat) �� ������� (eeg1-p4.dat)
fs = 100;
eeg1_p3 = load('eeg1-p3.dat');
eeg1_p3 = detrend(eeg1_p3);
eeg1_p4 = load('eeg1-p4.dat');
eeg1_p4 = detrend(eeg1_p4);
t_p3 = (0:length(eeg1_p3)-1)/fs;
t_p4 = (0:length(eeg1_p4)-1)/fs;
figure(19)
subplot (2, 1, 1); plot(t_p3, eeg1_p3);
subplot (2, 1, 2); plot(t_p4, eeg1_p4);

% �������� ����� ������� ��� �� t1 = 4,7 � �� t2 = 5,8 �
fs=100;
t1 =4,7; 
t2 =5,8;
t_p3 = (0:length(eeg1_p3)-1)/fs;
t_p4 = (0:length(eeg1_p4)-1)/fs;
n1 = (t1*fs) + 1; 
n2 = (t2*fs) + 1;
figure(20)
subplot (2, 1, 1); plot(t11(n1:n2),eeg1_p3(n1:n2));
subplot (2, 1, 2); plot(t22(n1:n2),eeg1_p4(n1:n2));

% ���������� ��� ������� ���
maxlag1=1;
vkf2 = xcorr(eeg1_o2,eeg1_c3,'biased');
figure(21)
plot(vkf2);

% ���������� ������ ����������� �������� ������� ���
Svkf2 = abs(fft(vkf2)/length(vkf2));
Axy2 = Svkf2';
Bxy2 = fftshift(Axy2);
N = length(vkf2);
f = (0:N-1)/N*fs;
figure(22)
plot(f,Bxy2), grid on;

%=== �������� #3.6 ===
% ������������ ������� ��� ���� (eeg1-f3.dat)
% ������������ ������� ��� ���� (eeg1-f4.dat)
fs = 100;
eeg1_f3 = load('eeg1-f3.dat');
eeg1_f3 = detrend(eeg1_f3);
eeg1_f4 = load('eeg1-f4.dat');
eeg1_f4 = detrend(eeg1_f4);
t_f3 = (0:length(eeg1_f3)-1)/fs;
t_f4 = (0:length(eeg1_f4)-1)/fs;
figure(23)
subplot (2, 1, 1); plot(t_f3, eeg1_f3);
subplot (2, 1, 2); plot(t_f4, eeg1_f4);



