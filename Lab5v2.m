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
n=512; %�����
p=3.1416;
t=(0:(n-1))/fs; % ������ ����
h=t*2*p*f0; %��������� �������
s =A*sin(h);
v = sqrt(D)*randn(1,length(t));% ��������� ����  
x=s+v;  %������+���

% ������� ��������� �������
figure(1)
subplot (3, 1, 1); plot(t, v), grid on;
title ('������ x(t)');
ylabel('��������'); xlabel('³����');
subplot (3, 1, 2); plot (t, s); grid on;
title ('������ s(t)');
ylabel('��������'); xlabel('³����');
subplot (3, 1, 3); plot (t, x); grid on;
title ('������ y(t)');
ylabel('��������'); xlabel('³����');

%=== �������� #1.2 ===
% ���������� ������ ������� ����, ������� ������� �� ��������� ���/������
fprintf('�������� ���� = %4.3g\n',var(v))
fprintf('�������� ���� = %4.3g\n',var(v))
fprintf('³�������� ������/��� = %4.3g\n',snr(s,v))

%=== �������� #1.3-1.4 ===
% ���������� �������� �� �� ������ ������ ��� �������������� �������
akf_biased_02 = xcorr(x,x,fix(0.2*length(x)),'unbiased');%������ �� �����
akf_biased_09 = xcorr(x,x,fix(0.9*length(x)),'unbiased');%������ �� �����
akf_unbiased_02 = xcorr(x,x,fix(0.2*length(x)),'biased');%������ �����
akf_unbiased_09 = xcorr(x,x,fix(0.9*length(x)),'biased');%������ �����
figure(2) %������� ��������� ������
subplot (4, 1, 1); plot (akf_biased_02); grid on;
title('������ ������ ������ �������(maxlag = 0.2)');
ylabel('��������'); xlabel('³����');
subplot (4, 1, 2); plot (akf_unbiased_02); grid on;
title('������ �� ������ ������ �������(maxlag = 0.9)');
ylabel('��������'); xlabel('³����');
subplot (4, 1, 3); plot (akf_biased_09); grid on
title('������ ������ ������ �������(maxlag = 0.2)');
ylabel('��������'); xlabel('³����');
subplot (4, 1, 4); plot (akf_unbiased_09); grid on
title('������ �� ������ ������ �������(maxlag = 0.9)');
ylabel('��������'); xlabel('³����');

%=== �������� #1.5 ===
% ���������� ������ ��� �������������� ������� ��� �������� ��������� �������
t1 = 0:100;% ������ ���� 1
t2 = 0:1000; % ������ ���� 2
h1=t1*2*p*f0;%��������� �������
s1 =sin(h1);
v1 = sqrt(D)*randn(1,length(t1));  %��������� ����
x1=s1+v1; %������������ ������ � ��������� 100  
h2=t2*2*p*f0;
s2 =sin(h2);
v2 = sqrt(D)*randn(1,length(t2));  
x2=s2+v2; %������������ ������ � ��������� 1000  
akf_un=xcorr(x1,x1,10, 'unbiased'); %������
akf_bias=xcorr(x2,x2,100, 'unbiased');%������
figure(3) %������� ��������� ������
subplot (2, 1, 1); plot (akf_un); grid on;
title('������ ������ �������(��������� 100)');
ylabel('��������'); xlabel('³����');
subplot (2, 1, 2); plot (akf_bias); grid on;
title('������ ������ �������(��������� 1000)');
ylabel('��������'); xlabel('³����');

%=== �������� #2.1 ===
% ������������ ������� ��� ���� (eeg1-p4.dat) 
eeg_p4 =load('eeg1-p4.dat'); % ������ EE�
eeg1_p4= detrend(eeg_p4); %����������� �������
figure(4)
plot(eeg1_p4);
title ('������ ���')
ylabel('��������'); xlabel('³����');

%=== �������� #2.2 ===
% �������� ����� ������� ��� �� t1 = 4,7 � �� t2 = 5,8 �
fs=100;% xfcnjnf lbcrhtnbpfws]
t = (0:length(eeg1_p4)-1)/fs; % ������ ����
n1 = (4.7*fs) + 1; 
n2 = (5.8*fs) + 1;
figure(5)
plot(t(n1:n2),eeg1_p4(n1:n2));
title ('����� ������� �� t1 = 4,7 � �� t2 = 5,8 �')
ylabel('��������'); xlabel('³����');

%== �������� #2.3 ===
% ���������� �������� ������ ��� ������� ���
akf = xcorr(eeg1_p4,eeg1_p4, fix(0.9*length(eeg1_p4)), 'unbiased'); %���������� ������
figure(6)
plot(akf);
title('������ �� ������ ������ ����� ������� eeg1-p4.dat');
ylabel('��������'); xlabel('³����');

%=== �������� #2.4 ===
% ���������� ����������� �������� �������
Sx = abs(fft(eeg1_p4)/length(eeg1_p4)); %���������� ������ �������� ������������ �����
B=Sx(1:length(eeg1_p4)/2+1);%��������� ������������� ����������� �������� 
B(2:end-1)=2*B(2:end-1);
f =fs*(0:(length(eeg1_p4)/2))/length(eeg1_p4);
figure(7)
plot(f,B), grid on;
title('����������� ��������');
ylabel('Sxx'); xlabel('³����');

%=== �������� #2.5 ===
% ������������ ������� ��� ���� (eeg1-f3.dat)
eeg_f3 = load('eeg1-f3.dat'); %������������ �������
eeg1_f3 = detrend(eeg_f3); %����������� 
t = (0:length(eeg1_f3)-1)/fs; % ������ ����
figure(8)
plot(t,eeg1_f3)%�������� �������
title ('������ ��� eeg1-f3.dat')
ylabel('��������'); xlabel('³����');

% �������� ����� ������� ��� �� t1 = 4,2 � �� t2 = 4,96 �
t_f3 = (0:length(eeg1_f3)-1)/fs;
n1 = (4.2*fs) + 1; %���������� �����
n2 = (4.96*fs) + 1; %���������� �����
figure(8)
plot(t_f3(n1:n2),eeg1_f3(n1:n2));
title ('����� ������� �� t1 = 4,2 � �� t2 = 4,96 �')
ylabel('��������'); xlabel('³����');

% ���������� �������� ������ ��� ������� ���
akf2 = xcorr(eeg1_f3,eeg1_f3,fix(0.2*length(eeg1_f3)) , 'unbiased'); %���������� ������
figure(10)
plot(akf2);
title('������ �� ������ ������ ������� ��� �������');
ylabel('��������'); xlabel('³����');

% ���������� ����������� �������� �������
Sxx1 = abs(fft(eeg1_f3))/length(eeg1_f3);% ���������� ������ ����������� ��������
Bx1=Sxx1(1:length(eeg1_f3)/2+1);%��������� ������������� ����������� �������� 
Bx1(2:end-1)=2*Bx1(2:end-1);
f =fs*(0:(length(eeg1_f3)/2))/length(eeg1_f3);
figure(11)
plot(f,Bx1), grid on;
title('����������� ��������');
ylabel('Sxx1'); xlabel('³����');

% ������������ �������  ��� ���� (eeg1-o2.dat)
% ������������ ������� ��� ���� (eeg1-c3.dat)
eeg_o2 = load('eeg1-o2.dat');  % ������������ �������  ��� ���� (eeg1-o2.dat)
eeg1_o2 = detrend(eeg_o2);  %����������� 
eeg_c3 = load('eeg1-c3.dat'); % ������������ ������� ��� ���� (eeg1-c3.dat)
eeg1_c3 = detrend(eeg_c3); %����������� 
t_o2 = (0:length(eeg1_o2)-1)/fs; % ������ ����
t_c3 = (0:length(eeg1_c3)-1)/fs; % ������ ����
figure(12)
subplot(2,1,1); plot(t_o2,eeg1_o2);
title ('������ ��� eeg1-o2.dat')
ylabel('��������'); xlabel('³����');
subplot(2,1,2); plot(t_c3,eeg1_c3);
title ('������ ��� eeg1_c3')
ylabel('��������'); xlabel('³����');

%=== �������� #3.2 ===
% �������� ����� ������� ��� �� t1 = 5,71 � �� t2 = 6,78 �
t_o2 = (0:length(eeg1_o2)-1)/fs; % ������ ����
t_c3 = (0:length(eeg1_c3)-1)/fs; % ������ ����
n1 = (5.71*fs) + 1;   %���������� �����
n2 = (6.78*fs) + 1;   %���������� �����
figure(13)
subplot (2, 1, 1); plot(t_o2(n1:n2),eeg1_o2(n1:n2));
title ('����� ������� eeg1_o2 �� t1 = 4,2 � �� t2 = 4,96 �')
ylabel('��������'); xlabel('³����');
subplot (2, 1, 2); plot(t_c3(n1:n2),eeg1_c3(n1:n2));
title ('����� ������� eeg1_c3 �� t1 = 4,2 � �� t2 = 4,96 �')
ylabel('��������'); xlabel('³����');

%=== �������� #3.3 ===
% ���������� ��� ������� ���
vkf1 = xcorr(eeg1_o2,eeg1_c3,'biased'); %���������� ���
figure(14)
plot(vkf1);
title ('��� ������� ��� eeg1_o2 �� eeg1_c3')
ylabel('��������'); xlabel('³����');

%=== �������� #3.4 ===
% ���������� ������ ����������� �������� ������� ���
Sxy = abs(fft(vkf1)/length(vkf1));  % ���������� ������ ����������� ��������
Bz=Sxy(1:length(vkf1)/2+1);%��������� ������������� ����������� �������� 
Bz(2:end-1)=2*Bz(2:end-1);
f =fs*(0:(length(vkf1)/2))/length(vkf1);
figure(15)
plot(f,Bz), grid on;
title ('������ ����������� ��������� ������� ������� ���')
ylabel('Bz'); xlabel('�������');

%=== �������� #3.5 ===
% ������������ ������� ��� ���� (eeg1-p3.dat) �� ������� (eeg1-p4.dat)
eeg1_p3 = load('eeg1-p3.dat'); % ������������
eeg1_p3 = detrend(eeg1_p3);
eeg1_p4 = load('eeg1-p4.dat');  % ������������
eeg1_p4 = detrend(eeg1_p4);
t_p3 = (0:length(eeg1_p3)-1)/fs; % ������ ����
t_p4 = (0:length(eeg1_p4)-1)/fs; % ������ ����
figure(16)
subplot (2, 1, 1); plot(t_p3, eeg1_p3);
title ('������ ������� eeg1-p3.dat')
subplot (2, 1, 2); plot(t_p4, eeg1_p4);
title ('������ ������� eeg1-p4.dat')

% �������� ����� ������� ��� �� t1 = 4,7 � �� t2 = 5,8 �
t_p3 = (0:length(eeg1_p3)-1)/fs;
t_p4 = (0:length(eeg1_p4)-1)/fs;
n1 = (4.7*fs) + 1; 
n2 = (5.8*fs) + 1;
figure(17)
subplot (2, 1, 1); plot(t_p3(n1:n2),eeg1_p3(n1:n2));
title ('����� ������� eeg1_p3 �� 4,7� �� 5,8�')
ylabel('��������'); xlabel('³����');
subplot (2, 1, 2); plot(t_p4(n1:n2),eeg1_p4(n1:n2));
title ('����� ������� eeg1-p4 �� 4,7� �� 5,8�')
ylabel('��������'); xlabel('³����');
% ���������� ��� ������� ���
vkf2 = xcorr(eeg1_p3,eeg1_p4); 
figure(18)
plot(vkf2);
title ('% ���������� ��� ������� ��� eeg1_p3 �� eeg1_p4')
ylabel('��������'); xlabel('³����');

% ���������� ������ ����������� �������� ������� ���
Svkf2 = abs(fft(vkf2)/length(vkf2)); %���������� ������ �������� ������������ �����
Bn=Svkf2(1:length(vkf2)/2+1);%��������� ������������� ����������� �������� 
Bn(2:end-1)=2*Bn(2:end-1);
f =fs*(0:(length(vkf2)/2))/length(vkf2);
figure(19)
plot(f,Bn), grid on;
title ('������ ����������� ��������� ������� eeg1_p3 �� eeg1-p4')
ylabel('��������'); xlabel('�������');

%=== �������� #3.6 ===
% ������������ ������� ��� ���� (eeg1-f3.dat)
% ������������ ������� ��� ���� (eeg1-f4.dat)
eeg1_f3 = load('eeg1-f3.dat');
eeg1_f3 = detrend(eeg1_f3);
eeg1_f4 = load('eeg1-f4.dat');
eeg1_f4 = detrend(eeg1_f4);
t_f3 = (0:length(eeg1_f3)-1)/fs;
t_f4 = (0:length(eeg1_f4)-1)/fs;
figure(20)
subplot (2, 1, 1); plot(t_f3, eeg1_f3);
title ('������ ������� eeg1-f3.dat')
ylabel('��������'); xlabel('³����');
subplot (2, 1, 2); plot(t_f4, eeg1_f4);
title ('������ ������� eeg1-f4.dat')
ylabel('��������'); xlabel('³����');

% �������� ����� ������� ��� �� t1 = 4,13 � �� t2 = 4,96 � t3 = 1,43 � �� t4 = 2,26 �
n1=(4.13*fs) + 1;%���������� �����
n2=(4.96*fs) + 1;%���������� �����
n3=(1.43*fs) + 1;%���������� �����
n4=(2.26*fs) + 1;%���������� �����
figure(21)
subplot (2, 1, 1); plot(t_f3(n1:n2),eeg1_f3(n1:n2));
title ('����� ������� eeg1_f3 �� 4,13� �� 4.96�')
ylabel('��������'); xlabel('³����');
subplot (2, 1, 2); plot(t_f4(n3:n4),eeg1_f4(n3:n4));
title ('����� ������� eeg1_f4 �� 1.43� �� 2.26�')
ylabel('��������'); xlabel('³����');

% ���������� ��� ������� ���
vkf3 = xcorr(eeg1_f3, eeg1_f4, 250); % ���������� ���
figure(22)
plot(vkf3);
title ('��� ���� ������� eeg1_f3 �� eeg1_f4')
ylabel('��������'); xlabel('³����');

% ���������� ������ ����������� �������� ������� ���
Sx = abs(fft(vkf3)/length(vkf3)); %���������� ������ �������� ������������ �����
Bq=Sx(1:length(vkf3)/2+1);%��������� ������������� ����������� �������� 
Bq(2:end-1)=2*Bq(2:end-1);
f =fs*(0:(length(vkf3)/2))/length(vkf3);
figure(23)
plot(f,Bq), grid on;
title ('������ ����������� ��������� ���� ������� eeg1_p3 �� eeg1-p4')
ylabel('��������'); xlabel('�������');



