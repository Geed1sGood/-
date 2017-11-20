fs = 20;
f0 = 5 ;
A = 0.5;
D = 0.25; 
n = 512; 
p=3.1416
t = [0:0.02:5];		% вектор часу
h=t*2*p*f0;
s =sin(h);
v = D*randn(1,length(t));  
x=s+v;  

t1=-5:0.02:5
maxlag = fix(0.1:0.3);
r = xcorr(x, maxlag)
[r,lags]=xcorr(x);
figure (2)
plot(t1,r)