%% Construct data matrix
clc
clear
close all

s1=exp(j*pi/4);
s2=exp(j*pi*3/4);
s3=exp(j*pi*5/4);
s4=exp(j*pi*7/4);
Collection=[s1,s2,s3,s4];
N=500;
P=8;
numbers = randi([1, 4], 1, N);

for i=1:length(numbers)
    s(i)=Collection(numbers(i)); % Construct s
end
sigma=0;
x=gendata_conv2(s,P,N,sigma);
for i=1:4*P
    for m=1:N-3
        X(i,m)=x((m-1)*P+i);
    end
end
rank(X)
%% Zero-forcing and Wiener 

s1=exp(j*pi/4);
s2=exp(j*pi*3/4);
s3=exp(j*pi*5/4);
s4=exp(j*pi*7/4);
Collection=[s1,s2,s3,s4];
N=500;
P=4;
L=1;
Ns=N-L+1;
numbers = randi([1, 4], 1, N);

for i=1:length(numbers)
    s(i)=Collection(numbers(i)); % Construct s
end
sigma=0.5;

x=gendata_conv2(s,P,N,sigma);
for i=1:2*P
    for m=1:Ns-1
        X_2(i,m)=x((m-1)*P+i);
    end
end

%h=[1;1;-1;-1;1;1;-1;-1];
h=[1;-1;1;-1];
H=[zeros(1,P),h';h',zeros(1,P)]';

w_zf=H*pinv(H'*H);
s_rec_zf=w_zf'*X_2;


w_wiener=pinv(H*H'+sigma^2*eye(size(H*H',1)))*H;
s_rec_wiener=w_wiener'*X_2;

figure
subplot(2,2,1)
plot(real(s_rec_zf(1,:)), imag(s_rec_zf(1,:)),'.')
title('Estimated symbols(1) by ZF equalizer')
subplot(2,2,2)
plot(real(s_rec_zf(2,:)), imag(s_rec_zf(2,:)),'.')
title('Estimated symbols(2) by ZF equalizer')
subplot(2,2,3)
plot(real(s_rec_wiener(1,:)), imag(s_rec_wiener(1,:)),'.')
title('Estimated symbols(1) by Wiener equalizer')
subplot(2,2,4)
plot(real(s_rec_wiener(2,:)), imag(s_rec_wiener(2,:)),'.')
title('Estimated symbols(2) by Wiener equalizer')