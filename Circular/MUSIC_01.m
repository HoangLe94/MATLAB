clear all; clc;

%Analyze automatically parameters of arrival
%Good

f = 25; %frequency in GHz
%t = 0:0.25:1; %time sampling in ns
%A = 1 ; %amplitude in m

c = 3e8 ; %speed of light in m/s

Pos = 0:1:10; %position of elements of antenna array

omega = 2*pi*f;
lambda = c/(f*10^9); %wavelength in m

Delta = lambda/2; %distance between 2 elements in m

s = [10];

theta_DOA = [30];
for i = 1:length(theta_DOA)
    mu_DOA(i) = 2*pi/lambda*Delta*sin(deg2rad(theta_DOA(i)));
    a_DOA(:,i) = transpose(exp(-1j*mu_DOA(i)*Pos));  
end

SNR = db2pow(inf+);

x = awgn(a_DOA*transpose(s), SNR, 'measured');

R = x*x';

[V,D] = eig(R);

q = diag(D);

%[B,I] = mink(q, length(Pos) - length(theta_DOA));

V(:,length(Pos)-length(theta_DOA)+1:length(Pos)) = [];
V = fliplr(V);

for i = 1:1000
    theta(i) = -pi/2 + (i-1)*pi/1000;
    theta_plot(i) = rad2deg(theta(i));
    mu(i) = -2*pi/lambda*Delta*sin(theta(i));
    
    a(:,i) = transpose(exp(1j*mu(i)*Pos));
    %P(i) = trace(a(:,i)*inv(a(:,i)'*a(:,i))*a(:,i)'*R);
    P(i) = 1/(a(:,i)'*V*V'*a(:,i));
end

[peaks ind] = findpeaks(abs(P));

[top_peaks ind_top] = maxk(peaks, length(theta_DOA));

DOA = theta_plot(ind(ind_top))


figure(1); clf;
f1 = plot(theta_plot, P);