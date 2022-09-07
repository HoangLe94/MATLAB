clear all; clc;

%Planar antenna on xy-plane
%Good

f = 25; %frequency in GHz
%t = 0:0.25:1; %time sampling in ns
%A = 1 ; %amplitude in m

c = 3e8 ; %speed of light in m/s

Pos_x = 0:1:5; %position of elements of antenna array
Pos_y = 0:1:7; %position of elements of antenna array

omega = 2*pi*f;
lambda = c/(f*10^9); %wavelength in m

Delta = lambda/2; %distance between 2 elements in m

s = [10];

theta_DOA = [30];
phi_DOA = [45];

for i = 1:length(theta_DOA)
    mu_DOA = -2*pi/lambda*Delta*cos(deg2rad(theta_DOA(i)))*sin(deg2rad(phi_DOA(i)))*Pos_x';
    a_mu_DOA(:,i) = transpose(exp(1j*mu_DOA));
    
    nu_DOA = -2*pi/lambda*Delta*sin(deg2rad(theta_DOA(i)))*sin(deg2rad(phi_DOA(i)))*Pos_y';
    a_nu_DOA(:,i) = transpose(exp(1j*nu_DOA)); 
end

SNR = db2pow(inf);

%x = awgn(a_DOA*transpose(s), SNR, 'measured');
x = reshape(awgn(s*a_mu_DOA*transpose(a_nu_DOA), SNR, 'measured'),[],1);

%w = reshape(mu_DOA*nu_DOA',[],1);
%y = reshape(a_mu_DOA*a_nu_DOA',[],1);

%z = mu_DOA'*Y*nu_DOA;

R = x*x';

%P = w'*R*w/(length(Pos_x)*length(Pos_y));


for i = 1:101
    theta(i) = -pi/2 + (i-1)*pi/100;
    theta_plot(i) = rad2deg(theta(i));
    
    
    
    for j = 1:101
        phi(j) = 0 + (j-1)*pi/100;
        phi_plot(j) = rad2deg(phi(j));
        
        mu(:,i,j) = -2*pi/lambda*Delta*cos(theta(i))*sin(phi(j))*Pos_x';
        a_mu(:,i,j) = transpose(exp(1j*mu(:,i,j)));
        
        nu(:,i,j) = -2*pi/lambda*Delta*sin(theta(i))*sin(phi(j))*Pos_y';
        a_nu(:,i,j) = transpose(exp(1j*nu(:,i,j)));
    end
end


for i = 1:101
    for j = 1:101
        w(:,i,j) = reshape(mu(:,i,j)*nu(:,i,j)',[],1);
        
        a(:,i,j) = reshape(a_mu(:,i,j)*transpose(a_nu(:,i,j)),[],1);
        
        %R = y(:,i,j)*y(:,i,j)';
        
        P(i,j) = w(:,i,j)'*R*w(:,i,j)/(length(Pos_x)*length(Pos_y));
        
        %P(i,j) = a(:,i,j)'*R*a(:,i,j)/(a(:,i,j)'*a(:,i,j));
        %P(i,j) = 1/(a(:,i,j)'*inv(R)*a(:,i,j));
    end
end

surf(theta_plot,phi_plot,abs(P));

%hold on
xlabel('Theta')
ylabel('Phi')
zlabel('P')


%[V,D] = eig(R);

%q = diag(D);

%[B,I] = mink(q, length(Pos) - length(theta_DOA));

%V(:,length(Pos_x)-length(theta_DOA)+1:length(Pos_x)) = [];
%V = fliplr(V);

%{
for i = 1:1000
    theta(i) = -pi/2 + (i-1)*pi/1000;
    theta_plot(i) = rad2deg(theta(i));
    mu(i) = -2*pi/lambda*Delta*sin(theta(i));
    
    a(:,i) = transpose(exp(1j*mu(i)*Pos_x));
    %P(i) = trace(a(:,i)*inv(a(:,i)'*a(:,i))*a(:,i)'*R);
    P(i) = 1/(a(:,i)'*V*V'*a(:,i));
end

[peaks ind] = findpeaks(abs(P));

[top_peaks ind_top] = maxk(peaks, length(theta_DOA));

DOA = theta_plot(ind(ind_top))


figure(1); clf;
f1 = plot(theta_plot, P);
%}