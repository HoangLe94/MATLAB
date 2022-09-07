clear all; clc;

%Planar antenna on zx-plane for several incident waves
%Good

f = 25; %frequency in GHz
%t = 0:0.25:1; %time sampling in ns
%A = 1 ; %amplitude in m

c = 3e8 ; %speed of light in m/s

Pos_x = 0:1:5; %position of elements of antenna array
Pos_z = 0:1:7; %position of elements of antenna array

omega = 2*pi*f;
lambda = c/(f*10^9); %wavelength in m

Delta = lambda/2; %distance between 2 elements in m

s = [10];
SNR = db2pow(-20);

theta_DOA = [-10];
phi_DOA = [58];

%{
for i = 1:length(Pos_x)
    for j = 1:length(Pos_z)
        A(i,j) = exp(-2*pi*1j/lambda*Delta*(Pos_x(i)*cos(deg2rad(theta_DOA)) + Pos_z(j)*sin(deg2rad(theta_DOA))*cos(deg2rad(phi_DOA))));
    end
end
%}

for i = 1:length(theta_DOA)
    mu_DOA(:,i) = -2*pi/lambda*Delta*sin(deg2rad(theta_DOA(i)))*Pos_z';
    a_mu_DOA(:,i) = exp(1j*mu_DOA(:,i));
    
    nu_DOA(:,i) = -2*pi/lambda*Delta*cos(deg2rad(theta_DOA(i)))*sin(deg2rad(phi_DOA(i)))*Pos_x';
    a_nu_DOA(:,i) = exp(1j*nu_DOA(:,i)); 
end

A = reshape(a_mu_DOA*transpose(a_nu_DOA),[],1);

x = awgn(A*s, SNR, 'measured');

R = x*x';

[V,D] = eig(R);

q = diag(D);

%[B,I] = mink(q, length(Pos) - length(theta_DOA));

V(:,length(x)-length(theta_DOA)+1:length(x)) = [];
V = fliplr(V);


for i = 1:181
    theta(i) = -pi/2 + (i-1)*pi/180;
    theta_plot(i) = rad2deg(theta(i));
     
    for j = 1:181
        phi(j) = -pi/2 + (j-1)*pi/180;
        phi_plot(j) = rad2deg(phi(j));
        
        mu(:,i,j) = -2*pi/lambda*Delta*sin(theta(i))*Pos_z';
        a_mu(:,i,j) = exp(1j*mu(:,i,j));
        
        nu(:,i,j) = -2*pi/lambda*Delta*cos(theta(i))*sin(phi(j))*Pos_x';
        a_nu(:,i,j) = exp(1j*nu(:,i,j));
    end
end

for i = 1:181
    for j = 1:181
        %w(:,i,j) = reshape(mu(:,i,j)*transpose(nu(:,i,j)),[],1);
        
        a(:,i,j) = reshape(a_mu(:,i,j)*transpose(a_nu(:,i,j)),[],1);
        
        %R = y(:,i,j)*y(:,i,j)';
        
        %P(i,j) = w(:,i,j)'*R*w(:,i,j)/(length(Pos_x)*length(Pos_z));
        
        %P(i,j) = a(:,i,j)'*R*a(:,i,j)/(a(:,i,j)'*a(:,i,j));
        %P(i,j) = 1/(a(:,i,j)'*inv(R)*a(:,i,j));
        
        P(i,j) = 1/(a(:,i,j)'*V*V'*a(:,i,j));
    end
end


surf(phi_plot,theta_plot,abs(P));

%hold on
xlabel('Phi')
ylabel('Theta')
zlabel('P')

max_index = islocalmax(abs(P));

local_max_P = [];
for i = 1:181
    for j = 1:181
        if max_index(i,j) == 1
            local_max_P = [local_max_P; [abs(P(i,j)) i j]];
        end
    end
end

[max_P index] = maxk(local_max_P(:,1), length(theta_DOA));

theta_estimate = theta_plot(local_max_P(index,2))
phi_estimate = theta_plot(local_max_P(index,3))
