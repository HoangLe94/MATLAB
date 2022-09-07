clear all; clc;

%{
t = 0:pi/50:10*pi;
st = sin(t);
ct = cos(t);
plot3(st,ct,t);
%}

A = 10*ones(10,10);
A(2:4,2:4) = 22;
A(3,3) = 21;
A(6:8,6:8) = 33; 
A(2,7) = 44;
A(3,8) = 45;
A(4,9) = 44;