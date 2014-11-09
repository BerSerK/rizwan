
clear all;
clc
format long
close all

global A Pr Nb Nt Le M S L

S  = 1;
A = 0.5;
M = 1; % stagnation
Pr = 0.7;
Nb = 0.3;
Nt = 0.7;
Le = 1.0;

L_range = [-3:0.2:1];
N_L = numel(L_range);

y_prime0_n = zeros(N_L, 1);
for j = 1:N_L
    L = L_range(j);
    [sol] = SHR;
    y_prime0 = sol.y;
    y_prime0_n(j) = y_prime0(3,1);
    display(L)
    
end
plot(L_range, y_prime0_n);
% hold on
% 
% 
% S = 1;
% y_prime0_n = zeros(N_L, 1);
% for j = 1:N_L
%     L = L_range(j);
%     [sol] = SHR;
%     y_prime0 = sol.y;
%     y_prime0_n(j) = y_prime0(3,1);
%     display(L)
%    
% end
% 
% plot(L_range, y_prime0_n);
% hold on
% 
% S  = 2;
% y_prime0_n = zeros(N_L, 1);
% for j = 1:N_L
%     L = L_range(j);
%     [sol] = SHR;
%     y_prime0 = sol.y;
%     y_prime0_n(j) = y_prime0(3,1);
%     display(L)
%     
% end
% 
% plot(L_range, y_prime0_n);
% 
% 
