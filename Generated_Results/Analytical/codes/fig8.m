clc;
close all;
clear all;

%defining the constant values
lambda = 100;
r = 45e-9;
d = 500e-9;
D = 4.265e-10;
delta_T = 9e-6;
T = 30*delta_T;
L = 5;
global sum_Cj;
sum=0;
sum1=0;
x = 0:1:60;
sum_Cj=0;

data = rand(1,L)>0.5;

%for calculation of sum values
for ii=1:1:5
    sum = sum + Cj_fun(ii);
    sum1 = sum1 + data(ii)*Cj_fun(ii);
end

avg = sum + lambda*T;
avg1 = avg + 20;
final_term = 100000;

for jj=1:1:length(x)
    %probability distribution for si=0
    factor = avg.^(x(jj));
    fact = factorial(x(jj));
    final(jj) = (exp(-avg)*(factor))/fact;
    
    %probability distribution for si=1
    factor1 = avg1.^(x(jj));
    final1(jj) = (exp(-avg1)*(factor1))/fact;
    
    lam1 = lambda*T + sum1;
    lam2 = lam1 + 20;
    first_term = Q_fun(lam1, x(jj));
    last_term = 1 - Q_fun(lam2, x(jj));
    final_term_next = 0.5*(first_term + last_term);
    
    %for finding the optimal threshold
    if(final_term_next < final_term)
        final_term = final_term_next;
        final_tao = x(jj);
    end
end

%threshold for the equibprobability condition
t3 = 20/avg;
t2 = log(1+t3); 
t1 = 20/t2;

%plotting the graphs
figure
stem(t1, 0.08, 'Marker','none');
hold on
stem(final_tao, 0.08, 'Marker', 'none');
semilogy(x,final,'r.-');
semilogy(x,final1,'k*-');
grid on                               %displaying the grid lines
ax = gca;                             % setting the font size
ax.FontSize = 11;                     %font-size
xlabel('Received Particles(r_i)')     % x-axis label name
ylabel('Probability')                 % y-axis label name
legend({'Threshold derived from equiprob','The optimal threshold','Aprx distr(si=0)', 'Aprx distr(si=1)' },'Location','northeast')   %for individual plot
title('Approximated distributions')   %title of the graph
%axis([0 60 0.001 0]);                %setting the x-axis and y-axis values
hold off

savefig('figure8.fig');

%defining the functions to be used
function sum_Cj = Cj_fun(j)
    Ntx = 10^2;
    global sum_Cj;
    sum_Cj = sum_Cj + Ntx*prob_j(j);
end

%function for defining the hitting probability
function y = prob_j(j)
    lambda = 100;
    r = 45e-9;
    d = 500e-9;
    D = 4.265e-10;
    delta_T = 9e-6;
    T = 30*delta_T;
    L = 5;
    
   x = r/d;
   t1 = ((d-r)/(sqrt(4*D*(j+1)*T)));
   t2 = ((d-r)/(sqrt(4*D*(j)*T)));
   y = x*(erfc(t1) - erfc(t2));  
end

%function for defining the Q-function
function ans = Q_fun(lambda, n)
    u1 = factorial(n);
    u3 = lambda.^n
    u2 = exp(-lambda)*u3;
    ans = u2/u1;
end