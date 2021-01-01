close all;
clear all;
clc; 

%range of SNR
snrdb =-5:1:55;
snr = 10.^(snrdb/10);

%defining the constants
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
x =0:1:65;
sum_Cj=0;
Co = zeros(1,length(snrdb));
s=0;
display(snr);

data = rand(1,L)>0.5;
display(data);

for ii=1:1:5
    sum = sum + Cj_fun(ii);
    display(sum);
    sum1 = sum1 + data(ii)*Cj_fun(ii);
    display(sum1);
end
avg = sum + lambda*T;
display(avg);

avgg = sum1 + lambda*T;
display(avgg);

avg1 = avg + 20;
final_term = 100000;

%loop for calculating the BER values
for ii=1:1:length(snrdb)
    Co(ii) = snr(ii)*2*lambda*T;
    a = sec_term(Co(ii), avg1);
    %calculating the sub-optimal BER
    Pe(ii)=0.5*a;
  
    %calculating the 
    for jj=1:1:length(x)     
        lam1 = lambda*T + sum1;
        lam2 = lam1 + Co(ii);
        first_term = 1-Q_fun(lam2, x(jj));
        last_term = Q_fun(lam1, x(jj));
        %first_term = (exp(lam1)*((lam1)^(jj)))/factorial(jj);
        %last_term = 1 - (exp(lam2)*((lam2)^(jj)))/factorial(jj);
        final_term_next = 0.5*(first_term + last_term);
       % s=s+final_term_next;
        if(final_term_next < final_term)
            final_term = final_term_next;
            final_tao = x(jj);
        end
    end
    
    b = sec_term1(Co(ii), avg1, final_tao);
    %Pe1(ii) = 0.5 - 0.5*(b);
    Pe1(ii)=0.5*b;
   
end

figure;
plot(snrdb, Pe, 'r.-');
hold on
plot(snrdb, Pe1, 'b*-');
axis([0 55 0 0.5]);
grid on                               %displaying the grid lines
ax = gca;                             % setting the font size
ax.FontSize = 11;                     %font-size
xlabel('SNR(dB)')                     % x-axis label name
ylabel('BER(Bit Error Rate')          % y-axis label name
legend({ 'Equiprobability zero-bit detector' ,'Optimal zero-bit detector'},'Location','northeast')   %for individual plot
title('Zero-bit reciever')            %title of the graph
%axis([20 120 0.5 0]);                %setting the x-axis and y-axis values
hold off

savefig('figure7.fig');

%defining the functions
%function for sub-optimal approach
function ans1 = sec_term(Co1, avg1)
    t3 = Co1/avg1;
    t2 = log(1+t3); 
    t1 = Co1/t2;  % calculating the value of tao
    
    ans1 = 0;
    lg0 = avg1;
    lg1 = avg1 + Co1;
    
    %calculating Pe as a function of tao
    for jj=1:1:t1
        d1 = lg0^(jj) * exp(-lg0);
        d2 = lg1^(jj) * exp(-lg1);
        ans1 = ans1 + (-d1+d2)/factorial(jj);
        
    end
    ans1 = 1 + ans1;
end

%function for optimal approach
function ans2 = sec_term1(Co1, avg1, final_tao)
    ans2 = 0;
    lg0 = avg1;
    lg1 = avg1 + Co1;
    
    %calculating Pe as a function of tao
    for jj=1:1:final_tao
        d1 = lg0^(jj) * exp(-lg0);
        d2 = lg1^(jj) * exp(-lg1);
        ans2 = ans2 + (-d1+d2)/factorial(jj);
    end
    ans2 = 1 + ans2;
end


function sum_Cj = Cj_fun(j)
    Ntx = 10^2;
    global sum_Cj;
    sum_Cj = sum_Cj + Ntx*prob_j(j);
end

%defining the hitting probability function
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

%defining the Q function 
function ans = Q_fun(lambda, n)
    ans=0;
    for j=1:1:n
         u1 = factorial(n);
         u3 = lambda.^n;
         u2 = exp(-lambda)*u3;
         ans = ans + u2/u1;
    end
   
end