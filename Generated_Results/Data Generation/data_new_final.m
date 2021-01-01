clc;
snrdb =0:1:54;
snr = 10.^(snrdb/10);

% basic parameters initialization 
lambda = 100;
r = 45e-9;
d = 500e-9;
D = 4.265e-10;
delta_T = 9e-6;
T = 30*delta_T;
L = 100;
global sum_Cj;
summation1=0;
x =0:1:65;
sum_Cj=0;


Co = zeros(1,length(snrdb));
s=0;
temp=0;

% Generating the data
data1 = rand(1,L)>0.5;
%display(data1);
P_0 = (r/d)*(erfc((d-r)/sqrt(4*D*T)));
%Ntx = 2.lambda.*T.(10.^(snrdb./10))./P_0;
%Ntx = 2.*lambda.*T.*(10.^(snr./10))./P_0;
Ntx = 2.*lambda.*T.*(snr)./P_0;

P_i1(1)=P_0;
display(P_i1(1));
final_term = 100000;
ri = zeros(55,100);
display(length(snrdb));

for ii=1:1:length(snrdb)
    for k = 2:1:L
        P_i1(k) = (r/d)*(erfc((d-r)/sqrt(4*D*k*T))-erfc((d-r)/sqrt(4*D*(k-1)*T))); % generating P_(i-1)
    end
    
    for j = 1:1:L
        check1(j) = Ntx(ii)*P_i1(j);
    end
    
    crack=data1.*check1;
    display(size(crack));
    summation1 = sum(crack);
    
    Co(ii) = snr(ii)*2*lambda*T;
    average1 = (lambda*T + summation1);
    average2 = average1 + data1.*Co(ii);
    
    ri(ii,:) = poissrnd(average2); % applying possion random varible...
end


disp('here is ri');
%disp(ri);
disp(size(ri))
ri=reshape(ri', [100*55,1]);
%disp(ri);
disp('size of new ri');
disp(size(ri));
data2=repmat(data1',55,1);
disp(size(data2));
data2=reshape(data2,[100*55,1]);
%display(data2);

% string in the csv file...
writematrix(data2,'F:/acadamic/sem-5/wirelss networking lab/project/data/data_test_5500.csv');
writematrix(ri,'F:/acadamic/sem-5/wirelss networking lab/project/data/ri_test_5500.csv');
