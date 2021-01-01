clear all;
clc;

Tau=[20:2:200];
SNR=30;
lambda0=100;
deltaT=0.000009;
slot_len=30*deltaT;
channel_len=[1:5];
R=45*(10^(-9));
d=500*(10^(-9));
D=4.625*(10^(-10));
L=5;

% for 30 delta
for i=1:length(channel_len)
    P(i)=(R/d)*(erfc((d-R)/sqrt(4*D*(i+1)*slot_len))-erfc((d-R)/sqrt(4*D*(i)*slot_len)));
    NTx(i)=(2*(lambda0)*(slot_len)*10^(SNR/10))/P(1);
end

for i=1:length(channel_len)
    C(i)=NTx(i)*P(i);
    sj=randi([0,1],1,i);
    for j=1:length(Tau)
        Petau1(j)=(1/2)*[gammainc(lambda0*slot_len+sum(sj.*C),(Tau(j))) + 1 - gammainc(lambda0*slot_len+sum(sj.*C)+C(1),(Tau(j)))];  
    end
end

% for 50 delta
slot_len=50*deltaT;
for i=1:length(channel_len)
    P1(i)=(R/d)*(erfc((d-R)/sqrt(4*D*(i+1)*slot_len))-erfc((d-R)/sqrt(4*D*(i)*slot_len)));
    NTx1(i)=(2*(lambda0)*(slot_len)*10^(SNR/10))/P1(1);
end

for i=1:length(channel_len)
    C1(i)=NTx1(i)*P1(i);
    sj1=randi([0,1],1,i);
    for j=1:length(Tau)
        Petau2(j)=(1/2).*[gammainc(lambda0.*slot_len+sum(sj1.*C1),Tau(j)) + 1 - gammainc(lambda0.*slot_len+sum(sj1.*C1)+C1(1),Tau(j))];  
    end
end

%ploting statements
close all;
figure
hold on
semilogy(Tau,Petau1,'k*-','LineWidth',1);
hold on 
semilogy(Tau,Petau2,'r*-','LineWidth',1);

grid on                               %displaying the grid lines
ax = gca;                             % setting the font size
ax.FontSize = 11;                     %font-size
xlabel('Threshold Value')                     % x-axis label name
ylabel('BER(Bit Error Rate)')         % y-axis label name
legend({'Slot Length 30 deltaT','Slot Length 50 deltaT'},'Location','northeast')   %for individual plot
title('BER as a function of Tau')   %title of the graph
axis([20 200 0 1]);                %setting the x-axis and y-axis values
hold off

savefig('fig1.fig');