clear all;
T=1;
A=10^(-4);
g=9.8;
R=6371000;
betta=5*10^(-8);
W(1)=[randi([-100,100])/100];
F=[ 1 -g*T 0 ; T/R 1 T ; 0 0 (1-betta*T) ];
X=[ 0; 0; 0];
G=[0 ; 0 ; 1];
H=[ 1  0  0];
Z(1)=0;
Zg(1)=0;
t=1:4000;
e=-1+2*rand;
Vk=round(e*100)/10000;
R=10^(-2);
Q(1,1)=[10^-13];
Xkd=[ 0; 0; 0 ];

Pk=[ 1000 0 0 ; 0 10^(-5) 0 ; 0 0 10^(-12) ];

% W = dlmread('WFile.txt');
% X = dlmread('X1File.txt');
% Z = dlmread('Z1File.txt');
for k=2:5000
    W=[W randi([-100,100])/100];
    X=[X (F*X(:,k-1)+G*W(k-1))];
    Z=[Z H*X(:,k)+((randi(201)-101)/100)*0.01];
   %Zg=[Zg H*X(:,k)+((randi(201)-101)/100)*0.000000001];
end

dlmwrite('WFile.txt',W,'delimiter','\t','precision',6);
dlmwrite('X1File.txt',X,'delimiter','\t','precision',6);
dlmwrite('Z1File.txt',Z,'delimiter','\t','precision',6);

SAF = 0;

for i=2:5000
    Pk_k_1 = F*Pk*F'+G*Q*G';
    Kk=(Pk_k_1*H')/(H*Pk_k_1*H'+R);
    Xkd=[Xkd (F*Xkd(:,(i-1))+Kk*(Z(i)-H*F*Xkd(:,(i-1))))];
    Pk=(eye(3)-Kk*H)*Pk_k_1;
    SAF = SAF + (X(3, i) - Xkd(3, i)).^2;
    
end
% figure('Name','Дрейф','NumberTitle','on');
% plot(t,Xkd(3,:),'y-',t,X(3,:),'m-')
% legend(['Xkd'],['X'],'FontSize',8);
% grid on;

%Обратный фильтр Калмана

XkdO(:, 5000) = [Xkd(1, 5000), Xkd(2, 5000), Xkd(3, 5000)];
PkO=Pk;
KkO=Kk;
for i = 4999:-1:1
    Pk_k_1O = F*PkO*F'+G*Q*G';
    KkO=(Pk_k_1O*H')/(H*Pk_k_1O*H'+R);
    XkdO(:, i) = (F*XkdO(:,(i+1))+KkO*(Z(i)-H*F*XkdO(:,(i+1))));
    PkO=(eye(3)-KkO*H)*Pk_k_1O;
    
end
reversed_XkdO = zeros(3, 5000);
reversed_XkdO(1, :) = flip(XkdO(1, :));
reversed_XkdO(2, :) = flip(XkdO(2, :));
reversed_XkdO(3, :) = flip(XkdO(3, :));

% SAF Search
SAF = SAF / 5000;
saf_min = 1e+10;
sdvig = floor(5000/3);
index = 0;
for i=1:sdvig
    saf_sum = 0;
    
    for j=1:i
        saf_sum = saf_sum + (Xkd(3, 5000 - j) - XkdO(3, 0 + j))^2;
    end
    saf = saf_sum / i;
    
    if saf < saf_min
        saf_min = saf;
        index = i;
    end
end
saf
saf_min
index

ideal_Xkd = zeros(1, 4500);
ideal_XkdO = zeros(1, 4500);
ideal = zeros(1, 4500);

for i = index+1:4500
    % ideal_kalm(i) = (XkdO(3, i+index) - Xkd(3, (i+index))) / 2;
    ideal_Xkd(i - index) = Xkd(3, i);
end

for i = 4500:-1:index+1
    % ideal_kalm(i) = (XkdO(3, i+index) - Xkd(3, (i+index))) / 2;
    ideal_XkdO(i) = XkdO(3, i-index);
end

for i=1:4500 
    ideal(i)=(ideal_Xkd(i)+ideal_XkdO(i))/2;
end
% ideal= zeros(1, 4500);
% ideal=(ideal_Xkd(3,:)-ideal_XkdO(3,:))/2;
% figure('Name','Ошибка по скорости','NumberTitle','off');
% plot(t,Z(1,t),'k-',t,Xkd(1,t),'y-',t,X(1,t),'m-',t,XkdO(1,t),'b-')
% grid on;
% legend(['Z'],['Xkd'],['X'],['XkdO'],'FontSize',8);
% figure('Name','Ошибка по углу','NumberTitle','off');
% plot(t,Xkd(2,t),'y-',t,X(2,t),'m-',t,reversed_XkdO(2,t),'b-')
% % % % % % % % grid on;
% legend(['Xkd'],['X'],['XkdO'],'FontSize',8);
figure('Name','Дрейф','NumberTitle','on');
plot(t,Xkd(3,t),'y-',t,X(3,t),'m-',t,XkdO(3, t),'b-',t,ideal(t),'g-');
legend(['Xkd'],['X'],['XkdO'],['IDEAl'], 'FontSize',8);
grid on;
figure();
plot(t,ideal_Xkd(t),'y-',t,X(3,t),'m-',t,ideal_XkdO(t),'b-',t,ideal(t),'g-');
legend(['Xkd'],['X'],['XkdO'],['IDEAl'], 'FontSize',8);
grid on;

index=66

ideal_Xkd2 = zeros(1, 4500);
ideal_XkdO2 = zeros(1, 4500);
ideal2 = zeros(1, 4500);

for i = index+1:4500
    ideal_kalm(i) = (XkdO(3, i+index) - Xkd(3, (i+index))) / 2;
    ideal_Xkd2(i - index) = Xkd(3, i);
end

for i = 4500:-1:index+1
    ideal_kalm(i) = (XkdO(3, i+index) - Xkd(3, (i+index))) / 2;
    ideal_XkdO2(i) = XkdO(3, i-index);
end

for i=1:4500 
    ideal2(i)=(ideal_Xkd2(i)+ideal_XkdO2(i))/2;
end

figure('Name','Дрейф','NumberTitle','on');
plot(t,X(3,t),'m-',t,ideal(t),'g-',t,ideal2(t),'r-');
legend(['X'],['IDEAl'],['IDEAl2'], 'FontSize',8);
grid on;
