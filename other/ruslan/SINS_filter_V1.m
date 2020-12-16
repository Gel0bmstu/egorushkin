%% загрузка данных
clear
cd 'c:\Users\Руслан\Desktop\Оценка точности'
load('dannie_V1')
%% вычисление оценки 
Pk=[ 1000 0 0;0 10^(-5) 0;0 0 10^(-16)];
GQG=[1e-12 0 0; 0 1e-12 0; 0 0 1e-12];
Rn=[10^-2 0;0 10^-16];
Xoc(1,1)= [ 0 ];
Xoc(2,1)= [ 0 ];
Xoc(3,1)= [ 0 ];
TAR=ones(3);
reTARd=0;
oc1=0;
oc2=0;
oc3=0;
SAF=0;
for k=2:5000
    Pkk1=F*Pk*F'+GQG;
    Kk=(Pkk1*H')/(H*Pkk1*H'+Rn);
    Xoc=[Xoc (F*Xoc(:,(k-1))+Kk*(Z(:,k)-H*F*Xoc(:,(k-1))))];
    Pk=(eye(3)-Kk*H)*Pkk1;
    TAR=TAR*((eye(3,3)-Kk*H)*F);
    oc1=[oc1 TAR(1,1)];
    oc2=[oc2 TAR(2,2)];
    oc3=[oc3 TAR(3,3)];
    SAF = SAF + (X(:, k) - Xoc(:, k)).^2;
end
SAF=SAF/5000;

%% проверка
figure('Name','Ошибка по скорости_V1','NumberTitle','off');
plot(t,Xoc(1,:),t,X(1,:))
figure('Name','Ошибка по углу_V1','NumberTitle','off');
plot(t,Xoc(2,:),t,X(2,:))
figure('Name','Дрейф_V1','NumberTitle','off');
plot(t,Xoc(3,:),t,X(3,:))
figure
plot(t,oc1)
figure
plot(t,oc2)
figure
plot(t,oc3)

