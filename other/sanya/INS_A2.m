clear
T=1;
A=10^(-4);
g=9.8;
R=6371000;
betta=5*10^(-8);
W(1)=[randi([-100,100])/100];
F=[ 1 -g*T 0 ; T/R 1 T ; 0 0 (1-betta*T) ];
X=[ 0; 0; 0];
G=[ 0; 0; 0; 0; 0; 0;0 ; 0 ; 1];
H=[ 1  0  0];
Z(1)=0;
Zg(1)=0;
t=(1:5000);
e=-1+2*rand;
Vk=round(e*100)/10000;
R=10^(-2);
Q(1,1)=[10^-13];
XA2=[ 0; 0; 0];
Pk=[ 1000 0 0 ; 0 10^(-5) 0 ; 0 0 10^(-12) ];

W = dlmread('WFile.txt');
X = dlmread('X1File.txt');
Z = dlmread('Z1File.txt');
% for k=2:5000
% %     W=[W randi([-100,100])/100];
%     X=[X (F*X(:,k-1)+G*W(k-1))];
%     Z=[Z H*X(:,k)+((randi(201)-101)/100)*0.01];
%    %Zg=[Zg H*X(:,k)+((randi(201)-101)/100)*0.000000001];
% end
% 
% dlmwrite('WFile.txt',W,'delimiter','\t','precision',6);
% dlmwrite('X1File.txt',X,'delimiter','\t','precision',6);
% dlmwrite('Z1File.txt',Z,'delimiter','\t','precision',6);

% Адаптивный фильтр второго рода
C_k2=[0];
NU_k2=[0];
Q_k=[0];
SAF_A2=[0];
Kk=0;
for i=2:5000
    NU_k2=[NU_k2 (Z(1,i)-H*XA2(:,i-1))];
    C_k2=[C_k2 ((i-1)/i*C_k2(i-1)+1/i*NU_k2*NU_k2')];
    Q_k=Kk*C_k2(i)*Kk'/G'/G;
    Pk_k_12 = F*Pk*F'+G.*Q_k.*G.';
    Kk2=(Pk_k_12*H.')/(H*Pk_k_12*H.'+R);
    XA2=[XA2 (F*XA2(:,(i-1))+Kk2*NU_k2(i))];
    Pk=(eye(3,3)-Kk2*H)*Pk_k_12;
    SAF_A2 = SAF_A2 + (X(3, i) - XA2(3, i)).^2;
end

SAF_A2 = SAF_A2 / 5000;
figure('Name','Ошибка по скорости','NumberTitle','off');
plot(t,Z(1,:),'k-',t,XA2(1,:),'r-',t,X(1,:),'m-')
grid on;
legend(['Z'],['XA2'],['X'],'FontSize',8);
figure('Name','Ошибка по углу','NumberTitle','off');
plot(t,XA2(2,:),'r-',t,X(2,:),'m-')
grid on;
legend(['XA2'],['X'],'FontSize',8);
figure('Name','Дрейф','NumberTitle','on');
plot(t,XA2(3,:),'r-',t,X(3,:),'m-')
legend(['XA2'],['X'],'FontSize',8);
grid on;