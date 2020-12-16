clear
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
t=(1:5000);
e=-1+2*rand;
Vk=round(e*100)/10000;
R=10^(-2);
Q(1,1)=[5*10^-12];
XA=[ 0; 0; 0];
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


%���������� ������ ������� ����
C_k=[0];
NU_k=[0];
R_k=[0];
SAF_A=0;
for i=2:5000
    Pk_k_1 = F*Pk*F'+G*Q*G';
    NU_k=[NU_k (Z(1,i)-H*XA(:,i-1))];
    C_k=[C_k ((i-1)/i*C_k(i-1)+1/i*NU_k*NU_k')];
    R_k=[R_k (C_k(i)-H*Pk_k_1*H')];
    if R_k(i)<0
        R_k(i)=0;
    end
    Kk=(Pk_k_1*H')/(H*Pk_k_1*H'+R_k(i));
    XA=[XA (F*XA(:,(i-1))+Kk*NU_k(i))];
    Pk=(eye(3,3)-Kk*H)*Pk_k_1;
    SAF_A = SAF_A + (X(3, i) - XA(3, i)).^2;
end


SAF_A = SAF_A / 5000;
figure('Name','������ �� ��������','NumberTitle','off');
plot(t,Z(1,:),'k-',t,XA(1,:),'b-',t,X(1,:),'m-')
grid on;
legend(['Z'],['XA'],['X'],'FontSize',8);
figure('Name','������ �� ����','NumberTitle','off');
plot(t,XA(2,:),'b-',t,X(2,:),'m-')
grid on;
legend(['XA'],['X'],'FontSize',8);
figure('Name','�����','NumberTitle','on');
plot(t,XA(3,:),'b-',t,X(3,:),'m-')
legend(['XA'],['X'],'FontSize',8);
grid on;