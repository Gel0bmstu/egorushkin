clear;
close('all');
T = 1;
A = 10^(-4);
g = 9.81;
R = 6371000;
betta = 5*10^(-8);

F = [1 -g*T 0; T/R 1 T; 0 0 (1-betta*T)];
G = [0; 0; T*A*sqrt(2*betta)];
H = [1 0 0];

number_of_steps = 5000;
t = (1:number_of_steps);

TAR = 1;
TAR_d = [0;0;0];

W = dlmread('W_values.txt');
X = dlmread('X_values.txt');
Z = dlmread('Z_values.txt');

% W = zeros(1,number_of_steps);
% W(1) = randi([-100,100])/100;
% X = zeros(3,number_of_steps);
% Z = zeros(1,number_of_steps);
% Vz = randi([-100,100])*1e-3;
% for k = 2:5000
%     W(k) = randi([-100,100])/100;
%     X(:,k) = F * X(:,k-1) + G * W(k-1);
%     Z(k) = H*X(:,k) + Vz;
% end

% figure('Name','�����','NumberTitle','on');
% plot(t,X(3,:));
% grid on;

dlmwrite('W_values.txt',W,'delimiter','\t','precision',6);
dlmwrite('X_values.txt',X,'delimiter','\t','precision',6);
dlmwrite('Z_values.txt',Z,'delimiter','\t','precision',6);


% -----------------------------------------------------------
% Simple kalman
R=10^(-2);
Q(1,1)=[3.3*10^-12];
Xkd=[ 0; 0; 0 ];
Pk=[ 1000 0 0 ; 0 10^(-5) 0 ; 0 0 10^(-12) ];
for i=2:5000
    Pk_k_1 = F*Pk*F'+Q;
    Kk=(Pk_k_1*H')/(H*Pk_k_1*H'+R);
    Xkd=[Xkd (F*Xkd(:,(i-1))+Kk*(Z(i)-H*F*Xkd(:,(i-1))))];
    Pk=(eye(3)-Kk*H)*Pk_k_1;
end

% Simple kalman
figure();
plot(t,X(3,:), t,Xkd(3, :));
grid on;
% xlim([t_b t_e]);
legend('X','FontSize',8);

figure();
plot(t,X(1,:), t,Xkd(1, :));
grid on;
% xlim([t_b t_e]);
legend('X','FontSize',8);

% -----------------------------------------------------------

% First-order adaptive filter
XA1 = zeros(3,number_of_steps);
Pk11 = [1e3 0 0; 0 1e-5 0; 0 0 1e-12];
Nu_1 = zeros(1,number_of_steps);
Ck_1 = zeros(1,number_of_steps);
Rk_1 = zeros(1,number_of_steps);
SAF_1 = 0;
Q = [0 0 0; 0 0 0; 0 0 3.3*10^-8];
for k = 2:number_of_steps
    Pkk11 = F * Pk11 * F.' + Q;
    Nu_1(k) = Z(1,k) - H*F*XA1(:,k-1);
    Ck_1(k) = (k-1)/k*Ck_1(k-1) + 1/k*Nu_1*Nu_1';
    if ((Ck_1(k) - H*Pkk11*H') > 0) 
        Rk_1(k) = Ck_1(k) - H*Pkk11*H';
    else
        Rk_1(k) =  0;
    end
    Kk = Pkk11*H.'/(H*Pkk11*H.' + Rk_1(k));
    XA1(:,k) = F*XA1(:,k-1) + Kk * Nu_1(k);
    Pk11 = (eye(3,3) - Kk*H)*Pkk11;
    SAF_1 = SAF_1 + (X(3, k) - XA1(3, k))^2;
    TAR = TAR * ((eye(3,3) - Kk*H)*F);
    TAR_d = [TAR_d [TAR(1,1);TAR(2,2);TAR(3,3)]];
end
SAF_1 = SAF_1 / 5000;

% ���������� ������ 2-�� ����
XA2 = zeros(3,number_of_steps);
Pk2 = [1e3 0 0; 0 1e-5 0; 0 0 1e-12];
R2 = 5*10^-5;
Nu_2 = zeros(1,number_of_steps);
Ck_2 = zeros(1,number_of_steps);
Rk2 = zeros(1,number_of_steps);
SAF_2 = 0;
Q_2 = [0];
K_2 = [0 0 1];
for k = 2:5000
    Nu_2(k) = Z(1,k) - H*F*XA2(:,k-1);
    Ck_2(k) = (k-1)/k*Ck_2(k-1) + 1/k*Nu_2*Nu_2';
    Q_2  = K_2 * Ck_2(k) * K_2';
    Pkk12 = F*Pk2*F.' + Q_2;
    K_2 = Pkk12*H.'/(H*Pkk12*H' + R2);
    XA2(:,k) = F*XA2(:,k-1) + K_2 * Nu_2(k);
    Pk2 = (eye(3,3) - K_2*H)*Pkk12;
    SAF_2 = SAF_2 + (X(3, k) - XA2(3, k))^2;
end
SAF_2 = SAF_2 / 5000;

t_b = 500;
t_e = 4000;
% figure('Name','������ �� ��������','NumberTitle','off');
% plot(t,X(1,:),'r-',t,XA1(1,:),'b-',t,XA2(1,:),'g-');
% grid on;
% legend(['X'],['XA'],['XA2'],'FontSize',8);
% % 
% figure('Name','������ �� ����','NumberTitle','off');
% plot(t,X(2,:),'r-',t,XA1(2,:),'b-',t,XA2(2,:),'g-');
% grid on;
% legend(['X'],['XA'],['XA2'],'FontSize',8);


Q_2 = [0 0 0; 0 0 0; 0 0 Q_2(3,3)*10^-2];
% 
% 
figure('Name','All adaptive','NumberTitle','on');
plot(t,X(3,:),'r-',t,XA1(3,:),'b-',t,XA2(3,:),'g-');
grid on;
xlim([t_b t_e]);
legend(['X'],['XA'],['XA2'],'FontSize',8);
% 
% %��������������, ���� ����� ������
% figure('Name','�����','NumberTitle','on');
% plot(t,X(3,:),'r-',t,XA1(3,:),'b-',t,XA2(3,:),'g-',t,Xkd(3,:));
% grid on;
% xlim([t_b t_e]);
% legend(['X'],['XA'],['XA2'],['Xkd'],'FontSize',8);

% figure('Name','TAR δVe измерение скорости','NumberTitle','off');
% plot(t,TAR_d(1,:));
% grid on;
% 
% figure('Name','TAR Фn измерение скорости','NumberTitle','off');
% plot(t,TAR_d(2,:));
% grid on;
% 
% figure('Name','TAR ωдр измерение скорости','NumberTitle','off');
% plot(t,TAR_d(3,:));
% grid on;



