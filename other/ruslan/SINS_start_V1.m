%% ��������� ������
% �� = �*�(�-1)+G*W(k-1)
% Zk = H*X(k) + V(k)
clear
T=1;
A=0.0001;
g=9.81;
R=6371000;
bi=5*0.00000001;
W(1)=(randi(201)-101)/100;
F=[ 1 -g*T 0 ; T/R 1 T ; 0 0 (1-bi*T) ];
X(1,1)= [ 0 ];
X(2,1)= [ 0 ];
X(3,1)= [ 0 ];
G=[0 ; 0 ; T*A*sqrt(2*bi)];
H=[ 1  0  0; 0 0 1];
Z=[0 ; 0];

t=(1:5000);

for k=2:5000
    W=[W (randi(201)-101)/100];
    X=[X (F*X(:,(k-1))+G*W((k-1)))];
    pogr=[((randi(201)-101)/100) * 10^(-2); 
    ((randi(201)-101)/100) * 10^(-8)];
    Z=[Z H*X(:,k)+pogr];
    
end
%% ���������� ��������
figure('Name','������ �� ��������','NumberTitle','off');
plot(t,X(1,:),t,Z(1,:))
figure('Name','������ �� ����','NumberTitle','off');
plot(t,X(2,:))
figure('Name','�����','NumberTitle','on');
plot(t,X(3,:),t,Z(2,:))
%% ���������� ����� � �������
cd 'c:\Users\������\Desktop\������ ��������'
save('dannie_V1')