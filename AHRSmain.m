%% ������
%{
LPMSB2_3.csv��ʽ:
1,2,3��Ϊ���������,ʱ��,���
4,5,6��Ϊaccx,accy,accz
7,8,9��Ϊgyrx,gyry,gyrz
10,11,12��Ϊmagx,magx,magz
13,14,15��Ϊpitch,roll,yaw
%}
clear;
clear AHRS NineAxisFusion EKF ifCalibrated;         %��պ����ڵľ�̬���� ��������Ҫ������
% data=load ('adata.txt');      %��txt��ʽ
% data=xlsread('adata.xls');    %��excel��ʽ
data=load('LPMSB2_3.csv');      %��csv��ʽ
[x,y,z]=AHRS(data);             %[yaw,pitch,roll]
figure(9);
%plot(x(240:end)-38,'b');	%4.3
x1=x-22;    % ֮ǰ��-24.2��֪��Ϊʲô�ĳ���24.2   9����-23.9
subplot(2,1,1);plot(x1,'k');title('�����Ư�� ( ���ľ�����̬�ں��㷨 ) ')
%title('Deviation curve of heading angle ( nine-axis fusion )')
axis([0 62000 -10 30]); % 67100
xlabel('t/0.01s');
ylabel('Ư����/��');
%ylabel('deviation (truncated)/��');
hold on;grid on;
%plot(-data(240:end,15)+180,'r');   % 4.3
data1=-data(:,15)-92.14;    % ����Ϊʲô-88.12   9����-88.91
subplot(2,1,2);plot(data1,'k');title('�����Ư�� ( LPMS-B2 GYR+ACC )')
%title('Deviation curve of heading angle ( LPMS-B2 GYR+ACC+MAG )')
axis([0 62000 -10 30]); %67100
xlabel('t/0.01s');
ylabel('Ư����/��');
grid on;
%ylabel('deviation (truncated)/��');
%plot(-data(:,15)+180,'r');
%legend('Six-axis Fusion','LPMSB2(GYR+ACC+MAG)');
%axis([0 68000 0 360]);
%xlabel('t/0.01s');ylabel('Yaw/��');
