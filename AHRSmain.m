%% 运行我
%{
LPMSB2_3.csv格式:
1,2,3列为传感器序号,时间,序号
4,5,6列为accx,accy,accz
7,8,9列为gyrx,gyry,gyrz
10,11,12列为magx,magx,magz
13,14,15列为pitch,roll,yaw
%}
clear;
clear AHRS NineAxisFusion EKF ifCalibrated;         %清空函数内的静态变量 ！！！必要！！！
% data=load ('adata.txt');      %读txt格式
% data=xlsread('adata.xls');    %读excel格式
data=load('LPMSB2_3.csv');      %读csv格式
[x,y,z]=AHRS(data);             %[yaw,pitch,roll]
figure(9);
%plot(x(240:end)-38,'b');	%4.3
x1=x-22;    % 之前是-24.2不知道为什么改成了24.2   9轴是-23.9
subplot(2,1,1);plot(x1,'k');title('航向角漂移 ( 本文九轴姿态融合算法 ) ')
%title('Deviation curve of heading angle ( nine-axis fusion )')
axis([0 62000 -10 30]); % 67100
xlabel('t/0.01s');
ylabel('漂移量/°');
%ylabel('deviation (truncated)/°');
hold on;grid on;
%plot(-data(240:end,15)+180,'r');   % 4.3
data1=-data(:,15)-92.14;    % 忘了为什么-88.12   9轴是-88.91
subplot(2,1,2);plot(data1,'k');title('航向角漂移 ( LPMS-B2 GYR+ACC )')
%title('Deviation curve of heading angle ( LPMS-B2 GYR+ACC+MAG )')
axis([0 62000 -10 30]); %67100
xlabel('t/0.01s');
ylabel('漂移量/°');
grid on;
%ylabel('deviation (truncated)/°');
%plot(-data(:,15)+180,'r');
%legend('Six-axis Fusion','LPMSB2(GYR+ACC+MAG)');
%axis([0 68000 0 360]);
%xlabel('t/0.01s');ylabel('Yaw/°');
