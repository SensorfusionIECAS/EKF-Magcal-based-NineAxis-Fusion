%% 主函数
%{
输入：data：
1,2,3列为传感器序号,时间,序号
4,5,6列为accx,accy,accz
7,8,9列为gyrx,gyry,gyrz
10,11,12列为magx,magx,magz
13,14,15列为pitch,roll,yaw
输出：航向角yaw，俯仰角pitch，翻滚角roll的列向量
%}
function [yaw,pitch,roll]=AHRS(data)

global normBc;global Index;global arrayBc;  % 校准磁容器，用于判定是否校准成功
normBc=0;Index=0;arrayBc=zeros(500,1);

time=data(:,2);

ax=data(:,4);ay=-data(:,5);az=-data(:,6);                         % accX,accY,accZ
gx=data(:,7)*pi/180;gy=-data(:,8)*pi/180;gz=-data(:,9)*pi/180;    % gyrX,gyrY,gyrZ,单位转换为rad/s
mx=data(:,10);my=-data(:,11);mz=-data(:,12);                      % magX,magY,magZ
                                                                  % yz取负是为了转换为正确的坐标系，
                                                                  % 即y正向朝前，x正向朝右，z正向朝上
[datalength,~]=size(time);
yaw=ones(datalength,1);pitch=ones(datalength,1);roll=ones(datalength,1);

for i=1:datalength
    [mcx,mcy,mcz] = EKF(time(i),gx(i),gy(i),gz(i),mx(i),my(i),mz(i));
    isCalibrated=ifCalibrated();
    [yaw(i),pitch(i),roll(i)]=NineAxisFusion(time(i),ax(i),ay(i),az(i),gx(i),gy(i),gz(i),mcx,mcy,mcz,1);
    % 这里，最后一个参数控制实际的算法类型
    % 0 六轴
    % 1 九轴
    % isCalibrated 混合九轴,即判定磁校准成功后才使用九轴
    
	% 使用原始磁，不使用EKF进行校准
%     [yaw(i),pitch(i),roll(i)]=NineAxisFusion(time(i),ax(i),ay(i),az(i),gx(i),gy(i),gz(i),mx(i),my(i),mz(i),1);
end
end

%% 九轴融合函数
%{
输入：时间，九轴数据，磁场校准成功标志位(isCalibrated)
输出：航向角yaw，俯仰角pitch，翻滚角roll的列向量
%}
function [yaw,pitch,roll]=NineAxisFusion(time,ax,ay,az,gx,gy,gz,mx,my,mz,isCalibrated)

%******** 定义并初始化静态变量 ********************************
persistent q_0 q_1 q_2 q_3;         %参与迭代的四元数
persistent curTime lastUpdate;      %当前时刻 上一时刻
persistent dcGx dcGy dcGz;          %陀螺仪三个轴上的静态偏差
persistent bStationary uStationaryCnt;      %陀螺仪是否静止的标志位 静止计数器
persistent norm normOld deltaNorm;          %当前模值 上次模值 模值差
persistent exInt eyInt ezInt;               %积分环节中的误差向量
% persistent euler0 euler1 euler2;            %分别对应yaw pitch roll，作为缓存用，一般正常的传感器不需要这部分
if isempty(q_0)                             %初始化静态变量
    q_0=1;q_1=0;q_2=0;q_3=0;
    curTime=0;lastUpdate=0;
    dcGx=0;dcGy=0;dcGz=0;
    bStationary=0;uStationaryCnt=0;
    norm=0;normOld=0;deltaNorm=0;
    euler0 = 0;euler1 = 0;euler2 = 0;
    exInt=0;eyInt=0;ezInt=0;
end
%**************************************************************

% curTime=time*1e-9;                %如果输入的time单位为ns，这里转为s
curTime=time;
if lastUpdate==0                    %第一组传进来的数据只取time传入lastUpdate，姿态角返回为0，退出函数
    lastUpdate=curTime;
    yaw=0;pitch=0;roll=0;
    return
end
% if lastUpdate==curTime            %如果两次运算之间传感器数据没有更新，姿态角返回上一次更新的结果，退出函数
%     yaw=euler0;pitch=euler1;roll=euler2;
%     lastUpdate=curTime;
%     return
% end
halfT=0.5*(curTime-lastUpdate);     %采样时间的一半，为后续计算方便
lastUpdate=curTime;
% halfT=0.01;

q0=q_0;q1=q_1;q2=q_2;q3=q_3;        %取上次迭代得到的四元数作为本次迭代的四元数

normAcc=sqrt(ax*ax+ay*ay+az*az);    %加速度数据归一化
ax1=ax/normAcc;ay1=ay/normAcc;az1=az/normAcc;

%************* 捕捉陀螺仪的静态漂移并消除 ***********************
norm=sqrt(gx * gx + gy * gy + gz * gz);                     %当前陀螺仪模值
deltaNorm = deltaNorm * 0.9 + abs(norm - normOld) * 0.1;    %当前模值与上一时刻模值的差值，这里做了低通滤波处理
normOld = norm;
if deltaNorm<0.001                                          %如果差值小于一个阈值
    if uStationaryCnt<100
        uStationaryCnt=uStationaryCnt+1;
    end
else
    bStationary=false;
    uStationaryCnt=0;
end
if uStationaryCnt==100          %如果差值小于阈值持续了一段时间
    bStationary=true;
    % isUsed=0;
end
if bStationary==true            %更新静止漂移，这里同样做了低通滤波处理
    dcGx=dcGx*0.99+gx*0.01;
    dcGy=dcGy*0.99+gy*0.01;
    dcGz=dcGz*0.99+gz*0.01;
end
gx1=gx-dcGx;
gy1=gy-dcGy;
gz1=gz-dcGz;
%*************************************************************

% 磁力计归一化
normB=sqrt(mx*mx+my*my+mz*mz);
mx1=mx/normB;my1=my/normB;mz1=mz/normB;

% 将磁力计测量值从传感器坐标系旋转到地球坐标系，得到地球系的磁场向量
hx = 2 * mx1 * (0.5 - q2*q2 - q3*q3) + 2 * my1 * (q1*q2 - q0*q3) + 2 * mz1 * (q1*q3 + q0*q2);
hy = 2 * mx1 * (q1*q2 + q0*q3) + 2 * my1 * (0.5 - q1*q1 - q3*q3) + 2 * mz1 * (q2*q3 - q0*q1);
hz = 2 * mx1 * (q1*q3 - q0*q2) + 2 * my1 * (q2*q3 + q0*q1) + 2 * mz1 * (0.5 - q1*q1 - q2*q2);

inclination=atan(hz/sqrt(hx*hx+hy*hy));

%{
计算地理坐标系下的磁场矢量bxyz（参考值）。
因为地理地磁水平夹角，我们已知是0度（抛去磁偏角的因素，固定向北），定义by指向正北，所以by=某值，bx=0
但地理参考地磁矢量在垂直面上也有分量bz，地球上每个地方都是不一样的。
我们无法得知，也就无法用来融合（有更适合做垂直方向修正融合的加速度计），所以直接从测量值hz上复制过来，bz=hz。
磁场水平分量，参考值和测量值的大小应该是一致的(bx*bx) + (by*by)) = ((hx*hx) + (hy*hy))。
因为bx=0，所以就简化成(by*by)  = ((hx*hx) + (hy*hy))。可算出by。这里修改by和bx指向可以定义哪个轴指向正北

设定bx=0,以Y轴为前进方向，得到另一套地球系的磁场向量bx，by，bz
%}

by = sqrt((hx * hx) + (hy * hy));
bz = hz;


%把重力向量从地球坐标系旋转到传感器坐标系，得到载体系的加速度参考矢量
vx=2*(q1*q3-q0*q2);
vy=2*(q0*q1+q2*q3);
vz=q0*q0-q1*q1-q2*q2+q3*q3;

%{
我们把地理坐标系上的磁场矢量bxyz，转到机体上来wxyz
因为bx=0，所以所有涉及到bx的部分都被省略了。同理by=0，所以所有涉及到by的部分也可以被省略，这根据自己定义那个轴指北有关
类似上面重力vxyz的推算，因为重力g的az=1，ax=ay=0，所以上面涉及到gx,gy的部分也被省略了
可以看看两个公式：wxyz的公式，把by换成ay（0），把bz换成az（1），就变成了vxyz的公式了（其中q0q0+q1q1+q2q2+q3q3=1）
设定by=0,以x轴为前进方向
wx = 2*bx*(0.5 - q2q2 - q3q3) + 2*bz*(q1q3 - q0q2);
wy = 2*bx*(q1q2 - q0q3) + 2*bz*(q0q1 + q2q3);
wz = 2*bx*(q0q2 + q1q3) + 2*bz*(0.5 - q1q1 - q2q2);
设定bx=0,以Y轴为前进方向，得到载体系的磁场向量
%}
wx = 2 * by * (q1*q2 + q0*q3) + 2 * bz * (q1*q3 - q0*q2);
wy = 2 * by * (0.5 - q1*q1 - q3*q3) + 2 * bz * (q0*q1 + q2*q3);
wz = 2 * by * (q2*q3 - q0*q1) + 2 * bz * (0.5 - q1*q1 - q2*q2);

%isCalibrated是判断磁场是否已经校准成功的标志位，为真时才使用从磁场得到的数据
if isCalibrated
    isUsed=1;
else
    isUsed=0;
end

%北京磁倾角57°51'，57.85*pi/180=1.0097
% if inclination<0.985 || inclination>1.02
%     isUsed=0;
% end

%从加速度计得到的控制误差
ex=ay1*vz-az1*vy;
ey=az1*vx-ax1*vz;
ez=ax1*vy-ay1*vx;

%从磁力计得到的控制误差，
exm = (my1 * wz - mz1 * wy)*isUsed;
eym = (mz1 * wx - mx1 * wz)*isUsed;
ezm = (mx1 * wy - my1 * wx)*isUsed;

% 比例环节
vKp=exp(-2*abs(normAcc-1));
mKp=1;
% mKp=exp(-0.01*abs(normB-45));
Ki=0;   %-0.001

% 积分环节
exInt=exInt+exm*Ki*halfT;
eyInt=eyInt+eym*Ki*halfT;
ezInt=ezInt+ezm*Ki*halfT;

% 误差补偿
refx=1.000*gx1+vKp*ex+mKp*exm+exInt;
refy=1.000*gy1+vKp*ey+mKp*eym+eyInt;
refz=1.000*gz1+vKp*ez+mKp*ezm+ezInt;

% 四元数积分
% qt=q(t-1)+[-1/2*(wt x q(t-1)]*dt, w is p([0 wx wy wz])	具体内容可参考论文《Keeping a...》page19316
temp_q0 = q0 + (-q1 * refx - q2 * refy - q3 * refz) * halfT;
temp_q1 = q1 + (q0 * refx + q2 * refz - q3 * refy) * halfT;
temp_q2 = q2 + (q0 * refy - q1 * refz + q3 * refx) * halfT;
temp_q3 = q3 + (q0 * refz + q1 * refy - q2 * refx) * halfT;

% 四元数归一化，以方便下一次运算
normq=sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
q0 = temp_q0 /normq;
q1 = temp_q1 /normq;
q2 = temp_q2 /normq;
q3 = temp_q3 /normq;

q_0=q0;q_1=q1;q_2=q2;q_3=q3;

%在计算欧拉角前对四元数做共轭，避免欧拉角的roll 和pitch互相影响
q1=-q1;q2=-q2;q3=-q3;

%计算欧拉角，乘以57.295780是为了将弧度转化为角度
yaw = atan2(2 * q1 * q2 - 2 * q0 * q3, 2 * q0 * q0 + 2 * q1 * q1 - 1) * 57.2957805;
if yaw < 0
    yaw = yaw + 360;
end
if yaw > 360
    yaw = yaw - 360;
end
pitch = -asin(2 * q1 * q3 + 2 * q0 * q2) * 57.2957805;
roll = atan2(2 * q2 * q3 - 2 * q0 * q1, 2 * q0 * q0 + 2 * q3 * q3 - 1) * 57.2957805;
% euler0=yaw;euler1=pitch;euler2=roll;    %将欧拉角存入静态变量，以方便dt=0时取用
% if time==210.7
%     normB
% end
end

%% EKF磁场校准函数
%{
输入：时间，陀螺仪数据，磁场数据
输出：校准的磁场数据向量
%}
function [mcx,mcy,mcz] = EKF(time,gx,gy,gz,mx,my,mz)
global normBc;global Index;global arrayBc;
persistent curTime lastUpdate;
persistent dcGx dcGy dcGz;
persistent bStationary uStationaryCnt;
persistent norm normOld deltaNorm;
if isempty(curTime)
    curTime=0;lastUpdate=0;
    dcGx=0;dcGy=0;dcGz=0;
    bStationary=0;uStationaryCnt=0;
    norm=0;normOld=0;deltaNorm=0;
end

%定义并初始化软磁矩阵W 硬磁矩阵V 校准磁场向量Bc 状态方程方差P
persistent W V Bc P;
if isempty(W)
    W = eye(3);
    V = [0 0 0]';
    Bc = [mx,my,mz]';
    P = [500*eye(3) zeros(3,6)  zeros(3);
         zeros(6,3) 1e-4*eye(6) zeros(6,3);
         zeros(3)   zeros(3,6)  500*eye(3)];
end

% curTime=time*1e-9;                  %输入的time单位为ns，这里转为s
curTime=time;
if lastUpdate==0                    %第一组传进来的数据只取time传入lastUpdate，校准值返回为0，退出函数
    lastUpdate=curTime;
    mcx=mx;mcy=my;mcz=mz;
    Bc=[mcx,mcy,mcz]';
    return
end
dt=curTime-lastUpdate;              %采样时间
lastUpdate=curTime;

% if dt==0                            %如果两次运算之间传感器数据没有更新，校准值返回上一次更新的结果，退出函数
%     mcx=Bc(1);mcy=Bc(2);mcz=Bc(3);
%     return;
% end

% dt=0.02;

% 捕捉陀螺仪的静态漂移并消除，与SixAxisFusion和NineAxisFusion中的部分相同
norm=sqrt(gx * gx + gy * gy + gz * gz);
deltaNorm = deltaNorm * 0.9 + abs(norm - normOld) * 0.1;
normOld = norm;
if deltaNorm<0.001
    if uStationaryCnt<100
        uStationaryCnt=uStationaryCnt+1;
    end
else
    bStationary=false;
    uStationaryCnt=0;
end
if uStationaryCnt==100
    bStationary=true;
end
if bStationary==true
    dcGx=dcGx*0.99+gx*0.01;
    dcGy=dcGy*0.99+gy*0.01;
    dcGz=dcGz*0.99+gz*0.01;
end
gx1=gx-dcGx;gy1=gy-dcGy;gz1=gz-dcGz;

% 计算旋转矩阵Ak
w = [0 -gz1 gy1;
     gz1 0 -gx1;
     -gy1 gx1 0];                   % 角速度矩阵
Ak = eye(3);
Ak = Ak + w * Ak * dt;
Ak = matrix_normal(Ak);             %把每一列归一化
Ak = Ak';                           %自身作转置

% 状态一步转移，h_p_pre为磁场部分的先验值，X_pre为状态量先验值
h_p_pre   = Ak * Bc;
X_pre     = [h_p_pre; W(1,1); W(2,2); W(3,3);
             W(1,2); W(1,3); W(2,3); V];

% 状态方程的方差，这里有两个参数是人为经验定义的
phi=10;                 %控制参数，表示对陀螺仪的利用程度（负相关），根据需要修改
gyrNoise=0.1/180*pi;     %陀螺仪的噪声，根据需要修改
Q_element = dt*phi*gyrNoise*[(Bc(2)+Bc(3)) 0 0
                            0 (Bc(1)+Bc(3)) 0
                            0 0 (Bc(1)+Bc(2))];
Qk        = [Q_element'*Q_element zeros(3,9)     %这里Q_element^2与Q_element'*Q_element是一样的
             zeros(9,3)  zeros(9,9)];
% 状态转移矩阵
Fk        = [Ak zeros(3,9);
             zeros(9,3) eye(9,9)]; 
% 方差矩阵
P_pre     = Fk * P * Fk' + Qk;

% 残差
yk        = [mx my mz]' - W * h_p_pre - V;

% 计算Hk
hWk       = [h_p_pre(1) 0 0 h_p_pre(2) h_p_pre(3) 0;
             0 h_p_pre(2) 0 h_p_pre(1) 0 h_p_pre(3);
             0 0 h_p_pre(3) 0 h_p_pre(1) h_p_pre(2)]; 
hk        = [W hWk eye(3)];

% 量测方差矩阵
magNoise=5;          %磁力计噪声，根据需要修改
Rk        = magNoise^2 * eye(3);

% 测量余量协方差
Sk        = hk * P_pre * hk' + Rk;

% 卡尔曼增益
Kk        = P_pre * hk' / Sk;

% 状态量后验估计
X_now     = X_pre + Kk * yk;

% 后验估计协方差
P_now     = (eye(12) - Kk * hk) * P_pre;        % 论文中Pk的计算公式是错的，没有最后边的Qk-1。

% 输出结果
h_p_cal   = X_now(1:3);
W_last    = [X_now(4) X_now(7) X_now(8)      
             X_now(7) X_now(5) X_now(9)
             X_now(8) X_now(9) X_now(6)];
% 归一化
balance   = abs(W_last(1,1) * W_last(2,2) * W_last(3,3))^(1/3);     % balance = trace(W_last)/3;
W = W_last / balance;
Bc = h_p_cal * balance;
V = X_now(10:12);
P = P_now;
mcx=Bc(1);mcy=Bc(2);mcz=Bc(3);

% 将校准磁场的数据存入数组，以供ifCalibrated()判断磁场是否已经校准完成
normBc=sqrt(mcx*mcx+mcy*mcy+mcz*mcz);
arrayBc(Index+1)=normBc;
Index=rem((Index+1),501);

end

%% 校准完成判断函数
%{
输出：校准成功标志
%}
function isCalibrated=ifCalibrated()
global normBc;global arrayBc;
persistent RMSE d2;     %RMSE为均方根误差 d2为计算辅助变量
if isempty(RMSE)
    RMSE=5;
    d2=0;
end
isnZero=true;
for i=1:500
    isnZero=isnZero&&(arrayBc(i)~=0);
end
if isnZero
    for i=1:500
        d2=d2+(arrayBc(i)-normBc)^2;
    end
RMSE=0.2*sqrt(d2/500)+0.8*RMSE;
d2=0;           % 不要忘了把静态变量d2归零！！！
end
isCalibrated=false;

%这里我们规定当RMSE小于一个阈值时判断磁场已经校准好了，这里应该根据实验来修改阈值，或添加其他条件
if RMSE<2
    isCalibrated=true;
end
end

%% 矩阵列归一化函数
%{
输入：矩阵
输出：列归一化矩阵
%}
function A = matrix_normal(B)
A = zeros(3,3);
A(:,1) = B(:,1) / norm(B(:,1));
A(:,2) = B(:,2) / norm(B(:,2));
A(:,3) = B(:,3) / norm(B(:,3));
end


