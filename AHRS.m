%% ������
%{
���룺data��
1,2,3��Ϊ���������,ʱ��,���
4,5,6��Ϊaccx,accy,accz
7,8,9��Ϊgyrx,gyry,gyrz
10,11,12��Ϊmagx,magx,magz
13,14,15��Ϊpitch,roll,yaw
����������yaw��������pitch��������roll��������
%}
function [yaw,pitch,roll]=AHRS(data)

global normBc;global Index;global arrayBc;  % У׼�������������ж��Ƿ�У׼�ɹ�
normBc=0;Index=0;arrayBc=zeros(500,1);

time=data(:,2);

ax=data(:,4);ay=-data(:,5);az=-data(:,6);                         % accX,accY,accZ
gx=data(:,7)*pi/180;gy=-data(:,8)*pi/180;gz=-data(:,9)*pi/180;    % gyrX,gyrY,gyrZ,��λת��Ϊrad/s
mx=data(:,10);my=-data(:,11);mz=-data(:,12);                      % magX,magY,magZ
                                                                  % yzȡ����Ϊ��ת��Ϊ��ȷ������ϵ��
                                                                  % ��y����ǰ��x�����ң�z������
[datalength,~]=size(time);
yaw=ones(datalength,1);pitch=ones(datalength,1);roll=ones(datalength,1);

for i=1:datalength
    [mcx,mcy,mcz] = EKF(time(i),gx(i),gy(i),gz(i),mx(i),my(i),mz(i));
    isCalibrated=ifCalibrated();
    [yaw(i),pitch(i),roll(i)]=NineAxisFusion(time(i),ax(i),ay(i),az(i),gx(i),gy(i),gz(i),mcx,mcy,mcz,1);
    % ������һ����������ʵ�ʵ��㷨����
    % 0 ����
    % 1 ����
    % isCalibrated ��Ͼ���,���ж���У׼�ɹ����ʹ�þ���
    
	% ʹ��ԭʼ�ţ���ʹ��EKF����У׼
%     [yaw(i),pitch(i),roll(i)]=NineAxisFusion(time(i),ax(i),ay(i),az(i),gx(i),gy(i),gz(i),mx(i),my(i),mz(i),1);
end
end

%% �����ںϺ���
%{
���룺ʱ�䣬�������ݣ��ų�У׼�ɹ���־λ(isCalibrated)
����������yaw��������pitch��������roll��������
%}
function [yaw,pitch,roll]=NineAxisFusion(time,ax,ay,az,gx,gy,gz,mx,my,mz,isCalibrated)

%******** ���岢��ʼ����̬���� ********************************
persistent q_0 q_1 q_2 q_3;         %�����������Ԫ��
persistent curTime lastUpdate;      %��ǰʱ�� ��һʱ��
persistent dcGx dcGy dcGz;          %�������������ϵľ�̬ƫ��
persistent bStationary uStationaryCnt;      %�������Ƿ�ֹ�ı�־λ ��ֹ������
persistent norm normOld deltaNorm;          %��ǰģֵ �ϴ�ģֵ ģֵ��
persistent exInt eyInt ezInt;               %���ֻ����е��������
% persistent euler0 euler1 euler2;            %�ֱ��Ӧyaw pitch roll����Ϊ�����ã�һ�������Ĵ���������Ҫ�ⲿ��
if isempty(q_0)                             %��ʼ����̬����
    q_0=1;q_1=0;q_2=0;q_3=0;
    curTime=0;lastUpdate=0;
    dcGx=0;dcGy=0;dcGz=0;
    bStationary=0;uStationaryCnt=0;
    norm=0;normOld=0;deltaNorm=0;
    euler0 = 0;euler1 = 0;euler2 = 0;
    exInt=0;eyInt=0;ezInt=0;
end
%**************************************************************

% curTime=time*1e-9;                %��������time��λΪns������תΪs
curTime=time;
if lastUpdate==0                    %��һ�鴫����������ֻȡtime����lastUpdate����̬�Ƿ���Ϊ0���˳�����
    lastUpdate=curTime;
    yaw=0;pitch=0;roll=0;
    return
end
% if lastUpdate==curTime            %�����������֮�䴫��������û�и��£���̬�Ƿ�����һ�θ��µĽ�����˳�����
%     yaw=euler0;pitch=euler1;roll=euler2;
%     lastUpdate=curTime;
%     return
% end
halfT=0.5*(curTime-lastUpdate);     %����ʱ���һ�룬Ϊ�������㷽��
lastUpdate=curTime;
% halfT=0.01;

q0=q_0;q1=q_1;q2=q_2;q3=q_3;        %ȡ�ϴε����õ�����Ԫ����Ϊ���ε�������Ԫ��

normAcc=sqrt(ax*ax+ay*ay+az*az);    %���ٶ����ݹ�һ��
ax1=ax/normAcc;ay1=ay/normAcc;az1=az/normAcc;

%************* ��׽�����ǵľ�̬Ư�Ʋ����� ***********************
norm=sqrt(gx * gx + gy * gy + gz * gz);                     %��ǰ������ģֵ
deltaNorm = deltaNorm * 0.9 + abs(norm - normOld) * 0.1;    %��ǰģֵ����һʱ��ģֵ�Ĳ�ֵ���������˵�ͨ�˲�����
normOld = norm;
if deltaNorm<0.001                                          %�����ֵС��һ����ֵ
    if uStationaryCnt<100
        uStationaryCnt=uStationaryCnt+1;
    end
else
    bStationary=false;
    uStationaryCnt=0;
end
if uStationaryCnt==100          %�����ֵС����ֵ������һ��ʱ��
    bStationary=true;
    % isUsed=0;
end
if bStationary==true            %���¾�ֹƯ�ƣ�����ͬ�����˵�ͨ�˲�����
    dcGx=dcGx*0.99+gx*0.01;
    dcGy=dcGy*0.99+gy*0.01;
    dcGz=dcGz*0.99+gz*0.01;
end
gx1=gx-dcGx;
gy1=gy-dcGy;
gz1=gz-dcGz;
%*************************************************************

% �����ƹ�һ��
normB=sqrt(mx*mx+my*my+mz*mz);
mx1=mx/normB;my1=my/normB;mz1=mz/normB;

% �������Ʋ���ֵ�Ӵ���������ϵ��ת����������ϵ���õ�����ϵ�Ĵų�����
hx = 2 * mx1 * (0.5 - q2*q2 - q3*q3) + 2 * my1 * (q1*q2 - q0*q3) + 2 * mz1 * (q1*q3 + q0*q2);
hy = 2 * mx1 * (q1*q2 + q0*q3) + 2 * my1 * (0.5 - q1*q1 - q3*q3) + 2 * mz1 * (q2*q3 - q0*q1);
hz = 2 * mx1 * (q1*q3 - q0*q2) + 2 * my1 * (q2*q3 + q0*q1) + 2 * mz1 * (0.5 - q1*q1 - q2*q2);

inclination=atan(hz/sqrt(hx*hx+hy*hy));

%{
�����������ϵ�µĴų�ʸ��bxyz���ο�ֵ����
��Ϊ����ش�ˮƽ�нǣ�������֪��0�ȣ���ȥ��ƫ�ǵ����أ��̶��򱱣�������byָ������������by=ĳֵ��bx=0
������ο��ش�ʸ���ڴ�ֱ����Ҳ�з���bz��������ÿ���ط����ǲ�һ���ġ�
�����޷���֪��Ҳ���޷������ںϣ��и��ʺ�����ֱ���������ںϵļ��ٶȼƣ�������ֱ�ӴӲ���ֵhz�ϸ��ƹ�����bz=hz��
�ų�ˮƽ�������ο�ֵ�Ͳ���ֵ�Ĵ�СӦ����һ�µ�(bx*bx) + (by*by)) = ((hx*hx) + (hy*hy))��
��Ϊbx=0�����Ծͼ򻯳�(by*by)  = ((hx*hx) + (hy*hy))�������by�������޸�by��bxָ����Զ����ĸ���ָ������

�趨bx=0,��Y��Ϊǰ�����򣬵õ���һ�׵���ϵ�Ĵų�����bx��by��bz
%}

by = sqrt((hx * hx) + (hy * hy));
bz = hz;


%�����������ӵ�������ϵ��ת������������ϵ���õ�����ϵ�ļ��ٶȲο�ʸ��
vx=2*(q1*q3-q0*q2);
vy=2*(q0*q1+q2*q3);
vz=q0*q0-q1*q1-q2*q2+q3*q3;

%{
���ǰѵ�������ϵ�ϵĴų�ʸ��bxyz��ת����������wxyz
��Ϊbx=0�����������漰��bx�Ĳ��ֶ���ʡ���ˡ�ͬ��by=0�����������漰��by�Ĳ���Ҳ���Ա�ʡ�ԣ�������Լ������Ǹ���ָ���й�
������������vxyz�����㣬��Ϊ����g��az=1��ax=ay=0�����������漰��gx,gy�Ĳ���Ҳ��ʡ����
���Կ���������ʽ��wxyz�Ĺ�ʽ����by����ay��0������bz����az��1�����ͱ����vxyz�Ĺ�ʽ�ˣ�����q0q0+q1q1+q2q2+q3q3=1��
�趨by=0,��x��Ϊǰ������
wx = 2*bx*(0.5 - q2q2 - q3q3) + 2*bz*(q1q3 - q0q2);
wy = 2*bx*(q1q2 - q0q3) + 2*bz*(q0q1 + q2q3);
wz = 2*bx*(q0q2 + q1q3) + 2*bz*(0.5 - q1q1 - q2q2);
�趨bx=0,��Y��Ϊǰ�����򣬵õ�����ϵ�Ĵų�����
%}
wx = 2 * by * (q1*q2 + q0*q3) + 2 * bz * (q1*q3 - q0*q2);
wy = 2 * by * (0.5 - q1*q1 - q3*q3) + 2 * bz * (q0*q1 + q2*q3);
wz = 2 * by * (q2*q3 - q0*q1) + 2 * bz * (0.5 - q1*q1 - q2*q2);

%isCalibrated���жϴų��Ƿ��Ѿ�У׼�ɹ��ı�־λ��Ϊ��ʱ��ʹ�ôӴų��õ�������
if isCalibrated
    isUsed=1;
else
    isUsed=0;
end

%���������57��51'��57.85*pi/180=1.0097
% if inclination<0.985 || inclination>1.02
%     isUsed=0;
% end

%�Ӽ��ٶȼƵõ��Ŀ������
ex=ay1*vz-az1*vy;
ey=az1*vx-ax1*vz;
ez=ax1*vy-ay1*vx;

%�Ӵ����Ƶõ��Ŀ�����
exm = (my1 * wz - mz1 * wy)*isUsed;
eym = (mz1 * wx - mx1 * wz)*isUsed;
ezm = (mx1 * wy - my1 * wx)*isUsed;

% ��������
vKp=exp(-2*abs(normAcc-1));
mKp=1;
% mKp=exp(-0.01*abs(normB-45));
Ki=0;   %-0.001

% ���ֻ���
exInt=exInt+exm*Ki*halfT;
eyInt=eyInt+eym*Ki*halfT;
ezInt=ezInt+ezm*Ki*halfT;

% ����
refx=1.000*gx1+vKp*ex+mKp*exm+exInt;
refy=1.000*gy1+vKp*ey+mKp*eym+eyInt;
refz=1.000*gz1+vKp*ez+mKp*ezm+ezInt;

% ��Ԫ������
% qt=q(t-1)+[-1/2*(wt x q(t-1)]*dt, w is p([0 wx wy wz])	�������ݿɲο����ġ�Keeping a...��page19316
temp_q0 = q0 + (-q1 * refx - q2 * refy - q3 * refz) * halfT;
temp_q1 = q1 + (q0 * refx + q2 * refz - q3 * refy) * halfT;
temp_q2 = q2 + (q0 * refy - q1 * refz + q3 * refx) * halfT;
temp_q3 = q3 + (q0 * refz + q1 * refy - q2 * refx) * halfT;

% ��Ԫ����һ�����Է�����һ������
normq=sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
q0 = temp_q0 /normq;
q1 = temp_q1 /normq;
q2 = temp_q2 /normq;
q3 = temp_q3 /normq;

q_0=q0;q_1=q1;q_2=q2;q_3=q3;

%�ڼ���ŷ����ǰ����Ԫ�����������ŷ���ǵ�roll ��pitch����Ӱ��
q1=-q1;q2=-q2;q3=-q3;

%����ŷ���ǣ�����57.295780��Ϊ�˽�����ת��Ϊ�Ƕ�
yaw = atan2(2 * q1 * q2 - 2 * q0 * q3, 2 * q0 * q0 + 2 * q1 * q1 - 1) * 57.2957805;
if yaw < 0
    yaw = yaw + 360;
end
if yaw > 360
    yaw = yaw - 360;
end
pitch = -asin(2 * q1 * q3 + 2 * q0 * q2) * 57.2957805;
roll = atan2(2 * q2 * q3 - 2 * q0 * q1, 2 * q0 * q0 + 2 * q3 * q3 - 1) * 57.2957805;
% euler0=yaw;euler1=pitch;euler2=roll;    %��ŷ���Ǵ��뾲̬�������Է���dt=0ʱȡ��
% if time==210.7
%     normB
% end
end

%% EKF�ų�У׼����
%{
���룺ʱ�䣬���������ݣ��ų�����
�����У׼�Ĵų���������
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

%���岢��ʼ����ž���W Ӳ�ž���V У׼�ų�����Bc ״̬���̷���P
persistent W V Bc P;
if isempty(W)
    W = eye(3);
    V = [0 0 0]';
    Bc = [mx,my,mz]';
    P = [500*eye(3) zeros(3,6)  zeros(3);
         zeros(6,3) 1e-4*eye(6) zeros(6,3);
         zeros(3)   zeros(3,6)  500*eye(3)];
end

% curTime=time*1e-9;                  %�����time��λΪns������תΪs
curTime=time;
if lastUpdate==0                    %��һ�鴫����������ֻȡtime����lastUpdate��У׼ֵ����Ϊ0���˳�����
    lastUpdate=curTime;
    mcx=mx;mcy=my;mcz=mz;
    Bc=[mcx,mcy,mcz]';
    return
end
dt=curTime-lastUpdate;              %����ʱ��
lastUpdate=curTime;

% if dt==0                            %�����������֮�䴫��������û�и��£�У׼ֵ������һ�θ��µĽ�����˳�����
%     mcx=Bc(1);mcy=Bc(2);mcz=Bc(3);
%     return;
% end

% dt=0.02;

% ��׽�����ǵľ�̬Ư�Ʋ���������SixAxisFusion��NineAxisFusion�еĲ�����ͬ
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

% ������ת����Ak
w = [0 -gz1 gy1;
     gz1 0 -gx1;
     -gy1 gx1 0];                   % ���ٶȾ���
Ak = eye(3);
Ak = Ak + w * Ak * dt;
Ak = matrix_normal(Ak);             %��ÿһ�й�һ��
Ak = Ak';                           %������ת��

% ״̬һ��ת�ƣ�h_p_preΪ�ų����ֵ�����ֵ��X_preΪ״̬������ֵ
h_p_pre   = Ak * Bc;
X_pre     = [h_p_pre; W(1,1); W(2,2); W(3,3);
             W(1,2); W(1,3); W(2,3); V];

% ״̬���̵ķ��������������������Ϊ���鶨���
phi=10;                 %���Ʋ�������ʾ�������ǵ����ó̶ȣ�����أ���������Ҫ�޸�
gyrNoise=0.1/180*pi;     %�����ǵ�������������Ҫ�޸�
Q_element = dt*phi*gyrNoise*[(Bc(2)+Bc(3)) 0 0
                            0 (Bc(1)+Bc(3)) 0
                            0 0 (Bc(1)+Bc(2))];
Qk        = [Q_element'*Q_element zeros(3,9)     %����Q_element^2��Q_element'*Q_element��һ����
             zeros(9,3)  zeros(9,9)];
% ״̬ת�ƾ���
Fk        = [Ak zeros(3,9);
             zeros(9,3) eye(9,9)]; 
% �������
P_pre     = Fk * P * Fk' + Qk;

% �в�
yk        = [mx my mz]' - W * h_p_pre - V;

% ����Hk
hWk       = [h_p_pre(1) 0 0 h_p_pre(2) h_p_pre(3) 0;
             0 h_p_pre(2) 0 h_p_pre(1) 0 h_p_pre(3);
             0 0 h_p_pre(3) 0 h_p_pre(1) h_p_pre(2)]; 
hk        = [W hWk eye(3)];

% ���ⷽ�����
magNoise=5;          %������������������Ҫ�޸�
Rk        = magNoise^2 * eye(3);

% ��������Э����
Sk        = hk * P_pre * hk' + Rk;

% ����������
Kk        = P_pre * hk' / Sk;

% ״̬���������
X_now     = X_pre + Kk * yk;

% �������Э����
P_now     = (eye(12) - Kk * hk) * P_pre;        % ������Pk�ļ��㹫ʽ�Ǵ�ģ�û�����ߵ�Qk-1��

% ������
h_p_cal   = X_now(1:3);
W_last    = [X_now(4) X_now(7) X_now(8)      
             X_now(7) X_now(5) X_now(9)
             X_now(8) X_now(9) X_now(6)];
% ��һ��
balance   = abs(W_last(1,1) * W_last(2,2) * W_last(3,3))^(1/3);     % balance = trace(W_last)/3;
W = W_last / balance;
Bc = h_p_cal * balance;
V = X_now(10:12);
P = P_now;
mcx=Bc(1);mcy=Bc(2);mcz=Bc(3);

% ��У׼�ų������ݴ������飬�Թ�ifCalibrated()�жϴų��Ƿ��Ѿ�У׼���
normBc=sqrt(mcx*mcx+mcy*mcy+mcz*mcz);
arrayBc(Index+1)=normBc;
Index=rem((Index+1),501);

end

%% У׼����жϺ���
%{
�����У׼�ɹ���־
%}
function isCalibrated=ifCalibrated()
global normBc;global arrayBc;
persistent RMSE d2;     %RMSEΪ��������� d2Ϊ���㸨������
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
d2=0;           % ��Ҫ���˰Ѿ�̬����d2���㣡����
end
isCalibrated=false;

%�������ǹ涨��RMSEС��һ����ֵʱ�жϴų��Ѿ�У׼���ˣ�����Ӧ�ø���ʵ�����޸���ֵ���������������
if RMSE<2
    isCalibrated=true;
end
end

%% �����й�һ������
%{
���룺����
������й�һ������
%}
function A = matrix_normal(B)
A = zeros(3,3);
A(:,1) = B(:,1) / norm(B(:,1));
A(:,2) = B(:,2) / norm(B(:,2));
A(:,3) = B(:,3) / norm(B(:,3));
end


