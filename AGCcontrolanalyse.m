% clear
load('SHDdata.mat')
load('SHDPCSdata1228deal')
data=SHDdata.data1228(:,4:6);
Agc=data(1:2:end,1);% AGC指令
P=data(1:2:end,2);% 机组出力
% Pbat=SDPCSdatadeal.all(:,37);
% Pall=Pbat+P;
Pall=data(1:2:end,3);% 联合出力
LineMax=length(Agc);% 一天的数据点的时长
Result=zeros(10,12);% 结果保存数据
% Result=[1         2        3          4          5      6         7          8              9               10             11            12]
% Result=[AGC 机组起始时刻值 t0 机组首次开始动作值 t1 反向记录  不调记录  动作缓慢记录  瞬间调节完成  调节效果（有1,无0）调节方向  机组本身到达死区]
ctrlNo=1;
Result(ctrlNo,1)=Agc(1);
Result(ctrlNo,2)=P(1);
Result(ctrlNo,3)=1;
detAGC=2;
Prate=600;
Erate=12;
T1=0;% 反调计时
Cdead=0.02;
T2i=0;% 不调起始时刻，为了保证连续一段时间内均是基本不动作
T2=0;% 不调计时
Tlen=40;
if Result(ctrlNo,1)>Result(ctrlNo,2)
    Result(ctrlNo,11)=1;
else
    Result(ctrlNo,11)=-1;
end
for i=1:LineMax
    if (Agc(i) > detAGC+Result(ctrlNo,1)) ||  (Agc(i) < Result(ctrlNo,1)-detAGC)
        % 新来一条指令，核算上调指令的情况
        if Result(ctrlNo,6)==0 && Result(ctrlNo,7)==0 
            % 没有反调或不调,以及机组未到达死区的情况下，计算机组的调节速率是否缓慢
            if (P(i)-Result(ctrlNo,4))/(i-Result(ctrlNo,5))*60<0.015*Prate
                % 若机组调节速率小于标准速率
                if Result(ctrlNo,12)==0
                    % 且机组没有到达调节死区范围内
                    Result(ctrlNo,8)=1;
                end
            end
        end
%         if abs(Pall(i)-Result(ctrlNo,1))<Cdead
%             % 联合出力到达调节范围内，记调节成功一次
%             Result(ctrlNo,10)=1;
%         end
        if Result(ctrlNo,5)==0 
            Result(ctrlNo,7)=1;
        end
        if i-Result(ctrlNo,3)>Tlen
            % 指令持续时长小于40秒，指令作废
            ctrlNo=ctrlNo+1;
            Result(ctrlNo,1)=Agc(i);
            Result(ctrlNo,2)=P(i);
            Result(ctrlNo,3)=i;
            if Result(ctrlNo,1)>Result(ctrlNo,2)
                Result(ctrlNo,11)=1;
            else
                Result(ctrlNo,11)=-1;
            end
            T1=0;
            T2i=0;
            T2=0;
        else
            Result(ctrlNo,:)=0;
            Result(ctrlNo,1)=Agc(i);
            Result(ctrlNo,2)=P(i);
            Result(ctrlNo,3)=i;
            if Result(ctrlNo,1)>Result(ctrlNo,2)
                Result(ctrlNo,11)=1;
            else
                Result(ctrlNo,11)=-1;
            end
            T1=0;
            T2i=0;
            T2=0;
        end
    else
        if abs(P(i)-Result(ctrlNo,2))>0.3 && Result(ctrlNo,5)==0
            % 判断机组开始响应指令时刻
            Result(ctrlNo,4)=P(i);
            Result(ctrlNo,5)=i;
        end
        if Result(ctrlNo,6)==0
            % 机组反向调节
            if i>2 && (P(i)-P(i-1))*(Result(ctrlNo,11))<0
                % 机组出力方向与指令方向相反
                T1=T1+1;
                if T1>=10 && (P(i)-P(i-T1))*Result(ctrlNo,11)<-2
                    % 持续10s即判断是反调，且始末的值要超过-2MW
                    Result(ctrlNo,6)=1;
                    T1=0;
                end
            end
        end
        if Result(ctrlNo,7)==0 && Result(ctrlNo,5)~=0
            % 机组不调记录
            if abs(Result(ctrlNo,1)-P(i))>Cdead*Result(ctrlNo,1)
                % 机组出力不在调节死区内才正常计数
                if T2i==0
                    if i>2 && abs(P(i)-P(i-1))<0.2
                        % 记录不调起始时刻，并开始计数
                        T2i=i;
                        T2=T2+1;
                    end
                else
                    if i>2 && abs(P(i)-P(i-1))<0.2 && i-T2i==T2
                        % 机组波动死区0.2MW,以免有小幅的波动，导致连续性破坏
                        T2=T2+1;
                        if T2>=Tlen/2 && abs(P(i)-P(i-T2))<0.4
                            % 达到20s后，若此时机组出力相较于起始出力仍在0.4MW内，则记为不调
                            Result(ctrlNo,7)=1;
                            T2=0;
                            T2i=0;
                        end
                    else
                        T2=0;
                        T2i=0;
                    end
                end
            end
        end
        if i-Result(ctrlNo,3)<=5
            % 3s内
            if abs(Pall(i)-Result(ctrlNo,1))<1  % && abs(Pall(i)-P(i))>Erate/2
                % 储能将出力拉至AGC处，并且功率大于0.5倍额定功率
                Result(ctrlNo,9)=1;
            end
        end
        if Result(ctrlNo,12)==0 && abs(Result(ctrlNo,1)-P(i))<Cdead*Result(ctrlNo,1)
            Result(ctrlNo,12)=1;
        end
        if abs(Pall(i)-Result(ctrlNo,1))<Cdead*Result(ctrlNo,1) && Result(ctrlNo,10)==0
            % 联合出力到达调节范围内，记调节成功一次
            Result(ctrlNo,10)=1;
        end
    end
end
%% 分析
N=length(Result(:,6));
M=0;
Op1=sum(Result(:,6))/N*100 % 机组反调比例
Op2=sum(Result(:,7))/N*100 % 机组长时间不动的比例
Op3=sum(Result(:,8))/N*100 % 机组调节速度小于标准调节速度的比例
Op4=sum(Result(:,9))/N*100 % 3s内将机组出力拉至AGC指令上
S=0;
for i=1:N
    if Result(i,6)==1 || Result(i,7)==1 || Result(i,8)==1
        M=M+1;
        if Result(i,10)==1
            S=S+1;
        end
    end
end
S %代表当前储能调节策略下，在机组出力存在问题的情况下，有效调节（可获Kp值）的次数
Sop=S/M*100  % 相应比例
MchuN=M/N*100 % 无效指令
%% 储能统计
Pdg=Pall-P;
M_point5=0;
M_1=0;
M_1point5=0;
M_2=0;
E_point5=0; 
E_1=0;
E_1point5=0;
E_2=0;
E_point5plus=0;
E_point5minus=0;
E_1plus=0;
E_1minus=0;
E_1point5plus=0;
E_1point5minus=0;
E_2plus=0;
E_2minus=0;
ResultPdg=zeros(10,5);
ctrl=1;
T=0;
T_time=0;
for i=1:length(Pdg)
    if abs(Pdg(i))>0  && abs(Pdg(i))<4.5
        M_point5=M_point5+1;
        E_point5=E_point5+abs(Pdg(i))*1/3600;% 按0.5C等效
        if Pdg(i)>0
            E_point5plus=E_point5plus+Pdg(i)*1/3600;
        else
            E_point5minus=E_point5minus+Pdg(i)*1/3600;
        end
    end
    if abs(Pdg(i))>4.5  && abs(Pdg(i))<9
        M_1=M_1+1;
        E_1=E_1+abs(Pdg(i))*1/3600;
        if Pdg(i)>0
            E_1plus=E_1plus+Pdg(i)*1/3600;% 按1C等效
        else
            E_1minus=E_1minus+Pdg(i)*1/3600;
        end
    end
    if abs(Pdg(i))>9  && abs(Pdg(i))<13.5
        M_1point5=M_1point5+1;
        E_1point5=E_1point5+abs(Pdg(i))*1/3600;
        if Pdg(i)>0
            E_1point5plus=E_1point5plus+Pdg(i)*1/3600; % 按1.5C等效
        else
            E_1point5minus=E_1point5minus+Pdg(i)*1/3600;
        end
    end
    if abs(Pdg(i))>13.5  && abs(Pdg(i))<18
        M_2=M_2+1;
        E_2=E_2+abs(Pdg(i))*1/3600;
        if Pdg(i)>0
            E_2plus=E_2plus+Pdg(i)*1/3600;% 按2C等效
        else
            E_2minus=E_2minus+Pdg(i)*1/3600;
        end
    end
    if abs(Pdg(i))>13.5
        if T==0
            ResultPdg(ctrl,1)=Pdg(i);% 起始时刻的功率
            ResultPdg(ctrl,2)=i;% 长输出的起始时刻
            ResultPdg(ctrl,5)=abs(Pdg(i))/3600;% 累计能量
            T=T+1;
        else
            ResultPdg(ctrl,5)=ResultPdg(ctrl,5)+abs(Pdg(i))/3600;
        end
        T_time=0;
    else
        T_time=T_time+1;
        if T_time<60
            % 连续60s以内低于2C的功率，不算中断2C
%             ResultPdg(ctrl,5)=ResultPdg(ctrl,5)+Pdg(i)/3600;
        else
            % 连续60s大于2C功率算中断2C持续时长
            if T~=0
                ResultPdg(ctrl,3)=i-T_time;% 刨除60s，计算长指令结束的时间
                ResultPdg(ctrl,4)=ResultPdg(ctrl,3)-ResultPdg(ctrl,2);% 长指令持续的时间长度
                if ResultPdg(ctrl,4)>300
                    ResultPdg(ctrl,6)=1;% 储能以2C功率连续出力在300s以上
                else
                    ResultPdg(ctrl,6)=0;
                end
                T=0;
                ctrl=ctrl+1;% 下一条长指令
                T_time=0;
            end
        end
    end
end
ctrl=1;
DetP=0;
DetN=0;
Agcstart=Agc(1);
% DetN=zeros(10,1);
% DetP=zeros(N,1);
for i=1:length(Agc)
    if (Agc(i) > detAGC+Agcstart) ||  (Agc(i) < Agcstart-detAGC)
        ctrl=ctrl+1;
        DetN(ctrl)=0;
        Agcstart=Agc(i);
    else
        DetP(i)=Agc(i)-P(i);
        DetN(ctrl)=DetP(i)*1+DetN(ctrl);% MW/s
    end
end
DetP=DetP';
DetN=DetN';
Pbat=Pall-P;
Numcycz = sum(Pbat(Pbat>0))*1/3600/9;
Numcycf = sum(Pbat(Pbat<0))*1/3600/9;


Q=0;
for i=1:288
    m=sum(Agc((i-1)*300+1:i*300));
    n=sum(Pall((i-1)*300+1:i*300));
    if n<m*0.98 || n>m*1.02
        Q=Q+1;
        i;
    end
end
Q