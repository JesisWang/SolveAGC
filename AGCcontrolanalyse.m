% clear
load('SHDdata.mat')
load('SHDPCSdata1228deal')
data=SHDdata.data1228(:,4:6);
Agc=data(1:2:end,1);% AGCָ��
P=data(1:2:end,2);% �������
% Pbat=SDPCSdatadeal.all(:,37);
% Pall=Pbat+P;
Pall=data(1:2:end,3);% ���ϳ���
LineMax=length(Agc);% һ������ݵ��ʱ��
Result=zeros(10,12);% �����������
% Result=[1         2        3          4          5      6         7          8              9               10             11            12]
% Result=[AGC ������ʼʱ��ֵ t0 �����״ο�ʼ����ֵ t1 �����¼  ������¼  ����������¼  ˲��������  ����Ч������1,��0�����ڷ���  ���鱾��������]
ctrlNo=1;
Result(ctrlNo,1)=Agc(1);
Result(ctrlNo,2)=P(1);
Result(ctrlNo,3)=1;
detAGC=2;
Prate=600;
Erate=12;
T1=0;% ������ʱ
Cdead=0.02;
T2i=0;% ������ʼʱ�̣�Ϊ�˱�֤����һ��ʱ���ھ��ǻ���������
T2=0;% ������ʱ
Tlen=40;
if Result(ctrlNo,1)>Result(ctrlNo,2)
    Result(ctrlNo,11)=1;
else
    Result(ctrlNo,11)=-1;
end
for i=1:LineMax
    if (Agc(i) > detAGC+Result(ctrlNo,1)) ||  (Agc(i) < Result(ctrlNo,1)-detAGC)
        % ����һ��ָ������ϵ�ָ������
        if Result(ctrlNo,6)==0 && Result(ctrlNo,7)==0 
            % û�з����򲻵�,�Լ�����δ��������������£��������ĵ��������Ƿ���
            if (P(i)-Result(ctrlNo,4))/(i-Result(ctrlNo,5))*60<0.015*Prate
                % �������������С�ڱ�׼����
                if Result(ctrlNo,12)==0
                    % �һ���û�е������������Χ��
                    Result(ctrlNo,8)=1;
                end
            end
        end
%         if abs(Pall(i)-Result(ctrlNo,1))<Cdead
%             % ���ϳ���������ڷ�Χ�ڣ��ǵ��ڳɹ�һ��
%             Result(ctrlNo,10)=1;
%         end
        if Result(ctrlNo,5)==0 
            Result(ctrlNo,7)=1;
        end
        if i-Result(ctrlNo,3)>Tlen
            % ָ�����ʱ��С��40�룬ָ������
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
            % �жϻ��鿪ʼ��Ӧָ��ʱ��
            Result(ctrlNo,4)=P(i);
            Result(ctrlNo,5)=i;
        end
        if Result(ctrlNo,6)==0
            % ���鷴�����
            if i>2 && (P(i)-P(i-1))*(Result(ctrlNo,11))<0
                % �������������ָ����෴
                T1=T1+1;
                if T1>=10 && (P(i)-P(i-T1))*Result(ctrlNo,11)<-2
                    % ����10s���ж��Ƿ�������ʼĩ��ֵҪ����-2MW
                    Result(ctrlNo,6)=1;
                    T1=0;
                end
            end
        end
        if Result(ctrlNo,7)==0 && Result(ctrlNo,5)~=0
            % ���鲻����¼
            if abs(Result(ctrlNo,1)-P(i))>Cdead*Result(ctrlNo,1)
                % ����������ڵ��������ڲ���������
                if T2i==0
                    if i>2 && abs(P(i)-P(i-1))<0.2
                        % ��¼������ʼʱ�̣�����ʼ����
                        T2i=i;
                        T2=T2+1;
                    end
                else
                    if i>2 && abs(P(i)-P(i-1))<0.2 && i-T2i==T2
                        % ���鲨������0.2MW,������С���Ĳ����������������ƻ�
                        T2=T2+1;
                        if T2>=Tlen/2 && abs(P(i)-P(i-T2))<0.4
                            % �ﵽ20s������ʱ��������������ʼ��������0.4MW�ڣ����Ϊ����
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
            % 3s��
            if abs(Pall(i)-Result(ctrlNo,1))<1  % && abs(Pall(i)-P(i))>Erate/2
                % ���ܽ���������AGC�������ҹ��ʴ���0.5�������
                Result(ctrlNo,9)=1;
            end
        end
        if Result(ctrlNo,12)==0 && abs(Result(ctrlNo,1)-P(i))<Cdead*Result(ctrlNo,1)
            Result(ctrlNo,12)=1;
        end
        if abs(Pall(i)-Result(ctrlNo,1))<Cdead*Result(ctrlNo,1) && Result(ctrlNo,10)==0
            % ���ϳ���������ڷ�Χ�ڣ��ǵ��ڳɹ�һ��
            Result(ctrlNo,10)=1;
        end
    end
end
%% ����
N=length(Result(:,6));
M=0;
Op1=sum(Result(:,6))/N*100 % ���鷴������
Op2=sum(Result(:,7))/N*100 % ���鳤ʱ�䲻���ı���
Op3=sum(Result(:,8))/N*100 % ��������ٶ�С�ڱ�׼�����ٶȵı���
Op4=sum(Result(:,9))/N*100 % 3s�ڽ������������AGCָ����
S=0;
for i=1:N
    if Result(i,6)==1 || Result(i,7)==1 || Result(i,8)==1
        M=M+1;
        if Result(i,10)==1
            S=S+1;
        end
    end
end
S %����ǰ���ܵ��ڲ����£��ڻ�������������������£���Ч���ڣ��ɻ�Kpֵ���Ĵ���
Sop=S/M*100  % ��Ӧ����
MchuN=M/N*100 % ��Чָ��
%% ����ͳ��
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
        E_point5=E_point5+abs(Pdg(i))*1/3600;% ��0.5C��Ч
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
            E_1plus=E_1plus+Pdg(i)*1/3600;% ��1C��Ч
        else
            E_1minus=E_1minus+Pdg(i)*1/3600;
        end
    end
    if abs(Pdg(i))>9  && abs(Pdg(i))<13.5
        M_1point5=M_1point5+1;
        E_1point5=E_1point5+abs(Pdg(i))*1/3600;
        if Pdg(i)>0
            E_1point5plus=E_1point5plus+Pdg(i)*1/3600; % ��1.5C��Ч
        else
            E_1point5minus=E_1point5minus+Pdg(i)*1/3600;
        end
    end
    if abs(Pdg(i))>13.5  && abs(Pdg(i))<18
        M_2=M_2+1;
        E_2=E_2+abs(Pdg(i))*1/3600;
        if Pdg(i)>0
            E_2plus=E_2plus+Pdg(i)*1/3600;% ��2C��Ч
        else
            E_2minus=E_2minus+Pdg(i)*1/3600;
        end
    end
    if abs(Pdg(i))>13.5
        if T==0
            ResultPdg(ctrl,1)=Pdg(i);% ��ʼʱ�̵Ĺ���
            ResultPdg(ctrl,2)=i;% ���������ʼʱ��
            ResultPdg(ctrl,5)=abs(Pdg(i))/3600;% �ۼ�����
            T=T+1;
        else
            ResultPdg(ctrl,5)=ResultPdg(ctrl,5)+abs(Pdg(i))/3600;
        end
        T_time=0;
    else
        T_time=T_time+1;
        if T_time<60
            % ����60s���ڵ���2C�Ĺ��ʣ������ж�2C
%             ResultPdg(ctrl,5)=ResultPdg(ctrl,5)+Pdg(i)/3600;
        else
            % ����60s����2C�������ж�2C����ʱ��
            if T~=0
                ResultPdg(ctrl,3)=i-T_time;% �ٳ�60s�����㳤ָ�������ʱ��
                ResultPdg(ctrl,4)=ResultPdg(ctrl,3)-ResultPdg(ctrl,2);% ��ָ�������ʱ�䳤��
                if ResultPdg(ctrl,4)>300
                    ResultPdg(ctrl,6)=1;% ������2C��������������300s����
                else
                    ResultPdg(ctrl,6)=0;
                end
                T=0;
                ctrl=ctrl+1;% ��һ����ָ��
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