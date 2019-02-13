function [BatPower,status] = ControlMethod1(Agc,Pdg,BatSoc)
% function [BatPower,status] = BatAgcMethod(Agc,Pdg,BatSoc,Verbose)
%
% ������ּ��ʵ�ִ���AGC�㷨��������AGC������ֵ�ͻ��鹦������财���ܳ����������
% ���룺
%	Agc��   ��������ʾAGC������ֵ���ɵ��ȸ�������λ��MW��
%	Pdg��   ��������ʾ�������ʵ�⹦��ֵ����λ��MW��
%   BatSoc��     ��������ؿ���������0��100����λ��%��
%   Verbose��    ��������ʾ�澯��ʾ�ȼ���0-9��9��ʾ���澯��0����ʾ�澯��
% �����
%	BatPower��	��������ʾ�����ܹ��ʣ���λ��MW���ŵ�Ϊ����
%	status��     ��������ʾ�����ķ���״̬��>=0��ʾ������<0��ʾ�쳣��
% �汾������޸��� 2016-09-13
% 2016-09-13 HB
% 2019-01-17 wby
% 1. ��������������������������Ƴ���׫дע�͡�
% ȫ�ֱ�������
global T             % ȫ��ʱ��
global PdgStart      % ��ʼ����
global AgcStart      % ��ʼAgc
global LastAgc       % ��һ�ε��ڵ�Agc
global Tstart        % ��ʼʱ��
global Pdg_adj_start % ����ĳ����������ʼ����
global fang          % ���ڷ���
global lastPdg       % ��һʱ�̵Ĵ��ܳ���
global T_fantiao     % ����ʱ��
global T_butiao      % ����ʱ��
global T_huantiao    % ����ʱ��
global lastPall      % ��һʱ�̵����ϳ���
global T_record      % ������¼ ��ʼʱ��
global Pdg_record    % ���鹦�ʣ���ʷ1Сʱ��¼
global flag          % ���鲻����¼��־
global SigFM         % һ�ε�Ƶ�ź�
global Agc_adj       % ʵ����Ӧָ��
global Pall_adj_start % �������ʱ����ʼ���Ϲ���
global PallStart     % ������ʼʱ�̵����Ϲ���

T01=5;% ��С����Ӧ����ʱ��
T02=3;% ��Ӧ�����ֶο��ԣ���������ٽ緶Χ���ϣ���ʹ���ȴ����ж�
ParaSOC=[50 15 85 10];
% ParaSoc=[����ֵ ���� ���� �ͻ���С]

deadK3=0.01;
Prgen=600;
Pmax=18;
Cdead=0.01;
Cdead_Res=0.005;
Flag_adj=1; % �Ƿ���Ӧ�����ı�־λ,1��Ӧ,0����Ӧ
Portion_adj=1; % ��Ӧ����,1->100%
% if T==1
%    ���ں�����
%     lastPdg=0;
%     lastPall=0;
% end
if abs(Agc-AgcStart)>1 % Cdead*Prgen % �̶���С:1 % �����С:Cdead*Prgen
    % ������һ��ָ���Ҫ���³�ʼ״̬
    PdgStart=Pdg;% ��ʼ�Ļ��鹦��
    LastAgc=AgcStart;
    PallStart=lastPall;
    AgcStart=Agc;% ��ʼ��AGC����
    Tstart=T;% ��ʼ�Ļ��鹦��
%     T_fantiao=0;
%     T_butiao=0;
    if PallStart<AgcStart
        fang=1;% ������
    else
        fang=-1;
    end
    Agc_adj=(AgcStart-PdgStart)*Portion_adj+PdgStart;
end
% DetP=AgcStart-Pdg;% �������
Vresp_ideal=PallStart*deadK3/T02;% �������Ӧ�����ٶ�
Vn=0.015*Prgen;% ��׼�����ٶ�MW/min
Vadj_ideal=5*Vn/60;% ����ĵ����ٶ�,MW/s
Socflag=0;
Ts=T-Tstart;
lastPbat=(lastPall-lastPdg);
if Ts <= T01
    if Ts<=T02
        % �ڶ�ʱ�Ĳ����£��ݲ����ǻ��鷴��������
        if Pdg-PallStart < fang*PallStart*deadK3 % ��2MW��  %��������1%�㣺PdgStart*0.01
            %�����鹦��δ�ﵽ��Ӧ������Χ��
            Pall_resp_ideal=PallStart+fang*Vresp_ideal*Ts;% ��������Ϲ���
            if fang>0 && Pdg>Pall_resp_ideal
                Pall_resp_ideal=Pdg;
            end
            if fang<0 && Pdg<Pall_resp_ideal
                Pall_resp_ideal=Pdg;
            end
            BatPower=Pall_resp_ideal-Pdg;% ���ܳ���
            if BatSoc<ParaSOC(2)
                BatPower=min(BatPower,0);% �����ޣ�ֻ�ܳ䲻�ܷ�
                Socflag=1;% Socά����־
            end
            if BatSoc>ParaSOC(3)
                BatPower=max(BatPower,0);% �����ޣ�ֻ�ܷŲ��ܳ�
                Socflag=1;% Socά����־
            end
        else
            Pall_resp_ideal=PallStart+(PallStart*deadK3+0.5)*fang;% ��������Ϲ���
            BatPower=Pall_resp_ideal-Pdg;% ���ܳ���
            if BatSoc<ParaSOC(2)
                BatPower=min(BatPower,0);% �����ޣ�ֻ�ܳ䲻�ܷ�
                Socflag=1;% Socά����־
            end
            if BatSoc>ParaSOC(3)
                BatPower=max(BatPower,0);% �����ޣ�ֻ�ܷŲ��ܳ�
                Socflag=1;% Socά����־
            end
        end
    else
        Pall_resp_ideal=PallStart+(PallStart*deadK3+0.5*(Ts-T02))*fang;% ��������Ϲ���
        BatPower=Pall_resp_ideal-Pdg;% ���ܳ���
        if BatSoc<ParaSOC(2)
            BatPower=min(BatPower,0);% �����ޣ�ֻ�ܳ䲻�ܷ�
            Socflag=1;% Socά����־
        end
        if BatSoc>ParaSOC(3)
            BatPower=max(BatPower,0);% �����ޣ�ֻ�ܷŲ��ܳ�
            Socflag=1;% Socά����־
        end
    end
else
    if Flag_adj==1
        if Ts == T01+1
            Pdg_adj_start=Pdg;
            Pall_adj_start=lastPall;
        end
        if lastPall>Agc_adj+fang*Cdead_Res*Prgen && lastPall<Agc_adj-fang*Cdead_Res*Prgen
            Pall_adj_ideal=Agc_adj;
        else
            Pall_adj_ideal=Pall_adj_start+fang*Vadj_ideal*(Ts-T01);
            Pall_adj_ideal=fang*min(abs(PallStart-Pall_adj_ideal),abs(PallStart-Agc))+PallStart;
        end
        BatPower=Pall_adj_ideal-Pdg;% ���ܳ���
        if BatSoc<ParaSOC(2)
            BatPower=min(BatPower,0);% �����ޣ�ֻ�ܳ䲻�ܷ�
            Socflag=1;% Socά����־
        end
        if BatSoc>ParaSOC(3)
            BatPower=max(BatPower,0);% �����ޣ�ֻ�ܷŲ��ܳ�
            Socflag=1;% Socά����־
        end
        if BatSoc<ParaSOC(3) && BatSoc>ParaSOC(2)
            % ����SOC����
            if abs(Agc_adj-Pdg)<Cdead_Res*Prgen
%                 % �ڲ�����Ӧ����ʱ�����鵽���趨ֵ�󣬴����˳�
%                 BatPower=0;
                % ������������������������SOC����ά��
                BatPower=Agc_adj-Pdg;
                if BatSoc<ParaSOC(1)-ParaSOC(4)
                    BatPower=min(BatPower,BatPower/2);% ��Ŀ���������£��������ٷ�
                    Socflag=1;% Socά����־
                end
                if BatSoc>ParaSOC(1)+ParaSOC(4)
                    BatPower=max(BatPower,BatPower/2);% ��Ŀ���������ϣ��������ٳ�
                    Socflag=1;% Socά����־
                end
            else
                % ����δ��������
                if abs(lastPall-Agc_adj)<Cdead_Res*Prgen
                    % ����δ������Ӧ�趨ֵ������һʱ�����ϵ����趨ֵ
%                     Pall_adj_ideal=Agc_adj-fang*Cdead_Res*Prgen;% �ﵽ�½�
                    Pall_adj_ideal=Agc_adj;% �趨Ϊ����ֵ
                    BatPower=Pall_adj_ideal-Pdg;% ���ܳ���
                end
                % ����δ������Ӧ�趨ֵ������Ҳδ�����趨ֵ
                % �������δ�ﵽ��������������Ҳδ��������
%                 if (Pdg-lastPdg)*fang<0
%                     % ���鷴��
%                     T_fantiao=T_fantiao+1;
%                     if T_fantiao>30
%                         % 30s����
%                         BatPower=0.8*(lastPall-lastPdg);% ��ÿ��20%�˳�
%                     end
%                 else
%                     if T_fantiao>0
%                         T_fantiao=T_fantiao-1;
%                     end
%                     if T_record==0
%                         if abs(Pdg-lastPdg)<Prgen*0.001
%                             T_record=T;
%                             T_butiao=T_butiao+1;
%                         end
%                     else
%                         if abs(Pdg-lastPdg)*fang<3.5 %0.003*Prgen/60
%                             % (��׼�������ʵ�20%)
%                             % ���鲻��������£���������
%                             T_butiao=T_butiao+1;
%                             a=find(Pdg_record~=0,1,'last');
%                             if T_butiao>200 && abs(Pdg_record(a)-Pdg_record(a-100))<3.75
%                                 % ʼĩ״̬��2MW�ڶ����ǻ�������
%                                 BatPower=0.95*(lastPall-lastPdg);% ��ÿ��5%�˳�
%                                 flag=1;
%                             end
% %                             if (Pdg-lastPdg)*fang<0.03*Prgen/60 %(��׼�������ʵ�20%)
% %                                 % ���鲻��������£���������
% %                                 T_butiao=T_butiao+1;
% %                                 if T_butiao>100
% %                                     BatPower=min(BatPower,Pmax);
% %                                     BatPower=max(BatPower,-Pmax);
% %                                     BatPower=0.98*BatPower;
% %                                     BatPower=0.98.^(T_butiao-100)*BatPower;
% %                                     BatPower=0.9*(lastPall-lastPdg);% ��ÿ��10%�����˳�
% %                                 end
%                         else
%                             T_record=0;
%                             T_butiao=0;
%                             flag=0;
%                         end
%                         if (Pdg-lastPdg)*fang<0.015*Prgen/60*0.75 && flag==0
%                             % ���鲻��������£�����������׼�����������75%��
%                             T_huantiao=T_huantiao+1;
%                             if T_huantiao>60
% %                                 BatPower=min(BatPower,Pmax);
% %                                 BatPower=max(BatPower,-Pmax);
%                                 BatPower=0.95*(lastPall-lastPdg);
% %                                 BatPower=0.95.^(T_huantiao-60)*BatPower;
%                             end
%                         else
% %                             T_record=0;
%                             T_huantiao=0;
% %                             T_butiao=0;
%                         end
%                     end
%                 end
            end
        end
    else
        BatPower=0.9*lastPbat;% ����Ӧ�󣬻����˳�
        if abs(BatPower)<0.1
            BatPower=0;
        end
    end
end
if T==4489
    ha=1;
end
if (Pdg-lastPdg)*fang<0
    % ���鷴��
    T_fantiao=T_fantiao+1;
    if T_fantiao>30
        % 30s����
        BatPower=0.8*lastPbat;% ��ÿ��20%�˳�
    end
else
    if T_fantiao>0
        T_fantiao=T_fantiao-1;
    end
    if T_record==0
        if abs(Pdg-lastPdg)<Prgen*0.001
            T_record=T;
            T_butiao=T_butiao+1;
        end
    else
        if abs(Pdg-lastPdg)*fang<3.5 %0.003*Prgen/60
            % (��׼�������ʵ�20%)
            % ���鲻��������£���������
            T_butiao=T_butiao+1;
            a=find(Pdg_record~=0,1,'last');
            if T_butiao>100 && abs(Pdg_record(a)-Pdg_record(a-50))<3.75
                % ʼĩ״̬��2MW�ڶ����ǻ�������
                BatPower=0.95*lastPbat;% ��ÿ��5%�˳�
                flag=1;
            end
%             if (Pdg-lastPdg)*fang<0.03*Prgen/60 %(��׼�������ʵ�20%)
%                 % ���鲻��������£���������
%                 T_butiao=T_butiao+1;
%                 if T_butiao>100
%                     BatPower=min(BatPower,Pmax);
%                     BatPower=max(BatPower,-Pmax);
%                     BatPower=0.98*BatPower;
%                     BatPower=0.98.^(T_butiao-100)*BatPower;
%                     BatPower=0.9*(lastPall-lastPdg);% ��ÿ��10%�����˳�
%                 end
        else
            T_record=0;
            T_butiao=0;
            flag=0;
        end
        if (Pdg-lastPdg)*fang<0.015*Prgen/60*0.75 && flag==0
            % ���鲻��������£�����������׼�����������75%��
            T_huantiao=T_huantiao+1;
            if T_huantiao>60
                % BatPower=min(BatPower,Pmax);
                % BatPower=max(BatPower,-Pmax);
                BatPower=0.95*lastPbat;
                % BatPower=0.95.^(T_huantiao-60)*BatPower;
            end
        else
            % T_record=0;
            T_huantiao=0;
            % T_butiao=0;
        end
    end
end
% if Ts>=0 && Ts<15
%     b=lastPall-lastPdg;
%     b=0.85*b;
%     if b>0 && b>BatPower
%         BatPower=b;
%     end
%     if b<0 && b<BatPower
%         BatPower=b;
%     end
% else
%     if abs(BatPower-(lastPall-lastPdg))>2.7
%         BatPower=(lastPall-lastPdg)+2.7*BatPower/abs(BatPower);
%     end
% end
%%% һ�ε�Ƶ�źŲ���Ӧ %%%
if SigFM==1
    BatPower=lastPbat;
end

BatPower=min(BatPower,lastPbat+2);
BatPower=max(BatPower,lastPbat-2);
status=Socflag;
lastPdg=Pdg;% ��һ���ӵĻ������
lastPall=BatPower+Pdg;% ��һ���е����ϳ���