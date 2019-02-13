function [BatPower,status] = BatAgcMethod2(AgcLimit,GenPower,Pall,BatSoc,Verbose)
% function [BatPower,status] = BatAgcMethod(AgcLimit,GenPower,BatSoc,Verbose)
%
% ������ּ��ʵ�ִ���AGC�㷨��������AGC������ֵ�ͻ��鹦������财���ܳ����������
% ���룺
%	AgcLimit��   ��������ʾAGC������ֵ���ɵ��ȸ�������λ��MW��
%	GenPower��   ��������ʾ�������ʵ�⹦��ֵ����λ��MW��
%   BatSoc��     ��������ؿ���������0��100����λ��%��
%   Verbose��    ��������ʾ�澯��ʾ�ȼ���0-9��9��ʾ���澯��0����ʾ�澯��
% �����
%	BatPower��	��������ʾ�����ܹ��ʣ���λ��MW���ŵ�Ϊ����
%	status��     ��������ʾ�����ķ���״̬��>=0��ʾ������<0��ʾ�쳣��
% �汾������޸��� 2016-09-13
% 2016-09-13 HB
% 1. ��������������������������Ƴ���׫дע�͡�

    % ȫ�ֱ�������
    global Tline;   % �����ã����ڼ���ʱ�䡣
    global AgcStart;    % ��ʼʱ�䡣
    global GenPower0;   % ������ʼ������
    global LastAgc;     % ��������ʾ��һ��AGC������ֵ���ɵ��ȸ�������λMW��
    global LastPbat;	% ��������ʾ��һ�δ���ָ��ֵ�����㷨��ã���λMW��
    global Para;        % ������1��14��t01\t12\Pmax\Pmin\Phold\SocTarget\SocZone1\SocZone2\SocMax\SocMin\Erate\Prgen\Vgen\DeadZone���㷨������
    global LastAgcLimit;
    global SOC0;
    global Pall0;
    global SocFlag;     % �ѿ�ʼSOCά��
    global FlagAGC;     % �������ﵽAGCָ��
    global VgP0;
    global Vg;
    global Qiflag;
    
    status = -1;    % ��ɳ�ʼ����״̬Ϊ-1

    % ������
    if (isempty(Verbose)||isnan(Verbose))
        Verbose = 0;
    end
    if (isempty(AgcLimit)||isempty(GenPower)||isempty(BatSoc)||isempty(Para)) || ...
       (isnan(AgcLimit)||isnan(GenPower)||isnan(BatSoc)||(sum(isnan(Para))>0))
        % ������ڿ������NAN��״̬Ϊ-2
        status = -2;
        WarnLevel = 1;
        if WarnLevel < Verbose
            fprintf('Input data can not be empty or NaN!');
        end
        return;
    elseif (length(Para) ~= 14)
        % �������ݸ�ʽ������Ҫ��״̬Ϊ-3
        status = -3;
        WarnLevel = 1;
        if WarnLevel < Verbose
            fprintf('Para data is not correct format!');
        end
        return;
    elseif AgcLimit <= 0
        % AGC��ֵС�ڵ���0��״̬Ϊ0
        status = 0;
        WarnLevel = 3;
        if WarnLevel < Verbose
            fprintf('AGC limit is 0!');
        end
        return;
    end
    
    Erate = Para(11);
    Prgen = Para(12);
    Vgen = Para(13);   % MW/min
    Cdead = Para(14);
    % AGC�㷨
    if (AgcLimit > LastAgc+0.01*Prgen) ||  (AgcLimit < LastAgc-0.01*Prgen)    % AGCָ��仯
    %if (AgcLimit > LastAgc+5) ||  (AgcLimit < LastAgc-5)    % AGCָ��仯
%         if rand>0.2 %AgcLimit - GenPower>25 
%             Qiflag=1;
%         else
%             Qiflag=0;
%         end
        FlagAGC = 0;
        SocFlag = 0;
        AgcStart = Tline;
        GenPower0 = GenPower;
        LastAgcLimit = LastAgc;
        LastAgc = AgcLimit;
        SOC0 = BatSoc;
        Pall0 = Pall;
        VgP0 = GenPower;
    end
    DetP = AgcLimit - GenPower;     % �������������
    DetT = Tline-AgcStart;
    if (DetT<=Para(1))                  % AGCָ���ʼt01�ڣ��޹����������֤K3ָ��
        BatPower = DetP+(Pall0-AgcLimit)*0.95;% ����õ�K3ָ��
        if (AgcLimit<GenPower)
            Vg = -2; % �����������
        else
            Vg = 2;
        end
        status = 1;
    else
        if (DetT==Para(1)+1)
            Pall0 = Pall;
        elseif (mod(DetT,60)==0) % ÿ60s����һ�λ����ٶ�
            Vg = GenPower-VgP0;
            VgP0 = GenPower;
        end
        
        if (BatSoc<=Para(10))                                % SOC������ά����SocMin
            Pall0 = Pall;
            BatPower = min(DetP,-Para(5)/2); % ���ά������DetP<0ָ������Ҫ�����������ԱȽ�Phold˭��С�����������Ҫ�ŵ磬��ô���ţ�����(���Ǳ���С���ʳ����)
            SocFlag = 1;
        elseif (BatSoc>Para(9))                             % SOC������ά����SocMax     
            Pall0 = Pall;
            BatPower = max(DetP,Para(5)/2); % �ŵ�ά��
            SocFlag = 1;
        elseif (BatSoc>Para(10)+Para(8))&&(BatSoc<=Para(9)-Para(8))   % SOC�ڷ�Χ�ڣ��ɽ���AGC���ڣ�SocMin+SocDead��SocMax-SocDead��>�ͻ�Ч��
            SocFlag = 0;
            if (abs(DetP)>Cdead*Prgen)          % ����AGC���ڹ��̣�����K1\K2
                if (FlagAGC>0)  % �������������
                    if (BatSoc>Para(6)+Para(7))
                        BatPower = max(DetP,Para(5)/2);
                    elseif (BatSoc<Para(6)-Para(7))
                        BatPower = min(DetP,-Para(5)/2);
                    else
                        BatPower = 0;
                    end
                    SocFlag=1;
                elseif (DetT>Para(2))   % ��ʱ�䲹��������������
                    if (AgcLimit>GenPower0)
                        if Vg<0 % ���鷴��
                            Pall0 = 0.95*Pall0+0.05*GenPower; % ��ÿ����5%�Ĵ��ܹ���
                        elseif Vg<Vgen/2 % ����ʱ�Ƚ���
                            Pall0 = 0.98*Pall0+0.02*GenPower; % ��ÿ����2%�Ĵ��ܹ���
                        end
                        BatPower = min(DetP,Pall0-GenPower);
                        BatPower = max(0,BatPower);
                    else
                        if Vg>0
                            Pall0 = 0.95*Pall0+0.05*GenPower;
                        elseif Vg>-Vgen/2
                            Pall0 = 0.98*Pall0+0.02*GenPower;
                        end
                        BatPower = max(DetP,Pall0-GenPower);
                        BatPower = min(0,BatPower);
                    end
                else                    % ��ʱ��������������
                    if (AgcLimit<GenPower0) % �������
                        if (((70-BatSoc)*3600*Erate/100)>(DetP/Vgen*60*Para(4)-0.5*Para(4)/Vgen*60*Para(4)))
                            Pall0 = 0.6*Pall0+0.4*AgcLimit;
                        end
                        BatPower = max(DetP,Pall0-GenPower);
                        BatPower = min(0,BatPower);
                    else
                        if (((BatSoc-30)*3600*Erate/100)>(DetP/Vgen*60*Para(3)-0.5*Para(3)/Vgen*60*Para(3)))
                            Pall0 = 0.6*Pall0+0.4*AgcLimit;
                        end
                        BatPower = min(DetP,Pall0-GenPower);
                        BatPower = max(0,BatPower);
                    end
                end
            else                                % ���������ϣ�ά��SOC
                FlagAGC = 1;
                if (BatSoc>Para(6)+Para(7))
                    BatPower = max(DetP,Para(5));
                elseif (BatSoc<Para(6)-Para(7))
                    BatPower = min(DetP,-Para(5));
                else
                    BatPower = 0;
                end
                SocFlag=1;
            end
        else     % ��������Pb����ֱ��SOC�ص��ɽ���AGC�ķ�Χ
            if (abs(DetP)<=Cdead*Prgen)
                Pall0 = Pall;
                FlagAGC = 1;
                if (BatSoc>Para(6)+Para(7))
                    BatPower = max(DetP,Para(5));
                elseif (BatSoc<Para(6)-Para(7))
                    BatPower = min(DetP,-Para(5));
                else
                    BatPower = 0;
                end
                SocFlag=1;
            else
                if (SocFlag==1) % ֮ǰ�ڽ���SOCά�����򱣳�P���򲻱䣬��Ӧͬ��AGC����
                    Pall0 = Pall;
                    if LastPbat>=0
                        BatPower = max(DetP,Para(5));
                    else
                        BatPower = min(DetP,-Para(5));
                    end
                else            % ֮ǰ�ڽ���AGC���ڣ�������𵴣����ܻ�����Ϊ0��������ﵽAGCָ������SOCά��
                    if (Pall0-GenPower>0)
                        Pall0 = 0.95*Pall0+0.05*GenPower;
                        BatPower = max(0,Pall0-GenPower);
                    else
                        Pall0 = 0.95*Pall0+0.05*GenPower;
                        BatPower = min(0,Pall0-GenPower);
                    end
                    if (BatSoc>Para(6)+Para(7))
                        BatPower = max(DetP,BatPower);
                    elseif (BatSoc<Para(6)-Para(7))
                        BatPower = min(DetP,BatPower);
                    else
                        BatPower = 0;
                    end
                end
            end
        end
    end
    if (BatPower>Para(3))
        BatPower = Para(3);
    elseif  (BatPower<Para(4))
        BatPower = Para(4);
    end
    if abs(BatPower-LastPbat)<0.2
        BatPower = LastPbat;
    end
%     if Qiflag==0
%         BatPower=0;
%     end
end

