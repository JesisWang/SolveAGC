num=size(AGC);
x = linspace(datenum('2018-08-13 07:00:00'),datenum('2018-08-14 06:59:59'),num(1));   %����360��ʱ��
y = AGC;    % ����ʱ������y
plot(x,y);            % ��ͼ
% ���X������
N = 12;   %����������ʾN���̶�
% x����ʾ�̶ȵ�ʱ�����䣬�Լ�������N
date_point = linspace(datenum('2018-08-13 06:00:00'),datenum('2018-08-14 08:00:00'),14);
set(gca,'xtick',date_point );  % x���m��С�̶�
date_point_str = datestr(date_point,'yyyy-mm-dd HH:MM:SS');   %X��̶��϶�Ӧ���ַ�
set(gca,'xticklabel',date_point_str)   %��ʾС�̶ȵ�ֵ
set(gca,'XTickLabelRotation',30)