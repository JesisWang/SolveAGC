num=size(AGC);
x = linspace(datenum('2018-08-13 07:00:00'),datenum('2018-08-14 06:59:59'),num(1));   %生成360个时间
y = AGC;    % 生成时间序列y
plot(x,y);            % 画图
% 设计X轴坐标
N = 12;   %坐标轴上显示N个刻度
% x轴显示刻度的时间区间，以及区间数N
date_point = linspace(datenum('2018-08-13 06:00:00'),datenum('2018-08-14 08:00:00'),14);
set(gca,'xtick',date_point );  % x轴分m个小刻度
date_point_str = datestr(date_point,'yyyy-mm-dd HH:MM:SS');   %X轴刻度上对应的字符
set(gca,'xticklabel',date_point_str)   %显示小刻度的值
set(gca,'XTickLabelRotation',30)