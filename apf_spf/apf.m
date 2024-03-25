%apf
%2023.10.09

clf;
close all;
clear;

%% 起始点位置
BegX = 1;                                   % 出发点位置
BegY = 8;
DesX = 12;                                    % 终点位置
DesY = 1;

% 文件导入障碍物
dataObstacle   = load('./apf_spf1009.txt');
ObX            = dataObstacle(:,2);
ObY            = dataObstacle(:,3);
Obtheta        = dataObstacle(:,4);             %障碍物运动方向
radius_warn    = dataObstacle(:,6);             % 障碍物影响范围半径
Ob_krep        = dataObstacle(:,7);             % 障碍物附加斥力系数
num_Ob         = length(ObX(:,1));              % 障碍物数量

% 参数
Katt = 1.2;
Krep = 1.95;
k_v  = 1.001;
ka   = 0.3;                                  %合成势能时引力势能的系数
kr   = 1.5;
steprate = 0.020;                            % 步长
steptime = 0.005;
max_interation = 3000;

%% 轨迹规划
current_x=BegX;                                                    % 初始化当前位置
current_y=BegY;
d_des=sqrt((current_x-DesX)^2 + (current_y-DesY)^2);               % 当前点和目标点的距离
num_iteration=0;
x_record=[];
y_record=[];

while(1)
    d_des=sqrt((current_x-DesX)^2 + (current_y-DesY)^2);               % 当前点和目标点的距离
    num_iteration=num_iteration+1;
    
    if(num_iteration>max_interation)
        fprintf('超时');
        break;
    end
    
    if(d_des<0.1)
        fprintf('到达');
        break;
    end
    
    % 引力
    [Fattx,Fatty] = Attractive(current_x,DesX,current_y,DesY,Katt);    % 引力计算
    
    % 斥力
    Frepx = zeros(1,num_Ob);
    Frepy = zeros(1,num_Ob);
    for i = 1:num_Ob
        [Frepx(1,i),Frepy(1,i)] = Repulsive(current_x,current_y,ObX(i),ObY(i),Krep + Ob_krep(i) ,radius_warn(i));      % 斥力计算
    end
    
    % 合力
    Fxsum = ka * Fattx + kr * sum(Frepx);
    Fysum = ka * Fatty + kr * sum(Frepy);
    v = k_v * sqrt(Fxsum^2 + Fysum^2);         %速度规划
    theta=atan2(Fysum,Fxsum);                  % 速度方向
%     current_x=current_x + steprate * Fxsum;  %更新位置
%     current_y=current_y + steprate * Fysum;
    current_x=current_x + steptime * v * cos(theta);  %更新位置
    current_y=current_y + steptime * v * sin(theta);
    
    x_record=[x_record;current_x];
    y_record=[y_record;current_y];
    
end

%% 画轨迹图
figure(1)

a = plot(BegX,BegY,'s','LineWidth',5,'color','g');hold on;
b = plot(DesX,DesY,'h','LineWidth',5','color',[0.8549,0.64706,0.12549]');hold on;
for i=1:num_Ob
    c = scatter(ObX(i),ObY(i),20*radius_warn(i),'k','filled');
    hold on;
end
d = plot(x_record,y_record,'b','LineWidth',2);hold on;
grid on;
xlabel("x / m");
ylabel("y / m");
xlim([1,12]);
ylim([-1,10]);
title('轨迹图');
% legend('起始点','终止点','障碍物','规划轨迹APF','location','best');
legend([a,b,c,d],{'StartPoint','EndPoint','Obstacles','APF Trajectory'},'location','best');

%% 画2D势场
figure(2)
x = 1:0.11:12;
y = 10:-0.11:-1;
lengthx=length(x);
lengthy=length(y);
[X,Y]=meshgrid(x,y);        %网格化

Uatt = 0.5 * Katt * ((X-DesX).^2 + (Y-DesY).^2);  %引力势能

% 斥力
Urepi = zeros(num_Ob,lengthx,lengthy);
for i = 1:num_Ob
    for j =1:lengthx
        for k =1:lengthy
%             Pobs = sqrt((x(j)-ObX(i))^2 + (y(k)-ObY(i))^2);  %Pobs是一个数
            Pobs = sqrt((X(j,k)-ObX(i))^2 + (Y(j,k)-ObY(i))^2);  %用矩阵不容易出错
            if(Pobs<=radius_warn(i))
                Urepi(i,j,k) = Urepi(i,j,k) + 0.5 * Krep * (1/Pobs - 1/radius_warn(i))^2; % 斥力势能
            else
                Urepi(i,j,k) = Urepi(i,j,k) + 0;
            end
        end
    end
end

%截断
random_limit = rand(num_Ob);
for i =1:num_Ob
    for j =1:lengthx
        for k =1:lengthy
%             if(Urepi(i,j,k)> 20 + 200 * random_limit(i))%不规则截断
%                 Urepi(i,j,k) = 20 + 200*random_limit(i);
                if(Urepi(i,j,k)> 7)
                Urepi(i,j,k) = 7;
            end
        end
    end
end

%所有障碍物斥力势能加起来
Urep = zeros(lengthx,lengthy);
for i =1:num_Ob
    for j =1:lengthx
        for k =1:lengthy
            Urep(j,k) = Urep(j,k)+Urepi(i,j,k);
        end
    end
end
U = ka * Uatt + kr * Urep;
% U = 1 * Urep;

p3 = plot(x_record,y_record,'b','LineWidth',2);hold on;
p1 = plot(BegX,BegY,'s','LineWidth',5,'color','g');hold on;
p2 = plot(DesX,DesY,'h','LineWidth',5,'color','[0.8549,0.64706,0.12549]');hold on;
grid on;
contour(X,Y,U);
% xlim([-1,5]);
% ylim([-2,4]);
xlabel("x / m");
ylabel("y / m");
legend([p1,p2,p3],{'StartPoint','EndPoint','APF Trajectory'},'location','best');
%% 画3D势场
figure(3)
mesh(X,Y,U);hold on;

limit = 0.1;  %需要选取合适的之，太小的话捕捉不到势能高度值，会为取到初始化值0
length_path=length(x_record);
U_3d=zeros(1,length_path);
U_3d_temp=zeros(1,length_path);

for i = 1:length_path
    for j =1:lengthx
        for k =1:lengthy
            if (abs(x_record(i)-X(j,k))<=limit && abs(y_record(i)-Y(j,k))<=limit)
                U_3d_temp(i) = U(j,k);
            end
            if U_3d(i) <= U_3d_temp(i)  %取接近值中最大的
                U_3d(i) = U_3d_temp(i);
            end
        end
    end
end
% p31 = plot3(x_record(1),y_record(1),U_3d(1),'go','LineWidth',6);hold on;%起点
% p32 = plot3(x_record(end),y_record(end),U_3d(end),'o','LineWidth',6,'color','[0.8549,0.64706,0.12549]');hold on;%终点
p33 = plot3(x_record,y_record,U_3d,'b-','LineWidth',2);hold on;
grid on
% xlim([-1,5]);
% ylim([-2,4]);
xlabel("x / m");
ylabel("y / m");
zlabel("U");
% legend([p31,p32,p33],{'起点','终点','APF轨迹'},'location','best');