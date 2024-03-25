%spf
%2023.10.09

clf;
close all;
clear;
%% 起始点位置
BegX = 1;                                   % 出发点位置
BegY = 8;
DesX = 10;                                    % 终点位置
DesY = 1;

% 文件导入障碍物
dataObstacle   = load('./apf_spf1009.txt');
ObX            = dataObstacle(:,2);
ObY            = dataObstacle(:,3);
Obtheta        = dataObstacle(:,4);             % 障碍物运动方向   ？？？为啥只能读出来4*pi/4的第一个数字？？？  
Obtheta = Obtheta*pi/180;    %朝向数据改成角度制，读入时处理一下
Obv            = dataObstacle(:,5);             % 障碍物速度
radius_warn    = dataObstacle(:,6);             % 障碍物影响范围半径
Ob_krep        = dataObstacle(:,7);             % 障碍物附加斥力系数
num_Ob         = length(ObX(:,1));              % 障碍物数量

% 参数
Katt = 0.4;
Krep = 2.95;
k_v  = 0.3;
ka   = 0.3;                                  %合成势能时引力势能的系数
kr   = 1.5;
k_normal = 0.001;                            %合成斥力时普通斥力的系数
steprate = 0.030;                            % 步长
steptime = 0.06;
max_interation = 4000;

%社交敏感参数
sigma_w   =1.0;%越大越胖
sigma_d   =1.5;%越大越圆
beta      =2;%越小越胖，且越高

%动量系数
alpha          = 0.7;        %越大越尖，刺
beta           = 0.7;        %越大波动越大
%% 轨迹规划
current_x=BegX;                                                    % 初始化当前位置
current_y=BegY;
theta    =0;                                                       %初始化
v        =0;
d_des=sqrt((current_x-DesX)^2 + (current_y-DesY)^2);               % 当前点和目标点的距离
num_iteration=0;
x_record=[];
y_record=[];
Fx_record=[];
Fy_record=[];
Fxsum_record=[];
Fysum_record=[];
Fx =0;
Fy =0;

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
    
    % 一样的引力
    [Fattx,Fatty] = Attractive(current_x,DesX,current_y,DesY,Katt);    % 引力计算
    
    % 不一样的斥力
    % 使用函数1
%     energy   = social_energy1(current_x,current_y,theta,v,ObX,ObY,Obtheta,Obv,num_Ob,sigma_w,sigma_d,beta);   %计算社交敏感斥力势能
%     energy_x = social_energy1(current_x+0.0001,current_y,theta,v,ObX,ObY,Obtheta,Obv,num_Ob,sigma_w,sigma_d,beta);
%     energy_y = social_energy1(current_x,current_y+0.0001,theta,v,ObX,ObY,Obtheta,Obv,num_Ob,sigma_w,sigma_d,beta);
    %     %使用函数2
        energy     = social_energy22(current_x,current_y,ObX,ObY,Obtheta,Obv);
        energy_x   = social_energy22(current_x+0.0001,current_y,ObX,ObY,Obtheta,Obv);
        energy_y   = social_energy22(current_x,current_y+0.0001,ObX,ObY,Obtheta,Obv);
    Fx_spf   = -(energy_x-energy)/0.0001;     % 势能求梯度得到力  -----
    Fy_spf   = -(energy_y-energy)/0.0001;
   
    %传统斥力
    Frepx = zeros(1,num_Ob);
    Frepy = zeros(1,num_Ob);
    for i = 1:num_Ob
        [Frepx(1,i),Frepy(1,i)] = Repulsive(current_x,current_y,ObX(i),ObY(i),Krep + Ob_krep(i) ,radius_warn(i));      % 斥力计算
    end
    
    % 合力
    Fxsum = ka * Fattx + kr * (Fx_spf + k_normal*sum(Frepx));
    Fysum = ka * Fatty + kr * (Fy_spf + k_normal*sum(Frepy));
    
%     v = k_v * sqrt(Fxsum^2 + Fysum^2);         % 机器人速度规划
%     if(v>=0.3)
%         v=0.3;%速度限制
%     end
%     theta=atan2(Fysum,Fxsum);                  % 机器人速度方向
%     current_x=current_x + steprate * Fxsum;  % 更新机器人位置
%     current_y=current_y + steprate * Fysum;
% %     current_x=current_x + steptime * v * cos(theta);  % 更新机器人位置
% %     current_y=current_y + steptime * v * sin(theta);

%%%%%%%%%%%%加入动量势场法%%%%%%%%%%%%
    Fx = alpha*Fxsum + beta*Fx;
    Fy = alpha*Fysum + beta*Fy;

    v = k_v * sqrt(Fx^2 + Fy^2);         % 机器人速度规划
    if(v>=0.3)
        v=0.3;%速度限制
    end
    theta=atan2(Fy,Fx);                  % 机器人速度方向
    current_x=current_x + steprate * Fx;  % 更新机器人位置
    current_y=current_y + steprate * Fy;
%%%%%%%%%%%%加入动量势场法%%%%%%%%%%%%   

    x_record=[x_record;current_x];
    y_record=[y_record;current_y];
    
    Fx_record=[Fx_record;Fx];
    Fy_record=[Fy_record;Fy]; 
    
    Fxsum_record=[Fxsum_record;Fxsum];
    Fysum_record=[Fysum_record;Fysum];
end
figure(4);
plot(Fx_record,Fy_record,'r','LineWidth',2);hold on;
plot(Fxsum_record,Fysum_record,'b','LineWidth',2);hold on;
%% 画轨迹图
figure(1)

a = plot(BegX,BegY,'s','LineWidth',5,'color','g');hold on;
b = plot(DesX,DesY,'h','LineWidth',5','color',[0.8549,0.64706,0.12549]');hold on;
for i=1:num_Ob
    c = scatter(ObX(i),ObY(i),20*radius_warn(i),'k','filled');
    hold on;
end
d = plot(x_record,y_record,'r','LineWidth',2);hold on;
grid on;
xlabel("x / m");
ylabel("y / m");
xlim([1,12]);
ylim([-1,10]);
title('轨迹图');
% legend('起始点','终止点','障碍物','规划轨迹APF','location','best');
legend([a,b,c,d],{'StartPoint','EndPoint','Obstacles','SPF Trajectory'},'location','best');

%% 画2D势场
figure(2)
x = 1:0.11:12;
y = 10:-0.11:-1;
lengthx=length(x);
lengthy=length(y);
[X,Y]=meshgrid(x,y);        %网格化

Uatt = 0.5 * Katt * ((X-DesX).^2 + (Y-DesY).^2);  %引力势能

% 社交敏感斥力
theta_robot = 0;   %设置机器人朝向
v_robot     = 1;   %设置机器人速度
Urepi = zeros(num_Ob,lengthx,lengthy);
Pobs  = zeros(num_Ob,lengthx,lengthy);
for i = 1:num_Ob
    for j =1:lengthx
        for k =1:lengthy
            Pobs(i,j,k) = sqrt((X(j,k)-ObX(i))^2 + (Y(j,k)-ObY(i))^2);  %Pobs是一个数，用矩阵不容易出错
            if(Pobs(i,j,k)<=radius_warn(i))
                %使用函数1
                %Urepi(i,j,k) = Urepi(i,j,k) + social_energy1(X(j,k),Y(j,k),theta_robot,v_robot,ObX(i),ObY(i),Obtheta(i),Obv(i),1,sigma_w,sigma_d,beta);   % 斥力势能
                %使用函数2
                Urepi(i,j,k) = Urepi(i,j,k) + social_energy22(X(j,k),Y(j,k),ObX(i),ObY(i),Obtheta(i),Obv(i));   %这里用了X（j,k)Y（j,k)没用x(j),y(k)真的不想说什么。。。
            else
                Urepi(i,j,k) = Urepi(i,j,k) + 0;
            end
        end
    end
end

% 普通斥力
Urepi_normal = zeros(num_Ob,lengthx,lengthy);
for i = 1:num_Ob
    for j =1:lengthx
        for k =1:lengthy
%             Pobs = sqrt((x(j)-ObX(i))^2 + (y(k)-ObY(i))^2);  %Pobs是一个数
            Pobs = sqrt((X(j,k)-ObX(i))^2 + (Y(j,k)-ObY(i))^2);  %用矩阵不容易出错
            if(Pobs<=radius_warn(i))
                Urepi_normal(i,j,k) = Urepi_normal(i,j,k) + 0.5 * Krep * (1/Pobs - 1/radius_warn(i))^2; % 斥力势能
            else
                Urepi_normal(i,j,k) = Urepi_normal(i,j,k) + 0;
            end
        end
    end
end

%社交斥力势能加起来
Urep = zeros(lengthx,lengthy);
for i =1:num_Ob
    for j =1:lengthx
        for k =1:lengthy
            Urep(j,k) = Urep(j,k)+Urepi(i,j,k);
        end
    end
end

%普通斥力势能加起来
Urep_normal = zeros(lengthx,lengthy);
for i =1:num_Ob
    for j =1:lengthx
        for k =1:lengthy
            Urep_normal(j,k) = Urep_normal(j,k)+Urepi_normal(i,j,k);
        end
    end
end

Urep_all = Urep + k_normal*Urep_normal;

%截断
random_limit = rand(num_Ob);
for j =1:lengthx
    for k =1:lengthy
        if(Urep_all(j,k)>1.2)
            Urep_all(j,k) = 1.2;
        end
    end
end

U = 1 * Urep_all;
U = ka * Uatt + kr * Urep_all;

p3 = plot(x_record,y_record,'r','LineWidth',2);hold on;
p1 = plot(BegX,BegY,'s','LineWidth',5,'color','g');hold on;
p2 = plot(DesX,DesY,'h','LineWidth',5,'color','[0.8549,0.64706,0.12549]');hold on;
grid on;
contour(X,Y,U);
% xlim([-1,5]);
% ylim([-2,4]);
xlabel("x / m");
ylabel("y / m");
legend([p1,p2,p3],{'StartPoint','EndPoint','SPF Trajectory'},'location','best');
%% 画3D势场
figure(3)
mesh(X,Y,U);hold on;

limit = 0.1;
length_path=length(x_record);
U_3d=zeros(1,length_path);
U_3d_temp=zeros(1,length_path);

for i = 1:length_path  %路径上的每个点
    for j =1:lengthx
        for k =1:lengthy
            if (abs(x_record(i)-X(j,k))<=limit && abs(y_record(i)-Y(j,k))<=limit)
                U_3d_temp(i) = U(j,k);
            end
            if U_3d(i) <= U_3d_temp(i)
                U_3d(i) = U_3d_temp(i);
            end
        end
    end
end
p31 = plot3(x_record(1),y_record(1),U_3d(1),'go','LineWidth',6);hold on;%起点
p32 = plot3(x_record(end),y_record(end),U_3d(end),'o','LineWidth',6,'color','[0.8549,0.64706,0.12549]');hold on;%终点
p34 = plot3(x_record,y_record,U_3d,'r-','LineWidth',2);hold on;
grid on
% xlim([-1,5]);
% ylim([-2,4]);
xlabel("x / m");
ylabel("y / m");
zlabel("U");
%legend([p31,p32,p33,p34],{'StartPoint','EndPoint','APF Trajectory','MPF Trajectory'},'location','best');
legend([p31,p32,p34],{'StartPoint','EndPoint','SPF Trajectory'},'location','best');