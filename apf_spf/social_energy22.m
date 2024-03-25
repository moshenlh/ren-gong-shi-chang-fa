function Energy = social_energy22(robot_posx,robot_posy,obs_posx,obs_posy,obs_dir,obs_v)
%论文《Mobile Robot Navigation for Human-Robot Social Interaction》

num_Group=length(obs_posx(:,1));%行人个数

%参数
%改变sigma的值可以让蛋形缩小
sigma_f=0.2+num_Group*0.2;%这个调大了会超出作用范围，导致边缘被截断产生毛刺，因此warningdis最好是跟着这个参数变化的？
%群体内个数越多，蛋形应该越大
sigma_s=0.4;      %sigma_f/sigma_s越大，椭圆越长
gamma=2;
A=1;
warningdis=2+num_Group*0.5;

Energy=0;%赋初值
X=robot_posx;
Y=robot_posy;

for k=1:num_Group
    dis_ob=sqrt((X-obs_posx(k))^2+(Y-obs_posy(k))^2);%距离,数
    theta=atan2(Y-obs_posy(k),X-obs_posx(k));%角度,数,atan2得到弧度制
    
    %计算front
    
    if(cos(theta-obs_dir(k))>=0)
        front=(dis_ob*cos(theta-obs_dir(k)))^2/(2*sigma_f^2);%如果机器人在行人背后，则身前能量给小一些
    else
        front=(dis_ob*cos(theta-obs_dir(k)))^2/(2*(sigma_f/(1+gamma*obs_v(k)))^2);%大一些
    end
    if(dis_ob>=warningdis)
        front=100000;%不是给0，而是给一个很大的值，因为下面的势能是取负数计算指数函数
    end
    
    
    %计算side
    side=(dis_ob*sin(theta-obs_dir(k)))^2/(2*sigma_s^2);
    if(dis_ob>=warningdis)
        side=100000;
    end
    
    
    %合成计算势能
    Energy=Energy+A*exp(-(side+front));
    % Energy=coefficient.*exp(-dis_ob);      %纯距离高斯函数
end
end

