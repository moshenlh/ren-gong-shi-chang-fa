function Energy = social_energy2(CurrentPos,obs_pos,obs_dir,obs_v_org,num_Group)
%输入：机器人位置，障碍物位置,速度大小，方向

%参数
sigma_f=0.6;%这个调大了会超出作用范围，导致边缘被截断产生毛刺，因此warningdis最好是跟着这个参数变化的？
sigma_s=0.4;      %sigma_f/sigma_s越大，椭圆越长
gamma=2;
A=1;
warningdis=2.3+num_Group*0.2;

obs_v=obs_v_org*100000;%放大速度，才会出现蛋形
obs_v=obs_v_org;
X=CurrentPos(1);
Y=CurrentPos(2);

Energy=0;

for i=1:num_Group
    dis_ob(i)=sqrt((X-obs_pos(i,1))^2+(Y-obs_pos(i,2))^2);%距离
    theta(i)=atan2(Y-obs_pos(i,2),X-obs_pos(i,1));%角度
    
    %计算front
    if(cos(theta(i)-obs_dir(i))>=0)
        front(i)=(dis_ob(i)*cos(theta(i)-obs_dir(i)))^2/(2*sigma_f^2);%如果机器人在行人背后，则身前能量给小一些
    else
        front(i)=(dis_ob(i)*cos(theta(i)-obs_dir(i)))^2/(2*(sigma_f/(1+gamma*obs_v(i)))^2);%大一些
    end
    if(dis_ob(i)>=warningdis)
        front(i)=100000;%不是给0，而是给一个很大的值，因为下面的势能是取负数计算指数函数
    end
    
    %计算side
    side(i)=(dis_ob(i)*sin(theta(i)-obs_dir(i)))^2/(2*sigma_s^2);
    if(dis_ob(i)>=warningdis)
        side(i)=100000;
    end
    
    %合成计算势能
    Energy=Energy+A*exp(-(side(i)+front(i)));
    % Energy=coefficient.*exp(-dis_ob);      %纯距离高斯函数
end

end

