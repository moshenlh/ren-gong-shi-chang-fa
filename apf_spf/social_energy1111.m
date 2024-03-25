function Dy_Rep_Energy = social_energy1111(CurrentPos,Theta,robot_v,Obs_Group,num_Group,sigma_w,sigma_d,beta)
%输入：机器人当前位置（xy)，机器人当前方向，速度(标量)，Obs_Group(12行人位置xy，3运动方向，4运动速度),行人个数，三个参数
%输出：行人产生的斥能

for obCount=1:num_Group
    dis_Group(obCount)=sqrt((CurrentPos(1)-Obs_Group(obCount,1))^2+(CurrentPos(2)-Obs_Group(obCount,2))^2);%距离，也是delta_p的模长
    delta_p(obCount,1:2)=[CurrentPos(1)-Obs_Group(obCount,1),CurrentPos(2)-Obs_Group(obCount,2)];%是向量[delta_x,delta_y]
    E_dir=-(delta_p(obCount,1)/dis_Group(obCount)*cos(Theta)+delta_p(obCount,2)/dis_Group(obCount)*sin(Theta));%两个向量相乘，等于对应项分别相乘再相加
    w(obCount)=exp(-dis_Group(obCount)^2/(2*sigma_w^2))*(0.5*(1+E_dir))^beta;
    delta_v_2(obCount)=(robot_v*cos(Theta)-Obs_Group(obCount,3)*cos(Obs_Group(obCount,4)))^2+(robot_v*sin(Theta)-Obs_Group(obCount,3)*sin(Obs_Group(obCount,4)))^2;
    d(obCount)=delta_p(obCount)-delta_p(obCount)*delta_v_2(obCount)/delta_v_2(obCount);
    Energy(obCount)=w(obCount)*exp(-d(obCount)^2/(2*sigma_d^2));
end
Dy_Rep_Energy=sum(Energy);%求和
end

