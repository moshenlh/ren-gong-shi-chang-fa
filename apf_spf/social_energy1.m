function social_energy = social_energy1(CurrentPosx,CurrentPosy,Theta,robot_v,Obs_Group_posx,Obs_Group_posy,obs_dir,obs_v,num_Group,sigma_w,sigma_d,beta)
%输入：机器人当前x位置，y位置,机器人当前方向，速度(大小)，行人x位置,y位置，运动方向，运动速度大小,行人个数，三个参数
%输出：行人产生的社交敏感斥能，能量
delta_p=cell(num_Group,1);%元胞数组
for i=1:num_Group
    delta_p{i,1}=[Obs_Group_posx(i)-CurrentPosx,Obs_Group_posy(i)-CurrentPosy];
    dis_ob(i)=sqrt((CurrentPosx-Obs_Group_posx(i))^2+(CurrentPosy-Obs_Group_posy(i))^2);%距离,也就是delta_p的模长 
    if(dis_ob(i)<=0.0001)
            dis_ob(i)=dis_ob(i)+0.0001;   %分母不能为0
    end 
    %计算w
    E_dir(i)=-dot(delta_p{i,1},[obs_v(i)*cos(obs_dir(i)),obs_v(i)*sin(obs_dir(i))])/(obs_v(i)*dis_ob(i));
    w(i)=exp((-dis_ob(i)^2)/(2*sigma_w^2))*(0.5*(1+E_dir(i)))^beta;     %w是数值
    %计算d
    robot_v_x=robot_v*cos(Theta);
    robot_v_y=robot_v*sin(Theta);
    obs_v_x(i)=obs_v(i)*cos(obs_dir(i));
    obs_v_y(i)=obs_v(i)*sin(obs_dir(i));
    delta_v=[robot_v_x-obs_v_x(i),robot_v_y-obs_v_y(i)];
    Norm_2(i)=(norm(delta_v))^2;
    d(i)=norm(delta_p{i,1}-dot(delta_p{i,1},delta_v)*delta_v/Norm_2(i));%d是向量的模长
    %合成能量
    Energy(i)=w(i)*exp(-d(i)^2/(2*sigma_d^2));
end
social_energy=sum(Energy);%求和

end

