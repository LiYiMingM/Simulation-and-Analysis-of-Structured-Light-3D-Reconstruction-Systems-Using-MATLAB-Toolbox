clc
clear all
clf

geometry_number=13500;

wide_left=[226:1:226];%%%%先让交叉点固定
 
wide_right=[30:1:30];

root='G:\2023博士课题\dataset\6.13_9_shape_13500_2700\train\';%根目录
root1='left\geometry\';%根目录
path_si_geometry_left=[root 'left\geometry_left\'];
path_si_geometry_right=[root 'right\geometry_right\'];
 if ~exist(path_si_geometry_left,'dir')
        mkdir(path_si_geometry_left)
 end
 if ~exist(path_si_geometry_right,'dir')
        mkdir(path_si_geometry_right)
 end
%%%%%%读低频叠加相移光栅
for k=1:geometry_number
        path_si=[root root1 sprintf('%06d', k) '.mat' ];
        si_left=load(path_si,'geometry');
        si_left1=si_left.geometry;
        geometry1=si_left1;
%%%%%%%%%%%保存左侧光栅照射范围内的物体几何图       

        si_left1(:,227:256)=0;
        geometry=si_left1;

        path=[path_si_geometry_left,num2str(k,'%06d'),'.mat'];%保存最后一步高频光栅,高频光栅的路径
        save(path, 'geometry');
%%%%%%%%%%%保存右侧光栅照射范围内的物体几何图

        geometry1(:,1:29)=0;
        geometry=geometry1;

        path=[path_si_geometry_right,num2str(k,'%06d'),'.mat'];%保存最后一步高频光栅,高频光栅的路径
        save(path, 'geometry');


end
disp("结束")