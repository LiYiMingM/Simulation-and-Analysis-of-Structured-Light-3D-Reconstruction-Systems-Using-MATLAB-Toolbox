clc
clear all
clf

geometry_number=2700;

root='G:\2023博士课题\dataset\仿真数据集\seperate_5.25_2700_270\train\left\';%根目录
root1='grating\';%根目录

%%%%%%此处是叠加相移光栅中提取低频条纹
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path_si_low=[root 'si_low_fre_grating\'];%保存低频光栅的路径%叠加相移光栅的话是si_low_fre_grating
% %%%%%%读低频叠加相移光栅
% for k=1:geometry_number
%         path_si=[root root1 sprintf('%06d', k) '\'];
%         si_low=load([path_si '1_6.mat'],'si');
%         si=si_low.si;
% %         figure
% %         imshow(si,[])
%         %%%保存低频叠加相移光栅
%          if ~exist(path_si_low,'dir')
%                 mkdir(path_si_low)
%          end
%         path=[path_si_low,num2str(k,'%06d'),'.mat'];%保存最后一步高频光栅,高频光栅的路径
%         save(path, 'si');
% 
% end


%%%%%%%%%%%%%%%%%此处是单侧光栅从光栅中提取低频条纹
path_si_low=[root 'low_fre_grating\'];%保存低频光栅的路径%叠加相移光栅的话是si_low_fre_grating
%%%%%%读低频叠加相移光栅
for k=1:geometry_number
        path_si=[root root1 sprintf('%06d', k) '\'];
        si_low=load([path_si '1_6.mat'],'grating_left');
        grating_left=si_low.grating_left;
%         figure
%         imshow(si,[])
        %%%保存低频叠加相移光栅
         if ~exist(path_si_low,'dir')
                mkdir(path_si_low)
         end
        path=[path_si_low,num2str(k,'%06d'),'.mat'];%保存最后一步高频光栅,高频光栅的路径
        save(path, 'grating_left');

end

disp('生成结束')
