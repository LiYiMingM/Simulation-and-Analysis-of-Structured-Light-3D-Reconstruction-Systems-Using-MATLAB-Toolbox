%%%6.13和5.22是这个代码
% clc
% clear all
% clf
% close all
%%%%%%%%%这个代码同时生成加噪声和不加噪声的多种形状的物体的数据%%%%%%%%%%%%%
function [elapsed_time] = dataset_gene(app_save_path, app_width, app_height, app_distance, app_length, app_pc_list, app_ps, app_freq, app_gn, app_gn_mean, app_gn_var, app_rn, app_rn_min, app_rn_max, app_mn, app_mn_order, app_mn_val, app_heitht_max, app_ball, app_ladder, app_ring, app_ring_tube, app_cone, app_cone_table, app_sc, app_gs, app_rm, app_prev, app_mask) 
tic;%开始计时
%%当列数为120时，T1=20，T2=24
scale=[app_height,app_width];
X = [1:1:1*scale(1)];Y = [1:1:1*scale(2)];
[x,y] = meshgrid(X,Y);
% 指定投影仪与相机的距离、相机距离基准面的高度
% 经高度调制后的球的相位
%%%%目前仿真的经验是被测物体的高度最好不超过200

D = app_distance;%相机和投影仪的距离
L = app_length;%相机到参考平面的距离mm
%%%%%%%四频六步相移法%%%%%%%%%
height=scale(1);%图片的高是矩阵的第一维
width=scale(2);%图片的宽是矩阵的第二维
T= app_pc_list; % 参数单位 mm
u=2*pi./T;
%%%%最终采用双频六步相移，频率比是4
%%%几频几步相移%%%%%%%%%
step=app_ps;
fre=app_freq;
k=0;%图片的序号

%%%%%%%%%%生成每个样品的次数%%%%%%
mm=1;%测试集20%合适300
%%%%%生成样本的总次数%%%%%%, app_ladder, app_ring, app_ring_tube, app_cone, app_cone_table, app_sc, app_gs
ball_k=app_ball;%生成半球的总次数
ladder_k=app_ladder;%生成阶梯的总次数
ring_k=app_ring;%生成ring的总次数
cone_k=app_cone_table;%生成锥台的总次数
ring_tube_k=app_ring_tube;%生成ring_tube的总次数
yuanzhui_k=app_cone;%生成圆锥的总次数
lk_yuanzhu_k=app_sc;%生成镂空圆柱的总次数
gaussian_surface_k=app_gs;%生成gaussian_surface的总次数
randon_matrix_k=app_rm;
app_total_count = app_ball+app_ladder+app_ring+app_cone_table+app_ring_tube+app_cone+app_sc+app_gs+app_rm;
%%%%%%%%%%%%%%%%%%%%被测物体最大高度%%%%%%%%%%%%%
height_all_max=app_heitht_max;

%%%%%%%%%%%%%%噪声值设置%%%%%%%%%%%
noisy=app_mn;%是否加噪，加噪为1，不加噪为0

k_h=app_mn_order;%高次谐波有几次
k_value=app_mn_val;

%%%%%%%%%%%%%%读写的路径%%%%%%%%%%%%%%
root=app_save_path;%根目录
%%%%%左投影目录
root1='\left\';%根目录
path_ori_geometory_l=[root root1 'geometry\'];
path_ori_in_l=[root root1 'grating\'];
path_ori_in_no_noisy_l=[root root1 'grating_no_noisy\'];
path_ori_hf_l=[root root1 'high_fre_grating\'];
path_ori_out_wrapped_low_l=[root root1 'phi_wrapped_low\'];
path_ori_out_wrapped_middle_l=[root root1 'phi_wrapped_middle\'];
path_ori_out_wrapped_high_l=[root root1 'phi_wrapped_high\'];
path_ori_out_unwrapped_l=[root root1 'phi_unwrapped\'];
path_ori_out_wrapped_low1_l=[root root1 'phi_wrapped_low_no_noisy\'];
path_ori_out_wrapped_middle1_l=[root root1 'phi_wrapped_middle_no_noisy\'];
path_ori_out_wrapped_high1_l=[root root1 'phi_wrapped_high_no_noisy\'];
path_ori_out_unwrapped1_l=[root root1 'phi_unwrapped_no_noisy\'];
path_ori_k_l=[root root1 'get_k_int\'];


% path generation生成投影文件夹
path_list = {path_ori_geometory_l, path_ori_hf_l, path_ori_in_l, ...
             path_ori_out_wrapped_low_l, path_ori_out_wrapped_high_l, path_ori_out_unwrapped_l, ...
             path_ori_in_no_noisy_l, path_ori_out_wrapped_low1_l, path_ori_out_wrapped_high1_l, ...
             path_ori_out_unwrapped1_l, path_ori_k_l,path_ori_out_wrapped_middle_l,path_ori_out_wrapped_middle1_l};

for i = 1:length(path_list)
    if ~exist(path_list{i},'dir')
        mkdir(path_list{i})
    end
end


%这一部分对球的相位进行编码，用四频六步相移法，加高斯噪声。
%但是作为标签的接包裹图，是由加噪声的相位光栅求解得到。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成球的几何体%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%这一部分对球的相位进行编码，用四频六步相移法，加高斯噪声。
%但是作为标签的接包裹图，是由加噪声的相位光栅求解得到。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成球的几何体%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%四个球的半径不同
%height_all_max=80;
ball_r=(10):(height_all_max);

ball_x1=(scale(1)/4-80):(scale(1)*3/4+80);%左上角，坐标原点在左上角，竖着的是x，横着的是y
ball_y1=(scale(2)/4-80):(scale(2)*3/4+80);



 for i=1:ball_k
      k=k+1;
      ballH = halfball(ball_r,ball_x1,ball_y1,scale); % 半径25mm 球心坐标（50,50）仿真假设一个像素为1mm
    %%把被测物体的形状图保存起来
    % 创建网格图
        if app_prev
            figure;
            subplot(3,2,1)
            mesh(ballH);
        title(strcat("geometry", num2str(k), "/", num2str(app_total_count)));
        colorbar
        view(-62,88);
        end
        
        file_name = [path_ori_geometory_l, num2str(k, '%06d'), '.png'];
        


      geometry=ballH/(height_all_max+10);%存成数据集的时候，数据进行归一化了，但是没有改变送进去包裹相位的高度
      path=[path_ori_geometory_l,num2str(k,'%06d'), '.mat'];
      save(path, 'geometry');


      geometry=ballH;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fringeModulation函数的参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        deltphi={};
        
        %%%%%左投影光栅保存噪声数据和非噪声数据
        path_in_left=[path_ori_in_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left,'dir')
                mkdir(path_in_left)
         end
         path_in_left_no_noisy=[path_ori_in_no_noisy_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left_no_noisy,'dir')
                mkdir(path_in_left_no_noisy)
         end
        if app_mask == 0
            mask = ones(height, width);
        else
            mask = generate_mask_with_black_circles(height, width);
        end
        if app_rn == 0
            noise = zeros(size(mask));
        else
            noise = -app_rn_min+ (app_rn_max + app_rn_max) * rand(size(mask));  % 范围在-0.05到0.05之间
        end
        
        for m=1:fre
            result= 2*pi*D*geometry./(L-geometry)/T(m); % T为正弦条纹的空间周期, 这里的./为点除
            deltphi{end+1} =result;
            u1=u(m);
            delta_phi=deltphi{m};
%             figure
%             imshow(delta_phi,[])
           for n=1:step
            phi=2*pi*(n-4)/step;
            [grating_left,img_left_nonoisy]=fringeModulation(u1,phi,delta_phi,height,width,noisy,k_h,k_value);%这一步已经对光栅进行归一化了

            img_left_nonoisy=img_left_nonoisy.*mask;         
            grating_gau_random=grating_left+ + noise;%加随机噪声
            if app_gn == 0
                grating_left=grating_gau_random;
            else 
                grating_left=imnoise(grating_gau_random, 'gaussian', app_gn_mean, app_gn_var); %加高斯噪声
            end
            grating_left=grating_left.*mask; 

            path=[path_in_left,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path,'grating_left');

            %%%%%%非噪声光栅保存
            path_nonoisy=[path_in_left_no_noisy,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path_nonoisy,'img_left_nonoisy');

           end

        end

        %%在最后一个频率的最后一步，保存加噪的高频光栅
            path=[path_ori_hf_l,num2str(k,'%06d'),'.mat'];%test的路径
            save(path,'grating_left');
            if app_prev
            subplot(3,2,2)
            imshow(grating_left,[])
            title("grating")
            end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解加噪过程中包裹相位和解包裹相位     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        [phi_wrapped_low,phi_wrapped_middle,phi_wrapped_high,phi_unwrapped,kh] = unwrap_mul_fre(fre,step,k,path_in_left,path_ori_out_wrapped_low_l, ...
            path_ori_out_wrapped_middle_l,path_ori_out_wrapped_high_l,path_ori_out_unwrapped_l) ;
        if app_prev
        subplot(3,1,3)
        imshow(phi_unwrapped,[])
        colorbar
        title('absolute phase')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解不加噪过程中的解包裹相位    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [phi_wrapped_low_no_noisy,phi_wrapped_middle_no_noisy,phi_wrapped_high_no_noisy,phi_unwrapped_no_noisy,kh] = unwrap_mul_fre(fre,step,k,path_in_left_no_noisy, ...
        path_ori_out_wrapped_low1_l,path_ori_out_wrapped_middle1_l,path_ori_out_wrapped_high1_l,path_ori_out_unwrapped1_l) ;
  if app_prev
        subplot(3,3,4)
        imshow(phi_wrapped_low_no_noisy,[])
        colorbar
         title('phi wrapped low')
        subplot(3,3,5)
        imshow(phi_wrapped_middle_no_noisy,[])
        colorbar
        title('phi wrapped middle')
        subplot(3,3,6)
        imshow(phi_wrapped_high_no_noisy,[])
        colorbar
        title('phi wrapped high')
        saveas(gcf, file_name);
  end
    end
  
disp('1 ball finished');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                   %%%%%   此部分是生成模拟阶梯   %%%%%%%%%%%

%球的部分生成之后图片的序数要k1=k+1；     
                                                                      

px_a=(20):(scale(1)/4);%竖着的是x轴，这代表x轴的起点。
py_a=(20):(scale(2)/4);%横着的事y轴，这是y轴的起点
step1_a=10:4:16;%步长
stepNum_a=1:6;%阶梯数，横着的是y轴
width_a=20:scale(2)-100;%x轴的范围

for i = 1:ladder_k

     k=k+1;

      ladderH = ladderFun(px_a,py_a,stepNum_a,step1_a,width_a,scale);
%%把被测物体的形状图保存起来
% 创建网格图

        if app_prev
            figure;
            subplot(3,2,1)
             mesh(ladderH);
    title(strcat("geometry", num2str(k), "/", num2str(app_total_count)));
    view(-62,88);
       colorbar;
        end
   
    % 保存图片
    file_name = [path_ori_geometory_l, num2str(k, '%06d'), '.png'];
    


      geometry=ladderH/(height_all_max+5);%存成数据集的时候，数据进行归一化了，但是没有改变送进去包裹相位的高度

      path=[path_ori_geometory_l,num2str(k,'%06d'), '.mat'];
      save(path, 'geometry');
      geometry=ladderH;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fringeModulation函数的参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%左投影光栅保存噪声数据和非噪声数据
        path_in_left=[path_ori_in_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left,'dir')
                mkdir(path_in_left)
         end
         path_in_left_no_noisy=[path_ori_in_no_noisy_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left_no_noisy,'dir')
                mkdir(path_in_left_no_noisy)
         end
        
        if app_mask == 0
            mask = ones(height, width);
        else
            mask = generate_mask_with_black_circles(height, width);
        end
        if app_rn == 0
            noise = zeros(size(mask));
        else
            noise = -app_rn_min+ (app_rn_max + app_rn_max) * rand(size(mask));  % 范围在-0.05到0.05之间
        end


   deltphi={};  
    for m=1:fre
        result= 2*pi*D*geometry./(L-geometry)/T(m); % T为正弦条纹的空间周期, 这里的./为点除
        deltphi{end+1} =result;
        u1=u(m);
        delta_phi=deltphi{m};
           for n=1:step
            phi=2*pi*(n-4)/step;
            [grating_left,img_left_nonoisy]=fringeModulation(u1,phi,delta_phi,height,width,noisy,k_h,k_value);%这一步已经对光栅进行归一化了

            img_left_nonoisy=img_left_nonoisy.*mask;         
            grating_gau_random=grating_left+ + noise;%加随机噪声
            if app_gn == 0
                grating_left=grating_gau_random;
            else 
                grating_left=imnoise(grating_gau_random, 'gaussian', app_gn_mean, app_gn_var); %加高斯噪声
            end
            grating_left=grating_left.*mask; 
            %%%%%噪声光栅保存
            path=[path_in_left,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path,'grating_left');

            %%%%%%非噪声光栅保存
            path_nonoisy=[path_in_left_no_noisy,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path_nonoisy,'img_left_nonoisy');

      end
    end
    %%在最后一个频率的最后一步，保存加噪的高频光栅
            path=[path_ori_hf_l,num2str(k,'%06d'),'.mat'];%test的路径
            save(path,'grating_left');
            if app_prev
            subplot(3,2,2)
            imshow(grating_left,[])
            title("grating")
            end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解加噪过程中包裹相位和解包裹相位     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        [phi_wrapped_low,phi_wrapped_middle,phi_wrapped_high,phi_unwrapped,kh] = unwrap_mul_fre(fre,step,k,path_in_left,path_ori_out_wrapped_low_l, ...
            path_ori_out_wrapped_middle_l,path_ori_out_wrapped_high_l,path_ori_out_unwrapped_l) ;
if app_prev
        subplot(3,1,3)
        imshow(phi_unwrapped,[])
        colorbar
        title('absolute phase')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解不加噪过程中的解包裹相位    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [phi_wrapped_low_no_noisy,phi_wrapped_middle_no_noisy,phi_wrapped_high_no_noisy,phi_unwrapped_no_noisy,kh] = unwrap_mul_fre(fre,step,k,path_in_left_no_noisy, ...
        path_ori_out_wrapped_low1_l,path_ori_out_wrapped_middle1_l,path_ori_out_wrapped_high1_l,path_ori_out_unwrapped1_l) ;
  if app_prev
        subplot(3,3,4)
        imshow(phi_wrapped_low_no_noisy,[])
        colorbar
         title('phi wrapped low')
        subplot(3,3,5)
        imshow(phi_wrapped_middle_no_noisy,[])
        colorbar
        title('phi wrapped middle')
        subplot(3,3,6)
        imshow(phi_wrapped_high_no_noisy,[])
        colorbar
        title('phi wrapped high')
        saveas(gcf, file_name);
  end
    end

  
disp('2 ladder finished');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成ring的几何体%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ring_rs=(20):(50);
ring_rl=(80):(150);
ring_x=(scale(1)/2-scale(1)*2/6):(scale(1)/2+scale(1)/4);
ring_y=(scale(2)/2-scale(2)*2/6):(scale(2)/2+scale(2)/4);
ring_h=10:height_all_max;

 for i=1:ring_k
      k=k+1;

      ringH=ring(ring_rl,ring_rs,ring_x,ring_y,ring_h,scale,height_all_max);

% 创建网格图
    
        if app_prev
            figure;
            subplot(3,2,1)
            mesh(ringH);
        title(strcat("geometry", num2str(k), "/", num2str(app_total_count)));
        view(-62,88);
         colorbar;
        end
        
        % 保存图片
        file_name = [path_ori_geometory_l, num2str(k, '%06d'), '.png'];
        
        % 关闭图形窗口
     
      geometry=ringH/(height_all_max+5);
      path=[path_ori_geometory_l,num2str(k,'%06d'), '.mat'];
      save(path, 'geometry');
      geometry=ringH;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fringeModulation函数的参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        deltphi={};

        %%%%%左投影光栅保存
        path_in_left=[path_ori_in_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left,'dir')
                mkdir(path_in_left)
         end
         path_in_left_no_noisy=[path_ori_in_no_noisy_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left_no_noisy,'dir')
                mkdir(path_in_left_no_noisy)
         end
        
        if app_mask == 0
            mask = ones(height, width);
        else
            mask = generate_mask_with_black_circles(height, width);
        end
        if app_rn == 0
            noise = zeros(size(mask));
        else
            noise = -app_rn_min+ (app_rn_max + app_rn_max) * rand(size(mask));  % 范围在-0.05到0.05之间
        end

        for m=1:fre
            result= 2*pi*D*geometry./(L-geometry)/T(m); % T为正弦条纹的空间周期, 这里的./为点除
            deltphi{end+1} =result;
            u1=u(m);
            delta_phi=deltphi{m};
%             figure
%             imshow(delta_phi,[])
           for n=1:step
            phi=2*pi*(n-4)/step;
            [grating_left,img_left_nonoisy]=fringeModulation(u1,phi,delta_phi,height,width,noisy,k_h,k_value);%这一步已经对光栅进行归一化了
            %非噪声数据乘以mask
            img_left_nonoisy=img_left_nonoisy.*mask; 
            %光栅加噪过程
            grating_gau_random=grating_left+ + noise;%加随机噪声
            if app_gn == 0
                grating_left=grating_gau_random;
            else 
                grating_left=imnoise(grating_gau_random, 'gaussian', app_gn_mean, app_gn_var); %加高斯噪声
            end
            grating_left=grating_left.*mask; 

            path=[path_in_left,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path,'grating_left');

            %%%%%%非噪声光栅保存
            path_nonoisy=[path_in_left_no_noisy,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path_nonoisy,'img_left_nonoisy');

           end

        end

        %%在最后一个频率的最后一步，保存加噪的高频光栅
            path=[path_ori_hf_l,num2str(k,'%06d'),'.mat'];%test的路径
            save(path,'grating_left');
            if app_prev
            subplot(3,2,2)
            imshow(grating_left,[])
            title("grating")
            end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解加噪过程中包裹相位和解包裹相位     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        [phi_wrapped_low,phi_wrapped_middle,phi_wrapped_high,phi_unwrapped,kh] = unwrap_mul_fre(fre,step,k,path_in_left,path_ori_out_wrapped_low_l, ...
            path_ori_out_wrapped_middle_l,path_ori_out_wrapped_high_l,path_ori_out_unwrapped_l) ;
if app_prev
        subplot(3,1,3)
        imshow(phi_unwrapped,[])
        colorbar
        title('absolute phase')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解不加噪过程中的解包裹相位    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [phi_wrapped_low_no_noisy,phi_wrapped_middle_no_noisy,phi_wrapped_high_no_noisy,phi_unwrapped_no_noisy,kh] = unwrap_mul_fre(fre,step,k,path_in_left_no_noisy, ...
        path_ori_out_wrapped_low1_l,path_ori_out_wrapped_middle1_l,path_ori_out_wrapped_high1_l,path_ori_out_unwrapped1_l) ;
  
  if app_prev
        subplot(3,3,4)
        imshow(phi_wrapped_low_no_noisy,[])
        colorbar
         title('phi wrapped low')
        subplot(3,3,5)
        imshow(phi_wrapped_middle_no_noisy,[])
        colorbar
        title('phi wrapped middle')
        subplot(3,3,6)
        imshow(phi_wrapped_high_no_noisy,[])
        colorbar
        title('phi wrapped high')
        saveas(gcf, file_name);
  end
    end
  
disp('3 ring finished');


%%%%%生成锥台%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成锥台的几何体%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        cone_rs=(10):(50);
        cone_rl=(60):(160);
        cone_x=(scale(1)/2-scale(1)*2/6):(scale(1)/2+scale(1)/4);
        cone_y=(scale(2)/2-scale(2)*2/6):(scale(2)/2+scale(2)/4);
        cone_h=10:height_all_max;

 for i=1:cone_k
      k=k+1;

      coneH=cone(cone_rl,cone_rs,cone_x,cone_y,cone_h,scale(1),scale(2));
         if app_prev
            figure;
            subplot(3,2,1)
            mesh(coneH);
    title(strcat("geometry", num2str(k), "/", num2str(app_total_count)));
    view(-62,88);
       colorbar;
        end
    
    % 保存图片
    file_name = [path_ori_geometory_l, num2str(k, '%06d'), '.png'];
    

      geometry=coneH/(height_all_max+5);
      path=[path_ori_geometory_l,num2str(k,'%06d'), '.mat'];
      save(path, 'geometry');
      geometry=coneH;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fringeModulation函数的参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        deltphi={};

        path_in_left=[path_ori_in_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left,'dir')
                mkdir(path_in_left)
         end
         path_in_left_no_noisy=[path_ori_in_no_noisy_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left_no_noisy,'dir')
                mkdir(path_in_left_no_noisy)
         end
        
        if app_mask == 0
            mask = ones(height, width);
        else
            mask = generate_mask_with_black_circles(height, width);
        end
        if app_rn == 0
            noise = zeros(size(mask));
        else
            noise = -app_rn_min+ (app_rn_max + app_rn_max) * rand(size(mask));  % 范围在-0.05到0.05之间
        end
        
        for m=1:fre
            result= 2*pi*D*geometry./(L-geometry)/T(m); % T为正弦条纹的空间周期, 这里的./为点除
            deltphi{end+1} =result;
            u1=u(m);
            delta_phi=deltphi{m};
%             figure
%             imshow(delta_phi,[])
            for n=1:step
            phi=2*pi*(n-4)/step;
            [grating_left,img_left_nonoisy]=fringeModulation(u1,phi,delta_phi,height,width,noisy,k_h,k_value);%这一步已经对光栅进行归一化了

            img_left_nonoisy=img_left_nonoisy.*mask;         
            grating_gau_random=grating_left+ + noise;%加随机噪声
            if app_gn == 0
                grating_left=grating_gau_random;
            else 
                grating_left=imnoise(grating_gau_random, 'gaussian', app_gn_mean, app_gn_var); %加高斯噪声
            end
            grating_left=grating_left.*mask; 
            %%%%%噪声光栅保存
            path=[path_in_left,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path,'grating_left');

            %%%%%%非噪声光栅保存
            path_nonoisy=[path_in_left_no_noisy,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path_nonoisy,'img_left_nonoisy');

      end
    end
    %%在最后一个频率的最后一步，保存加噪的高频光栅
            path=[path_ori_hf_l,num2str(k,'%06d'),'.mat'];%test的路径
            save(path,'grating_left');
if app_prev
            subplot(3,2,2)
            imshow(grating_left,[])
            title("grating")
            end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解加噪过程中包裹相位和解包裹相位     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        [phi_wrapped_low,phi_wrapped_middle,phi_wrapped_high,phi_unwrapped,kh] = unwrap_mul_fre(fre,step,k,path_in_left,path_ori_out_wrapped_low_l, ...
            path_ori_out_wrapped_middle_l,path_ori_out_wrapped_high_l,path_ori_out_unwrapped_l) ;
if app_prev
        subplot(3,1,3)
        imshow(phi_unwrapped,[])
        colorbar
        title('absolute phase')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解不加噪过程中的解包裹相位    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [phi_wrapped_low_no_noisy,phi_wrapped_middle_no_noisy,phi_wrapped_high_no_noisy,phi_unwrapped_no_noisy,kh] = unwrap_mul_fre(fre,step,k,path_in_left_no_noisy, ...
        path_ori_out_wrapped_low1_l,path_ori_out_wrapped_middle1_l,path_ori_out_wrapped_high1_l,path_ori_out_unwrapped1_l) ;
  
  if app_prev
        subplot(3,3,4)
        imshow(phi_wrapped_low_no_noisy,[])
        colorbar
         title('phi wrapped low')
        subplot(3,3,5)
        imshow(phi_wrapped_middle_no_noisy,[])
        colorbar
        title('phi wrapped middle')
        subplot(3,3,6)
        imshow(phi_wrapped_high_no_noisy,[])
        colorbar
        title('phi wrapped high')
        saveas(gcf, file_name);
  end
    end

  
disp('4 cone finished');


%%%%%%镂空半球%%%%%%%%%%%%
%这一部分对球的相位进行编码，用四频六步相移法，加高斯噪声。
%但是作为标签的接包裹图，是由加噪声的相位光栅求解得到。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5生成镂空半球的几何体%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ring_tube_rs=(10):(30);
ring_tube_rl=(40):(height_all_max);
ring_tube_x=(scale(1)/2-scale(1)*2/6):(scale(1)/2+scale(1)/4);
ring_tube_y=(scale(2)/2-scale(2)*2/6):(scale(2)/2+scale(2)/4);

   for i=1:ring_tube_k
      k=k+1;
      ring_tubeH=ring_tube(ring_tube_rl,ring_tube_rs,ring_tube_x,ring_tube_y,scale(1),scale(2));
            if app_prev
            figure;
            subplot(3,2,1)
            mesh(ring_tubeH);
    title(strcat("geometry", num2str(k), "/", num2str(app_total_count)));
    view(-62,88);
       colorbar;
        end
    
    % 保存图片
    file_name = [path_ori_geometory_l, num2str(k, '%06d'), '.png'];

      geometry=ring_tubeH/(height_all_max+5);
      path=[path_ori_geometory_l,num2str(k,'%06d'), '.mat'];
      save(path, 'geometry');
      geometry=ring_tubeH;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fringeModulation函数的参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        deltphi={};

        path_in_left=[path_ori_in_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left,'dir')
                mkdir(path_in_left)
         end
         path_in_left_no_noisy=[path_ori_in_no_noisy_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left_no_noisy,'dir')
                mkdir(path_in_left_no_noisy)
         end
        
        if app_mask == 0
            mask = ones(height, width);
        else
            mask = generate_mask_with_black_circles(height, width);
        end
        if app_rn == 0
            noise = zeros(size(mask));
        else
            noise = -app_rn_min+ (app_rn_max + app_rn_max) * rand(size(mask));  % 范围在-0.05到0.05之间
        end
     
        for m=1:fre
            result= 2*pi*D*geometry./(L-geometry)/T(m); % T为正弦条纹的空间周期, 这里的./为点除
            deltphi{end+1} =result;
            u1=u(m);
            delta_phi=deltphi{m};
%             figure
        for n=1:step
            phi=2*pi*(n-4)/step;
            [grating_left,img_left_nonoisy]=fringeModulation(u1,phi,delta_phi,height,width,noisy,k_h,k_value);%这一步已经对光栅进行归一化了

            img_left_nonoisy=img_left_nonoisy.*mask;         
            grating_gau_random=grating_left+ + noise;%加随机噪声
            if app_gn == 0
                grating_left=grating_gau_random;
            else 
                grating_left=imnoise(grating_gau_random, 'gaussian', app_gn_mean, app_gn_var); %加高斯噪声
            end
            grating_left=grating_left.*mask; 
            %%%%%噪声光栅保存
            path=[path_in_left,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path,'grating_left');

            %%%%%%非噪声光栅保存
            path_nonoisy=[path_in_left_no_noisy,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path_nonoisy,'img_left_nonoisy');

      end
    end
    %%在最后一个频率的最后一步，保存加噪的高频光栅
            path=[path_ori_hf_l,num2str(k,'%06d'),'.mat'];%test的路径
            save(path,'grating_left');
if app_prev
            subplot(3,2,2)
            imshow(grating_left,[])
            title("grating")
            end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解加噪过程中包裹相位和解包裹相位     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        [phi_wrapped_low,phi_wrapped_middle,phi_wrapped_high,phi_unwrapped,kh] = unwrap_mul_fre(fre,step,k,path_in_left,path_ori_out_wrapped_low_l, ...
            path_ori_out_wrapped_middle_l,path_ori_out_wrapped_high_l,path_ori_out_unwrapped_l) ;
if app_prev
        subplot(3,1,3)
        imshow(phi_unwrapped,[])
        colorbar
        title('absolute phase')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解不加噪过程中的解包裹相位    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [phi_wrapped_low_no_noisy,phi_wrapped_middle_no_noisy,phi_wrapped_high_no_noisy,phi_unwrapped_no_noisy,kh] = unwrap_mul_fre(fre,step,k,path_in_left_no_noisy, ...
        path_ori_out_wrapped_low1_l,path_ori_out_wrapped_middle1_l,path_ori_out_wrapped_high1_l,path_ori_out_unwrapped1_l) ;
  
  if app_prev
        subplot(3,3,4)
        imshow(phi_wrapped_low_no_noisy,[])
        colorbar
         title('phi wrapped low')
        subplot(3,3,5)
        imshow(phi_wrapped_middle_no_noisy,[])
        colorbar
        title('phi wrapped middle')
        subplot(3,3,6)
        imshow(phi_wrapped_high_no_noisy,[])
        colorbar
        title('phi wrapped high')
        saveas(gcf, file_name);
  end
    end

  
disp('5 hollow_sphere finished');

%%%%%%生成圆锥%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成圆锥的几何体%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%height_all_max
yuanzhui_h=10:height_all_max;
yuanzhui_rl=(scale/10):(scale/2-40);
yuanzhui_x=(scale(1)/2-scale(1)*2/6):(scale(1)/2+scale(1)/4);
yuanzhui_y=(scale(2)/2-scale(2)*2/6):(scale(2)/2+scale(2)/4);



 for i=1:yuanzhui_k
      k=k+1;
      yuanzhuiH=yuanzhui(yuanzhui_rl,yuanzhui_x,yuanzhui_y,yuanzhui_h,scale(1),scale(2));
    if app_prev
            figure;
            subplot(3,2,1)
            mesh(yuanzhuiH);
    title(strcat("geometry", num2str(k), "/", num2str(app_total_count)));
    view(-62,88);
    colorbar;
    end
    
    % 保存图片
    file_name = [path_ori_geometory_l, num2str(k, '%06d'), '.png'];



      geometry=yuanzhuiH/(height_all_max+5);
      path=[path_ori_geometory_l,num2str(k,'%06d'), '.mat'];
      save(path, 'geometry');
      geometry=yuanzhuiH;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fringeModulation函数的参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        deltphi={};

        path_in_left=[path_ori_in_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left,'dir')
                mkdir(path_in_left)
         end
         path_in_left_no_noisy=[path_ori_in_no_noisy_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left_no_noisy,'dir')
                mkdir(path_in_left_no_noisy)
         end
        
        if app_mask == 0
            mask = ones(height, width);
        else
            mask = generate_mask_with_black_circles(height, width);
        end
        if app_rn == 0
            noise = zeros(size(mask));
        else
            noise = -app_rn_min+ (app_rn_max + app_rn_max) * rand(size(mask));  % 范围在-0.05到0.05之间
        end
        
        for m=1:fre
            result= 2*pi*D*geometry./(L-geometry)/T(m); % T为正弦条纹的空间周期, 这里的./为点除
            deltphi{end+1} =result;
            u1=u(m);
            delta_phi=deltphi{m};
%             figure
%             imshow(delta_phi,[])
            for n=1:step
            phi=2*pi*(n-4)/step;
            [grating_left,img_left_nonoisy]=fringeModulation(u1,phi,delta_phi,height,width,noisy,k_h,k_value);%这一步已经对光栅进行归一化了

            img_left_nonoisy=img_left_nonoisy.*mask;         
            grating_gau_random=grating_left+ + noise;%加随机噪声
            if app_gn == 0
                grating_left=grating_gau_random;
            else 
                grating_left=imnoise(grating_gau_random, 'gaussian', app_gn_mean, app_gn_var); %加高斯噪声
            end
            grating_left=grating_left.*mask; 
            %%%%%噪声光栅保存
            path=[path_in_left,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path,'grating_left');

            %%%%%%非噪声光栅保存
            path_nonoisy=[path_in_left_no_noisy,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path_nonoisy,'img_left_nonoisy');

      end
    end
    %%在最后一个频率的最后一步，保存加噪的高频光栅
            path=[path_ori_hf_l,num2str(k,'%06d'),'.mat'];%test的路径
            save(path,'grating_left');
if app_prev
            subplot(3,2,2)
            imshow(grating_left,[])
            title("grating")
            end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解加噪过程中包裹相位和解包裹相位     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        [phi_wrapped_low,phi_wrapped_middle,phi_wrapped_high,phi_unwrapped,kh] = unwrap_mul_fre(fre,step,k,path_in_left,path_ori_out_wrapped_low_l, ...
            path_ori_out_wrapped_middle_l,path_ori_out_wrapped_high_l,path_ori_out_unwrapped_l) ;
if app_prev
        subplot(3,1,3)
        imshow(phi_unwrapped,[])
        colorbar
        title('absolute phase')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解不加噪过程中的解包裹相位    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [phi_wrapped_low_no_noisy,phi_wrapped_middle_no_noisy,phi_wrapped_high_no_noisy,phi_unwrapped_no_noisy,kh] = unwrap_mul_fre(fre,step,k,path_in_left_no_noisy, ...
        path_ori_out_wrapped_low1_l,path_ori_out_wrapped_middle1_l,path_ori_out_wrapped_high1_l,path_ori_out_unwrapped1_l) ;
  
  if app_prev
        subplot(3,3,4)
        imshow(phi_wrapped_low_no_noisy,[])
        colorbar
         title('phi wrapped low')
        subplot(3,3,5)
        imshow(phi_wrapped_middle_no_noisy,[])
        colorbar
        title('phi wrapped middle')
        subplot(3,3,6)
        imshow(phi_wrapped_high_no_noisy,[])
        colorbar
        title('phi wrapped high')
        saveas(gcf, file_name);
  end
    end

  
disp('6 cone finished');

%%%%镂空圆柱
                                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成镂空半圆柱的几何体%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lk_yuanzhu_l=(20):2:(scale(2)/2);%半圆柱的长
lk_yuanzhu_r=(30):(height_all_max);
lk_yuanzhu_rs=(20):(50);


for i=1:lk_yuanzhu_k
      k=k+1;
      lk_yuanzhuH=lk_yuanzhu(lk_yuanzhu_l,lk_yuanzhu_r,lk_yuanzhu_rs,scale(1),scale(2));
            if app_prev
            figure;
            subplot(3,2,1)
             mesh(lk_yuanzhuH);
    title(strcat("geometry", num2str(k), "/", num2str(app_total_count)));
    view(-62,88);
       colorbar;
        end
   
    % 保存图片
    file_name = [path_ori_geometory_l, num2str(k, '%06d'), '.png'];
    


      geometry=lk_yuanzhuH/(height_all_max+5);
      path=[path_ori_geometory_l,num2str(k,'%06d'), '.mat'];
      save(path, 'geometry');
      geometry=lk_yuanzhuH;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fringeModulation函数的参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         deltphi={};

        path_in_left=[path_ori_in_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left,'dir')
                mkdir(path_in_left)
         end
         path_in_left_no_noisy=[path_ori_in_no_noisy_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left_no_noisy,'dir')
                mkdir(path_in_left_no_noisy)
         end
        
        if app_mask == 0
            mask = ones(height, width);
        else
            mask = generate_mask_with_black_circles(height, width);
        end
        if app_rn == 0
            noise = zeros(size(mask));
        else
            noise = -app_rn_min+ (app_rn_max + app_rn_max) * rand(size(mask));  % 范围在-0.05到0.05之间
        end
        
        for m=1:fre
            result= 2*pi*D*geometry./(L-geometry)/T(m); % T为正弦条纹的空间周期, 这里的./为点除
            deltphi{end+1} =result;
            u1=u(m);
            delta_phi=deltphi{m};
%             figure
%             imshow(delta_phi,[])
            for n=1:step
            phi=2*pi*(n-4)/step;
            [grating_left,img_left_nonoisy]=fringeModulation(u1,phi,delta_phi,height,width,noisy,k_h,k_value);%这一步已经对光栅进行归一化了

            img_left_nonoisy=img_left_nonoisy.*mask;         
            grating_gau_random=grating_left+ + noise;%加随机噪声
            if app_gn == 0
                grating_left=grating_gau_random;
            else 
                grating_left=imnoise(grating_gau_random, 'gaussian', app_gn_mean, app_gn_var); %加高斯噪声
            end
            grating_left=grating_left.*mask; 
            %%%%%噪声光栅保存
            path=[path_in_left,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path,'grating_left');

            %%%%%%非噪声光栅保存
            path_nonoisy=[path_in_left_no_noisy,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path_nonoisy,'img_left_nonoisy');

      end
    end
    %%在最后一个频率的最后一步，保存加噪的高频光栅
            path=[path_ori_hf_l,num2str(k,'%06d'),'.mat'];%test的路径
            save(path,'grating_left');
if app_prev
            subplot(3,2,2)
            imshow(grating_left,[])
            title("grating")
            end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解加噪过程中包裹相位和解包裹相位     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        [phi_wrapped_low,phi_wrapped_middle,phi_wrapped_high,phi_unwrapped,kh] = unwrap_mul_fre(fre,step,k,path_in_left,path_ori_out_wrapped_low_l, ...
            path_ori_out_wrapped_middle_l,path_ori_out_wrapped_high_l,path_ori_out_unwrapped_l) ;
if app_prev
        subplot(3,1,3)
        imshow(phi_unwrapped,[])
        colorbar
        title('absolute phase')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解不加噪过程中的解包裹相位    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [phi_wrapped_low_no_noisy,phi_wrapped_middle_no_noisy,phi_wrapped_high_no_noisy,phi_unwrapped_no_noisy,kh] = unwrap_mul_fre(fre,step,k,path_in_left_no_noisy, ...
        path_ori_out_wrapped_low1_l,path_ori_out_wrapped_middle1_l,path_ori_out_wrapped_high1_l,path_ori_out_unwrapped1_l) ;
  
  if app_prev
        subplot(3,3,4)
        imshow(phi_wrapped_low_no_noisy,[])
        colorbar
         title('phi wrapped low')
        subplot(3,3,5)
        imshow(phi_wrapped_middle_no_noisy,[])
        colorbar
        title('phi wrapped middle')
        subplot(3,3,6)
        imshow(phi_wrapped_high_no_noisy,[])
        colorbar
        title('phi wrapped high')
        saveas(gcf, file_name);
  end
    end

  
disp('7 hollow_cylinder finished');

%这个部分生成高斯曲面，曲面的峰值数是随机的，g_num_k_a=2:6;%%高斯函数的峰值个数。
                                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成高斯曲面%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g_num_k = 10;%%高斯函数的峰值个数。
%gaussian_surfaceH=0;
mux_a=1:(scale(1));%高斯曲面的中心位置
muy_a=(1):(scale(2));
sigmax_a=scale(1)/8:scale(1)*2/8;%高斯曲面的范围
sigmay_a=scale(2)/8:scale(2)*2/8;
amplitude_a=10:20;

 for i=1:gaussian_surface_k
      k=k+1;
      g_sH=0;
  for i=1:g_num_k%十五个高斯函数累加
         g_s=gaussian_surface(mux_a,muy_a,sigmax_a,sigmay_a,amplitude_a,scale(1),scale(2) );
         
         g_sH=g_sH+g_s;
   end
     for i = 1:size(g_sH,1)
       for j = 1:size(g_sH,2) % 如果元素的值大于60，用50-60之间的随机数代替
        if g_sH(i,j) > height_all_max % 生成50-60之间的随机数并替换该元素
            g_sH(i,j) = height_all_max;
        end
       end
     end

           if app_prev
            figure;
            subplot(3,2,1)
            mesh(g_sH);
    title(strcat("geometry", num2str(k), "/", num2str(app_total_count)));
    view(-23,71);
       colorbar;
        end
    
    % 保存图片
    file_name = [path_ori_geometory_l, num2str(k, '%06d'), '.png'];
  gaussian_surfaceH=g_sH;
  geometry=gaussian_surfaceH/(height_all_max+5);
  path=[path_ori_geometory_l,num2str(k,'%06d'), '.mat'];
  save(path, 'geometry');
   geometry=gaussian_surfaceH;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fringeModulation函数的参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  deltphi={};

        path_in_left=[path_ori_in_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left,'dir')
                mkdir(path_in_left)
         end
         path_in_left_no_noisy=[path_ori_in_no_noisy_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left_no_noisy,'dir')
                mkdir(path_in_left_no_noisy)
         end
        
        if app_mask == 0
            mask = ones(height, width);
        else
            mask = generate_mask_with_black_circles(height, width);
        end
        if app_rn == 0
            noise = zeros(size(mask));
        else
            noise = -app_rn_min+ (app_rn_max + app_rn_max) * rand(size(mask));  % 范围在-0.05到0.05之间
        end
        
        for m=1:fre
            result= 2*pi*D*geometry./(L-geometry)/T(m); % T为正弦条纹的空间周期, 这里的./为点除
            deltphi{end+1} =result;
            u1=u(m);
            delta_phi=deltphi{m};
%             figure
%             imshow(delta_phi,[])
            for n=1:step
            phi=2*pi*(n-4)/step;
            [grating_left,img_left_nonoisy]=fringeModulation(u1,phi,delta_phi,height,width,noisy,k_h,k_value);%这一步已经对光栅进行归一化了

            img_left_nonoisy=img_left_nonoisy.*mask;         
            grating_gau_random=grating_left+ + noise;%加随机噪声
            if app_gn == 0
                grating_left=grating_gau_random;
            else 
                grating_left=imnoise(grating_gau_random, 'gaussian', app_gn_mean, app_gn_var); %加高斯噪声
            end
            grating_left=grating_left.*mask; 
            %%%%%噪声光栅保存
            path=[path_in_left,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path,'grating_left');

            %%%%%%非噪声光栅保存
            path_nonoisy=[path_in_left_no_noisy,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path_nonoisy,'img_left_nonoisy');

      end
    end
    %%在最后一个频率的最后一步，保存加噪的高频光栅
            path=[path_ori_hf_l,num2str(k,'%06d'),'.mat'];%test的路径
            save(path,'grating_left');
if app_prev
            subplot(3,2,2)
            imshow(grating_left,[])
            title("grating")
            end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解加噪过程中包裹相位和解包裹相位     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        [phi_wrapped_low,phi_wrapped_middle,phi_wrapped_high,phi_unwrapped,kh] = unwrap_mul_fre(fre,step,k,path_in_left,path_ori_out_wrapped_low_l, ...
            path_ori_out_wrapped_middle_l,path_ori_out_wrapped_high_l,path_ori_out_unwrapped_l) ;
if app_prev
        subplot(3,1,3)
        imshow(phi_unwrapped,[])
        colorbar
        title('absolute phase')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解不加噪过程中的解包裹相位    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [phi_wrapped_low_no_noisy,phi_wrapped_middle_no_noisy,phi_wrapped_high_no_noisy,phi_unwrapped_no_noisy,kh] = unwrap_mul_fre(fre,step,k,path_in_left_no_noisy, ...
        path_ori_out_wrapped_low1_l,path_ori_out_wrapped_middle1_l,path_ori_out_wrapped_high1_l,path_ori_out_unwrapped1_l) ;
  if app_prev
        subplot(3,3,4)
        imshow(phi_wrapped_low_no_noisy,[])
        colorbar
         title('phi wrapped low')
        subplot(3,3,5)
        imshow(phi_wrapped_middle_no_noisy,[])
        colorbar
        title('phi wrapped middle')
        subplot(3,3,6)
        imshow(phi_wrapped_high_no_noisy,[])
        colorbar
        title('phi wrapped high')
        saveas(gcf, file_name);
  end
       
    end

  
disp('8 gaussian_surface finised');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成随机矩阵%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


     
 for i=1:randon_matrix_k
        k=k+1;
%         disp(k);
        randon_matrixH=randon_matrix(height_all_max,scale(1),scale(2)) ;
                if app_prev
            figure;
            subplot(3,2,1)
             mesh(randon_matrixH);
        title(strcat("geometry", num2str(k), "/", num2str(app_total_count)));
        view(-23,71);
        colorbar;
        end
       
        file_name = [path_ori_geometory_l, num2str(k, '%06d'), '.png'];

        geometry=randon_matrixH/(height_all_max+5);
        path=[path_ori_geometory_l,num2str(k,'%06d'), '.mat'];
        save(path, 'geometry');
        geometry=randon_matrixH;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fringeModulation函数的参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        
       deltphi={};

        path_in_left=[path_ori_in_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left,'dir')
                mkdir(path_in_left)
         end
         path_in_left_no_noisy=[path_ori_in_no_noisy_l,sprintf('%06d', k),filesep];
         if ~exist(path_in_left_no_noisy,'dir')
                mkdir(path_in_left_no_noisy)
         end
        
        if app_mask == 0
            mask = ones(height, width);
        else
            mask = generate_mask_with_black_circles(height, width);
        end
        if app_rn == 0
            noise = zeros(size(mask));
        else
            noise = -app_rn_min+ (app_rn_max + app_rn_max) * rand(size(mask));  % 范围在-0.05到0.05之间
        end
        
        for m=1:fre
            result= 2*pi*D*geometry./(L-geometry)/T(m); % T为正弦条纹的空间周期, 这里的./为点除
            deltphi{end+1} =result;
            u1=u(m);
            delta_phi=deltphi{m};
%             figure
%             imshow(delta_phi,[])
            for n=1:step
            phi=2*pi*(n-4)/step;
            [grating_left,img_left_nonoisy]=fringeModulation(u1,phi,delta_phi,height,width,noisy,k_h,k_value);%这一步已经对光栅进行归一化了

            img_left_nonoisy=img_left_nonoisy.*mask;         
            grating_gau_random=grating_left+ + noise;%加随机噪声
            if app_gn == 0
                grating_left=grating_gau_random;
            else 
                grating_left=imnoise(grating_gau_random, 'gaussian', app_gn_mean, app_gn_var); %加高斯噪声
            end
            grating_left=grating_left.*mask; 
            %%%%%噪声光栅保存
            path=[path_in_left,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path,'grating_left');

            %%%%%%非噪声光栅保存
            path_nonoisy=[path_in_left_no_noisy,num2str(m), '_',num2str(n),'.mat'];%test的路径
            save(path_nonoisy,'img_left_nonoisy');

      end
    end
    %%在最后一个频率的最后一步，保存加噪的高频光栅
            path=[path_ori_hf_l,num2str(k,'%06d'),'.mat'];%test的路径
            save(path,'grating_left');
            if app_prev
            subplot(3,2,2)
            imshow(grating_left,[])
            title("grating")
            end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解加噪过程中包裹相位和解包裹相位     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        [phi_wrapped_low,phi_wrapped_middle,phi_wrapped_high,phi_unwrapped,kh] = unwrap_mul_fre(fre,step,k,path_in_left,path_ori_out_wrapped_low_l, ...
            path_ori_out_wrapped_middle_l,path_ori_out_wrapped_high_l,path_ori_out_unwrapped_l) ;
if app_prev
        subplot(3,1,3)
        imshow(phi_unwrapped,[])
        colorbar
        title('absolute phase')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解不加噪过程中的解包裹相位    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [phi_wrapped_low_no_noisy,phi_wrapped_middle_no_noisy,phi_wrapped_high_no_noisy,phi_unwrapped_no_noisy,kh] = unwrap_mul_fre(fre,step,k,path_in_left_no_noisy, ...
        path_ori_out_wrapped_low1_l,path_ori_out_wrapped_middle1_l,path_ori_out_wrapped_high1_l,path_ori_out_unwrapped1_l) ;
  
  if app_prev
        subplot(3,3,4)
        imshow(phi_wrapped_low_no_noisy,[])
        colorbar
         title('phi wrapped low')
        subplot(3,3,5)
        imshow(phi_wrapped_middle_no_noisy,[])
        colorbar
        title('phi wrapped middle')
        subplot(3,3,6)
        imshow(phi_wrapped_high_no_noisy,[])
        colorbar
        title('phi wrapped high')
        saveas(gcf, file_name);
  end
    end

elapsed_time = toc;
fprintf('finished in %.4f minutes\n', elapsed_time/60);
  

