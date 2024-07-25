
clc
clear all
clf
%%%%%%%%%有噪声的数据，此代码和无噪声一起使用此代码产生阶梯的仿真

%%当列数为120时，T1=20，T2=24
scale = 512;
X = [1:1:1*scale];Y = [1:1:1*scale];

[x,y] = meshgrid(X,Y);
% 指定投影仪与相机的距离、相机距离基准面的高度
% 经高度调制后的球的相位
%%%%目前仿真的经验是被测物体的高度最好不超过200

D = 200;
L = 300;
%%%%%%%四频六步相移法%%%%%%%%%
height=scale;%图片的高是矩阵的第一维
width=scale;%图片的宽是矩阵的第二维
T= [width,width/2,width/8,width/32]; % 参数单位 mm
u=[2*pi/T(1),2*pi/T(2),2*pi/T(3),2*pi/T(4)];
%%%几频几步相移%%%%%%%%%
step=6;
fre=4;
k=0;%图片的序号

%这一部分对球的相位进行编码，用四频六步相移法，加高斯噪声。
%但是作为标签的接包裹图，是由加噪声的相位光栅求解得到。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成球的几何体%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%这一部分对球的相位进行编码，用四频六步相移法，加高斯噪声。
%但是作为标签的接包裹图，是由加噪声的相位光栅求解得到。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成高斯曲面%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gaussian_surface_k=2;%生成gaussian_surface的总次数
g_num_k_a=2:6;%%高斯函数的峰值个数。
g_num_k = g_num_k_a(randi(length(g_num_k_a)));
% g_num_k =1;
gaussian_surfaceH=0;

gaussian_mux_a=(100):(scale-100);%高斯曲面的中心位置
gaussian_muy_a=(100):(scale-100);
gaussian_sigmax_a=10:80;%高斯曲面的范围
gaussian_sigmay_a=10:80;
gaussian_amplitude_a=20:60;


 for i=1:gaussian_surface_k
      
       for i=1:g_num_k
      k=k+1;
      gaussian_mux = gaussian_mux_a(randi(length(gaussian_mux_a)));
      gaussian_muy = gaussian_muy_a(randi(length(gaussian_muy_a)));
      gaussian_sigmax = gaussian_sigmax_a(randi(length(gaussian_sigmax_a)));
      gaussian_sigmay = gaussian_sigmay_a(randi(length(gaussian_sigmay_a)));
      gaussian_amplitude = gaussian_amplitude_a(randi(length(gaussian_amplitude_a)));
  

      g_s=gaussian_surface(gaussian_mux,gaussian_muy,gaussian_sigmax,gaussian_sigmay,gaussian_amplitude,scale );
      gaussian_surfaceH=gaussian_surfaceH+g_s;
 

      figure;
      mesh(x,y,gaussian_surfaceH);
%  path=['E:\liyimingPCL\博士课题\实验记录-源代码+过程\12.08UNET相位预测仿真\' ...
%      'Simulation of Reconstruction\simu_result\half_ball_test_1\depth_noisy\',num2str(k),'.mat']; 
%  save(path, 'gaussian_surfaceH');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fgaussian_surfaceeModulation函数的参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        noisy=1;%是否加噪，加噪为1，不加噪为0
        noisy_value_guassian=0.001;%噪声值
        k_h=6;%高次谐波有几次
        k_value=0.1;

        deltphi={};
        
        
      for m=1:fre
            result= 2*pi*D*gaussian_surfaceH./(L-gaussian_surfaceH)/T(m); % T为正弦条纹的空间周期, 这里的./为点除
            deltphi{end+1} =result;
            u1=u(m);
            delta_phi=deltphi{m};
         for n=1:step
            phi=2*pi*(n-4)/step;
          
            high_fre_grating=fringeModulation(u1,phi,delta_phi,height,width,noisy,noisy_value_guassian,k_h,k_value);
%           figure;
%             mesh(x,y,high_fre_grating);
            path=['E:\liyimingPCL\博士课题\实验记录-源代码+过程\3.14-基本工件几何体仿真\code\new\test\',num2str(m), '_',num2str(n),'.mat'];
            save(path, 'high_fre_grating');
        %  save(path, 'First');
          end
       end
          path=['E:\liyimingPCL\博士课题\实验记录-源代码+过程\3.14-基本工件几何体仿真\code\new\high_fre_grating\',num2str(k), '.mat'];
         save(path, 'high_fre_grating');
        
        % figure(1);imshow(First);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解加噪过程中包裹相位和解包裹相位     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        path_in='E:\liyimingPCL\博士课题\实验记录-源代码+过程\3.14-基本工件几何体仿真\code\new\test\';
        path_out_wrapped_low='E:\liyimingPCL\博士课题\实验记录-源代码+过程\3.14-基本工件几何体仿真\code\new\phi_wrapped_low\';
        path_out_wrapped_high='E:\liyimingPCL\博士课题\实验记录-源代码+过程\3.14-基本工件几何体仿真\code\new\phi_wrapped_high\';
        path_out_unwrapped='E:\liyimingPCL\博士课题\实验记录-源代码+过程\3.14-基本工件几何体仿真\code\new\phi_unwrapped\';
        
        [phi_wrapped_low,phi_wrapped_high,phi_unwrapped] = unwrap_mul_fre(4,6,k,path_in,path_out_wrapped_low,path_out_wrapped_high,path_out_unwrapped) ;
        %返回低频/高频包裹条纹和解包裹条纹
        % m为频率数，n为步数，k为图片保存的序列
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    求解不加噪过程中的解包裹相位    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        noisy=0;%是否加噪，加噪为1，不加噪为0
        
        
        deltphi={};
        for m=1:fre
            result= 2*pi*D*gaussian_surfaceH./(L-gaussian_surfaceH)/T(m); % T为正弦条纹的空间周期, 这里的./为点除
            deltphi{end+1} =result;
            u1=u(m);
            delta_phi=deltphi{m};
           for n=1:step
            phi=2*pi*(n-4)/step;
            high_fre_grating=fringeModulation(u1,phi,delta_phi,height,width,noisy,noisy_value_guassian,k_h,k_value);
            path_unwrapped_no_noisy=['E:\liyimingPCL\博士课题\实验记录-源代码+过程\3.14-基本工件几何体仿真\code\new\test_no_noisy\',num2str(m), '_',num2str(n),'.mat'];
            save(path_unwrapped_no_noisy, 'high_fre_grating');
        
           end
        end
        
        
        % figure(1);imshow(First);
        path_in='E:\liyimingPCL\博士课题\实验记录-源代码+过程\3.14-基本工件几何体仿真\code\new\test_no_noisy\';
        path_out_wrapped_low1='E:\liyimingPCL\博士课题\实验记录-源代码+过程\3.14-基本工件几何体仿真\code\new\phi_wrapped_low_no_noisy\';
        path_out_wrapped_high1='E:\liyimingPCL\博士课题\实验记录-源代码+过程\3.14-基本工件几何体仿真\code\new\phi_wrapped_high_no_noisy\';
        path_out_unwrapped1='E:\liyimingPCL\博士课题\实验记录-源代码+过程\3.14-基本工件几何体仿真\code\new\phi_unwrapped_no_noisy\';
        
        [phi_wrapped_low_no_noisy,phi_wrapped_high_no_noisy,phi_unwrapped_no_noisy] = unwrap_mul_fre(4,6,k,path_in,path_out_wrapped_low1,path_out_wrapped_high1,path_out_unwrapped1) ;
        %返回低频/高频包裹条纹和解包裹条纹

       end
 end
 
  
disp('gaussian_surface的部分已经生成完毕');
