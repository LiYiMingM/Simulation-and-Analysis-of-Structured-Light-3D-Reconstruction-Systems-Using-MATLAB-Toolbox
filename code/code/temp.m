clc
clear all
clf

scale = 256;
X = [1:1:1*scale];Y = [1:1:1*scale];
[x,y] = meshgrid(X,Y);
height_all_max=60;
g_num_k = 15;%%高斯函数的峰值个数。
%gaussian_surfaceH=0;
mux_a=(20):(236);%高斯曲面的中心位置
muy_a=(20):(236);
sigmax_a=10:40;%高斯曲面的范围
sigmay_a=10:40;
amplitude_a=10:20;
gaussian_surface_k=3;
k=0;
 for i=1:gaussian_surface_k
     
   %  for i=1:gaussian_surface_k
       k=k+1;
        g_sH=0;
        for i=1:g_num_k
         g_s=gaussian_surface(mux_a,muy_a,sigmax_a,sigmay_a,amplitude_a,scale );
         g_sH=g_sH+g_s;
        end
      for i = 1:size(g_sH,1)
       for j = 1:size(g_sH,2)
        
        % 如果元素的值大于60
        if g_sH(i,j) > 60
            
            % 生成50-60之间的随机数并替换该元素
            g_sH(i,j) = randi([40 60]);
            
        end
        
       end
  end
      geometry=g_sH/(height_all_max+5);

figure
imshow(g_sH,[])
figure
mesh(x, y, g_sH);
xlabel('x')
ylabel('y')
zlabel('z')
title('高斯曲面')
figure
mesh(x, y, geometry);
xlabel('x')
ylabel('y')
zlabel('z')
title('高斯曲面')
 end