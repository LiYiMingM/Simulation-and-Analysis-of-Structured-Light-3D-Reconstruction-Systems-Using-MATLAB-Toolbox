function output = gaussian_surface(mux_a,muy_a,sigmax_a,sigmay_a,amplitude_a,scale1,scale2 )
%这个高斯曲面由四个唯一峰值的高斯曲面相加而得，视觉效果上就是一个多峰的高斯曲面。
%每个高峰值都有5个参数，分别是mu_x，代表高斯函数在x轴的中心位置
%mu_y第一个高斯函数在y轴的中心位置
%sigma_x高斯函数在x轴方向上的标准差，用于控制函数在x轴方向上的宽度
%sigma_yy轴方向上的标准差，用于控制函数在y轴方向上的宽度
%amplitude高斯函数的振幅
      mu_x = mux_a(randi(length(mux_a)));
      mu_y = muy_a(randi(length(muy_a)));
      sigma_x = sigmax_a(randi(length(sigmax_a)));
      sigma_y = sigmay_a(randi(length(sigmay_a)));
      amplitude = amplitude_a(randi(length(amplitude_a)));

f1 = @(x, y) amplitude * exp(-((x-mu_x).^2/(2*sigma_x^2) + (y-mu_y).^2/(2*sigma_y^2)));
        
% 在范围 [0, 512] 中创建网格点


temp1 = zeros(scale1,scale2);

for x = 1:scale1
    for y = 1:scale2-20
     temp1(x,y)=f1(x,y);
   end
 end
% 使用 mesh 函数绘制高斯曲面




temp1=single(temp1);
output = temp1;
end

