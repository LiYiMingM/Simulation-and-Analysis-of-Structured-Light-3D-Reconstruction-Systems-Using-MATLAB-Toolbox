function output = randon_matrix(height_max,scale1,scale2)
%%%生成随机矩阵,noise_max=0表示噪音的最大值为0,真实值的最小值为height_min，最大值为height_max+height_min
% parameters

scale=[scale1,scale2];
%height_min = 10;height_max = 40;noise_max=0;
size_min = 5; % minimum size of initial matrix
size_max = 10; % maximum size of initial matrix

  
% initial absolute phase
size_xy = randi([size_min,size_max]); % get size of initial matrix
initial_matrix = rand(size_xy,size_xy); % get initial matrix#每个元素的值都是0-1
%         figure;
%         mesh(initial_matrix);
%        view(-23,71);
%        title('随机小矩阵')
%        colorbar;
scale_pre =[scale(2) * 1.25,scale(2) * 1.25]; % get size of pre-enlarged absolute phase
pre_enlarged_ap = imresize(initial_matrix,scale_pre);  % get pre-enlarged absolute phase
%         figure;
%         mesh(pre_enlarged_ap);
%        view(-23,71);
%        title('扩充后的矩阵')
%        colorbar;

%%%%裁剪扩充后的矩阵
% 计算中心位置
center_row = ceil(scale_pre(1) / 2);
center_col = ceil(scale_pre(2) / 2);

% 计算裁剪区域的起始和结束索引
start_row = center_row - floor(scale1 / 2);
end_row = start_row + scale1 - 1;

start_col = center_col - floor(scale2 / 2);
end_col = start_col + scale2  - 1;    

initial_ap = pre_enlarged_ap(start_row:end_row, start_col:end_col);
%         figure;
%         mesh(initial_ap);
%        view(-23,71);
%        title('裁剪后的矩阵')
%        colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%最大值的范围在30-60之间%%%%%%%%%%%%%%%
ymax=unifrnd(10,height_max, 1, 1);    %%生成均匀分布的随机数   

%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%   

ap=zeros(scale(1),scale(2));
for i = 1:scale(1)
    for j = 1:scale(2)
     ap(i,j) = (ymax)*(initial_ap(i,j));%高度重新映射
    end
end
%         figure;
%         mesh(ap);
%        view(-23,71);
%        title('高度重映射后的矩阵')
%        colorbar;

gt = single(ap);%single类型可以减少内存提升计算速度

output=gt;

end