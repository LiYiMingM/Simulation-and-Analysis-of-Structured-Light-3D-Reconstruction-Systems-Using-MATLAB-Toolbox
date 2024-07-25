function mask = generate_mask_with_black_circles(height, width)
    % 生成一个全为1的掩膜
    mask = ones(height, width);
    
    % 随机生成1到6个圆
    num_circles = randi([1, 10]);
    
    for n = 1:num_circles
        % 随机生成圆心位置
        center_x = randi([1, width]);
        center_y = randi([1, height]);
        
        % 随机生成半径，范围为5到10个像素
        radius = randi([5, 30]);
        
        % 在掩膜上绘制黑色圆形区域
        for x = 1:width
            for y = 1:height
                if (x - center_x)^2 + (y - center_y)^2 <= radius^2
                    mask(y, x) = 0; % 设置为黑色
                end
            end
        end
    end
end


