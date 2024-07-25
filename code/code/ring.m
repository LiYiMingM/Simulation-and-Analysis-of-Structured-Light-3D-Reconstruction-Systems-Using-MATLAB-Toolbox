function output = ring(rl,rs,x,y,h,scale,height_all_max)
% 圆环绘制函数,rl为半球的大环半径,rs为半球的小环半径
% x,y为球心坐标,h为圆环的高度，scale图形的长宽尺度
temp = zeros(scale);
num_circles = randi([1, 3]);
for i = 1:num_circles
    % 随机选择圆环的参数
    h_val = h(randi(length(h)));
    rs_val = rs(randi(length(rs)));
    rl_val = rl(randi(length(rl)));
    x_val = x(randi(length(x)));
    y_val = y(randi(length(y)));
    
    % 在矩阵中生成圆环
    for m = 1:scale
        for n = 1:scale
            if rs_val^2 < (m - x_val)^2 + (n - y_val)^2 && (m - x_val)^2 + (n - y_val)^2 < rl_val^2
                temp(m, n) = h_val;
            end
        end
    end
end


output = temp;
output(output > height_all_max) = height_all_max;

output=single(output);
end