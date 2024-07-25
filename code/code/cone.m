function output = cone(rl,rs,x,y,h,scale1,scale2)
% 锥台，rl是外环半径，rs是内环半径，x,y是圆心的坐标，h是锥台的高度，


temp = zeros(scale1, scale2);  % 初始化temp矩阵
num_objects=randi([1, 4]);

for obj = 1:num_objects
    % 随机生成每个物体的参数
      rs1 = rs(randi(length(rs)));
      rl1 = rl(randi(length(rl)));
      x1 = x(randi(length(x)));
      y1 = y(randi(length(y)));
      h1 = h(randi(length(h)));

    % 处理每个物体的参数
    for m = 1:scale1
        for n = 1:scale2
            if (rs1^2 < ((m-x1)^2 + (n-y1)^2) && ((m-x1)^2+(n-y1)^2) < rl1^2)
                temp(m,n) = h1 * (rl1 - sqrt((m-x1)^2 + (n-y1)^2)) / rl1;
            elseif (((m-x1)^2 + (n-y1)^2) <= rs1^2)
                temp(m,n) = h1 - h1 * rs1 / rl1;
            end
        end
    end
end

temp=single(temp);
output = temp;
end