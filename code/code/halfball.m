function output = halfball(ball_r,ball_x1,ball_y1,scale)
% 半球绘制函数,r为半球半径,x,y为球心坐标,scale图形的长宽尺度
%一张图中生成四个大小不一的球体，按照顺时针提供四个球心的坐标(x1,y1),
% (x2,y2),(x3,y3),(x4,y4)
%右侧的球半径/2是为了防止右侧超范围


num_object = randi([1, 6]);
temp = zeros(scale(1),scale(2));



for i=1:num_object
   r = ball_r(randi(length(ball_r)));
   x1 = ball_x1(randi(length(ball_x1)));
   y1 = ball_y1(randi(length(ball_y1)));

    for mx = 1:scale(1)
        for ny = 1:scale(2)
            if((mx-x1)^2 + (ny-y1)^2 < r^2)
                temp(mx,ny) = sqrt(r^2 - (x1-mx)^2 - (y1-ny)^2);
            end
        end
    end


temp=single(temp);
output = temp;

end
% figure 
% mesh(output)
% max(max(output))
