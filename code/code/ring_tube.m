function output = ring_tube(rl1,rs1,x1,y1,scale1,scale2)

% 半球绘制函数,r为半球半径,x,y为球心坐标,scale图形的长宽尺度

num_objects = randi([1, 5]);
temp = zeros(scale1,scale2);
for obj = 1:num_objects
x = x1(randi(length(x1)));
y = y1(randi(length(y1)));

rs = rs1(randi(length(rs1)));
rl = rl1(randi(length(rl1)));
    for m = 1:scale1
       for n = 1:scale2-20
           if((m-x)^2 + (n-y)^2 < rl^2 && (m-x)^2 + (n-y)^2 > rs^2)
                temp(m,n) = sqrt(rl^2 - (x-m)^2 - (y-n)^2);
            end
        end
       
    end
end

temp=single(temp);
output = temp;
% figure 
% mesh(output)
end