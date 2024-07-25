function output = yuanzhui(rl1,x1,y1,h1,scale1,scale2)
% 圆环绘制函数,r为半球半径,x,y为球心坐标,scale图形的长宽尺度
temp = zeros(scale1,scale2);  
random_number = randi([1, 4]);
for obj = 1:random_number
rl = rl1(randi(length(rl1)));
x = x1(randi(length(x1)));
y = y1(randi(length(y1)));
h = h1(randi(length(h1)));

    for m = 1:scale1
        for n = 1:scale2
          if(((m-x)^2+(n-y)^2)<=rl^2)
              temp(m,n) = h*(rl-sqrt((m-x)^2 + (n-y)^2))/rl;
       
          end
        end
    end
temp=single(temp);
output = temp;
end
end