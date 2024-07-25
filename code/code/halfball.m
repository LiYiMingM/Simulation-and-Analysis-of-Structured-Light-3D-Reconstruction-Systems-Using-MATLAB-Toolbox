function output = halfball(ball_r,ball_x1,ball_y1,scale)
% ������ƺ���,rΪ����뾶,x,yΪ��������,scaleͼ�εĳ���߶�
%һ��ͼ�������ĸ���С��һ�����壬����˳ʱ���ṩ�ĸ����ĵ�����(x1,y1),
% (x2,y2),(x3,y3),(x4,y4)
%�Ҳ����뾶/2��Ϊ�˷�ֹ�Ҳ೬��Χ


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
