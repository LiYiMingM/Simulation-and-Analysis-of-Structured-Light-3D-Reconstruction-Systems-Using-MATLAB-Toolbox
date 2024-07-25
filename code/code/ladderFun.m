
function output = ladderFun(px_a,py_a,stepNum_a,step_a,width_a,scale )
% 阶梯绘制函数,x,y指定台阶放置位置step步长,stepNum台阶数,width台阶的宽度，scale图形的长宽尺度

 %width 台阶的总宽度
%length 台阶得总长度，
%step,一个台阶得高度

% scale = 256;
% px_a=(20):(40);%竖着的是x轴，这代表x轴的起点。
% py_a=(20):(25);%横着的事y轴，这是y轴的起点
% step_a=5:1:7;%步长
% stepNum_a=3:8;%阶梯数，横着的是y轴
% width_a=100:200;%x轴的范围

px = px_a(randi(length(px_a)));
py = py_a(randi(length(py_a)));
temp = zeros(scale(1),scale(2));

%两种阶梯形状每次随机选择一个
random_number = randi([1, 2]);
if random_number == 1
%阶梯平面和条纹平行
        width = width_a(randi(length(width_a)));
        step = step_a(randi(length(step_a)));
        stepNum = stepNum_a(randi(length(stepNum_a)));
            for x = px:(px+width)
                 for y = (py):(py+step*stepNum)   
                    k=floor((y-1-(py))/step+1);
                   temp(x,y)=step*k;
                 end 
            end
        width = width_a(randi(length(width_a)));
        step = step_a(randi(length(step_a)));
        stepNum = stepNum_a(randi(length(stepNum_a)));   
             for x = px:(px+width)
                for y = (py+185):(py+step*stepNum+185)    
                  k=floor((y-1-(py+185))/step+1);
                   temp(x,y)=step*k;
                             
                end
             end
        
        step = step_a(randi(length(step_a)));
        stepNum = stepNum_a(randi(length(stepNum_a)));
        width = width_a(randi(length(width_a)));
              for x = px:(px+width)
                for y = (py+244):(py+step*stepNum+244)    
                  k=floor((y-1-(py+244))/step+1);
                   temp(x,y)=step*k;
               end
              end

%选择生成阶梯平面和条纹垂直的形状
else

        width = width_a(randi(length(width_a)));
        step = step_a(randi(length(step_a)));
        stepNum = stepNum_a(randi(length(stepNum_a)));
            for  y= py:(py+width)
                 for x = (px):(px+step*stepNum)   
                    k=floor((x-1-(px))/step+1);
                   temp(x,y)=step*k;
                 end 
            end
        width = width_a(randi(length(width_a)));
        step = step_a(randi(length(step_a)));
        stepNum = stepNum_a(randi(length(stepNum_a)));   
             for y = py:(py+width)
                for x = (px+185):(px+step*stepNum+185)    
                  k=floor((x-1-(px+185))/step+1);
                   temp(x,y)=step*k;
                             
                end
             end
        
        step = step_a(randi(length(step_a)));
        stepNum = stepNum_a(randi(length(stepNum_a)));
        width = width_a(randi(length(width_a)));
              for y = py:(py+width)
                for x = (px+244):(px+step*stepNum+244)    
                  k=floor((x-1-(px+244))/step+1);
                   temp(x,y)=step*k;
               end
              end






end
temp=single(temp);
output = temp;
end


