
function output = ladderFun(px_a,py_a,stepNum_a,step_a,width_a,scale )
% ���ݻ��ƺ���,x,yָ��̨�׷���λ��step����,stepNum̨����,width̨�׵Ŀ�ȣ�scaleͼ�εĳ���߶�

 %width ̨�׵��ܿ��
%length ̨�׵��ܳ��ȣ�
%step,һ��̨�׵ø߶�

% scale = 256;
% px_a=(20):(40);%���ŵ���x�ᣬ�����x�����㡣
% py_a=(20):(25);%���ŵ���y�ᣬ����y������
% step_a=5:1:7;%����
% stepNum_a=3:8;%�����������ŵ���y��
% width_a=100:200;%x��ķ�Χ

px = px_a(randi(length(px_a)));
py = py_a(randi(length(py_a)));
temp = zeros(scale(1),scale(2));

%���ֽ�����״ÿ�����ѡ��һ��
random_number = randi([1, 2]);
if random_number == 1
%����ƽ�������ƽ��
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

%ѡ�����ɽ���ƽ������ƴ�ֱ����״
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


