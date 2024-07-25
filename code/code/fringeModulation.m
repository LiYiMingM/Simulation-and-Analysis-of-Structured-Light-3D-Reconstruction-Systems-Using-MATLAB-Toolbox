% 物体产生相移条纹函数(x轴频率/y轴频率，原条纹初相，因被测物体，高度调制的相移量，尺度,是否有噪声，噪声值)


% function [grating_left] = fringeModulation(fre,phi,delta_phi,height,width,noisy,k_h,k_value) 
% 
% 
%  if noisy==1
% % amp_harm = [k_value/1, k_value/2, k_value/5, k_value/10, k_value/20,k_value/100]; % 四次谐波振幅
% % freq_harm = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]; % 四次谐波频率
% % phase_harm = rand(size(amp_harm))*2*pi; % 四次谐波随机相位
% amp_harm = [k_value*2, k_value*0.01, k_value*0.04, k_value*0.02, k_value*0.01,k_value*0.01]; % 四次谐波振幅
% 
%         %%%%%%%%%左投影光栅%%%%%%%%%%%
% 
%              img_left = 0.5*ones(height,width);
%                 for a = 1:width
%                     for b = 1:height
%                         for i=1:k_h
%                         img_left(b,a) = img_left(b,a)+amp_harm(i)*cos(i*(fre*(a)+phi-delta_phi(b,a)));%amp_harm(i) 是第 i 个谐波的振幅
%             
%                         end
%                     end
%                 end
%            
%        
%         
% 
%  else
% 
% %加入高斯噪声和高次谐波
% % A=1;
%             
%             %%%%%%%无噪音的输出
%             %%%%%%%%%%%%左投影%%%%%%%%%%%%
%                  img_left = zeros(height,width);
%                   for a = 1:width
%                         for b = 1:height
%                          img_left(b,a) = (1+cos(fre*(a)+phi-delta_phi(b,a)))/2;
%                         end
%                   end 
% 
%    
%  end
% 
% 
% 
% 
% 
% grating_left=single(img_left);
% 
% end



function [grating_left,img_left_nonoisy] = fringeModulation(fre,phi,delta_phi,height,width,noisy,k_h,k_value) 


 if noisy==1
% amp_harm = [k_value/1, k_value/2, k_value/5, k_value/10, k_value/20,k_value/100]; % 四次谐波振幅
% freq_harm = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]; % 四次谐波频率
% phase_harm = rand(size(amp_harm))*2*pi; % 四次谐波随机相位
amp_harm = [0.5, k_value*1, k_value*0.4, k_value*0.2, k_value*0.1,k_value*0.05]; % 四次谐波振幅

        %%%%%%%%%左投影光栅%%%%%%%%%%%

             img_left = 0.5*ones(height,width);
                for a = 1:width
                    for b = 1:height
                        for i=1:k_h
                        img_left(b,a) = img_left(b,a)+amp_harm(i)*cos(i*(fre*(a)+phi-delta_phi(b,a)));%amp_harm(i) 是第 i 个谐波的振幅
            
                        end
                    end
                end
           
       
        

%加入高斯噪声和高次谐波
% A=1;
            
            %%%%%%%无噪音的输出
            %%%%%%%%%%%%左投影%%%%%%%%%%%%
                 img_left_nonoisy = zeros(height,width);
                  for a = 1:width
                        for b = 1:height
                         img_left_nonoisy(b,a) = (1+cos(fre*(a)+phi-delta_phi(b,a)))/2;
                        end
                  end 

   
 end



img_left_nonoisy=single(img_left_nonoisy);

grating_left=single(img_left);

end









