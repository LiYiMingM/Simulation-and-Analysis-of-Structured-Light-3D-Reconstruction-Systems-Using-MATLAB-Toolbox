function [phi_wrapped_low,phi_wrapped_middle,phi_wrapped_high,phi_unwrapped,kh] = unwrap_mul_fre(m,n,k,path_in,path_out_wrapped_low,path_out_wrapped_middle,path_out_wrapped_high,path_out_unwrapped) 
%%%%%%%%m为频率数，n为步数，k为图片保存的序列，path1为'E:\liyimingPCL\博士课题\实验记录-源代码+过程\3.14-基本工件几何体仿真\code\new\test\'
% % m=4;  %频率数
% % n=6;  %相移数

HIGHrLOW = [4,8];%m-1频率之比1 4 16 32
Img_total=cell(m,n);


%%%数据读取

    for i=1:m
        for j=1:n
            path=[path_in, num2str(i), '_', num2str(j),'.mat'];
            Img_total{i,j}=load(path);
            fieldname = fieldnames(Img_total{i,j});
            Img_total{i,j}=im2double(Img_total{i,j}.(fieldname{1}));
        end
    end
    %%%提取包裹相位
    for i=1:m
        numerator=0;
        denominator=0;
        for j=1:n
            numerator=numerator+Img_total{i,j}*sin(2*(j-1)*pi/n);
            denominator=denominator+Img_total{i,j}*cos(2*(j-1)*pi/n);
        end
        phi(:,:,i)=-atan2(numerator,denominator);
    end
      phi_wrapped = phi; 
      phi_wrapped_low_nogy=phi_wrapped(:,:,1);
      phi_wrapped_middle_nogy=phi_wrapped(:,:,2);     
      phi_wrapped_high_nogy=phi_wrapped(:,:,3);
      %%%%%%对包裹相位进行归一化
      phi_wrapped_low=(phi_wrapped_low_nogy+pi)/(2*pi);
      phi_wrapped_middle=(phi_wrapped_middle_nogy+pi)/(2*pi);
      phi_wrapped_high=(phi_wrapped_high_nogy+pi)/(2*pi);
      phi_wrapped_low=single(phi_wrapped_low);
      phi_wrapped_middle=single(phi_wrapped_middle);     
      phi_wrapped_high=single(phi_wrapped_high);
      
      path1=[path_out_wrapped_low,num2str(k,'%06d'),'.mat'];  
      save(path1, 'phi_wrapped_low');
      path2=[path_out_wrapped_middle,num2str(k,'%06d'),'.mat'];  
      save(path2, 'phi_wrapped_middle');
      path3=[path_out_wrapped_high,num2str(k,'%06d'),'.mat'];  
      save(path3, 'phi_wrapped_high');
    
    %%相位展开
    phl = phi(:,:,1)+pi;
    for i=1:m-1
        phh = phi(:,:,i+1)+pi;
        kh = round((HIGHrLOW(i)*phl-phh)/(2*pi));
        khlist(:,:,i) = kh;
        phl = phh + kh*2*pi;
        phllist(:,:,i) = phl;
    end
      phi_unwrapped_nogy = phl; 
      %%%%%对解包裹相位进行归一化
      phi_unwrapped=phi_unwrapped_nogy/(pi*(1+2*(HIGHrLOW(1)*HIGHrLOW(2))+5));%%归一化除的是（pi+2k+1）
      phi_unwrapped=single(phi_unwrapped);%提高运算速度
      path4=[path_out_unwrapped,num2str(k,'%06d'),'.mat'];  
      save(path4, 'phi_unwrapped');
      
      kh=single(kh);

end



