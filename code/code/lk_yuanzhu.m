function output = lk_yuanzhu(l1,r1,rs1,scale1,scale2)

temp = zeros(scale1,scale2);
num_object=randi([1,2]);

for i=1:num_object

l = l1(randi(length(l1)));
r = r1(randi(length(r1)));
rs = rs1(randi(length(rs1)));
x=randi([1, scale1/2+20]);%半圆柱的x起点
y=randi([10, scale2/2]);%半圆柱的左侧y起点
xl=randi([x+10,x+2*r-10]);%镂空圆柱的x起点
yl=randi([y+10,y+l-10]);%镂空圆柱的y起点

 for m = x:x+2*r
      for n =y:(y+l)
           if m<=x+r
                 temp(m,n) = sqrt((r^2 - (r+x-m)^2));
                  if(xl)<m && m<(xl+rs) && (yl)<n && n<(yl+rs)
                  temp(m,n)=0;
                  end

           else
                     temp(m,n) = sqrt((r^2 - (m-x-r)^2));
                  if(xl)<m && m<(xl+rs) && (yl)<n && n<(yl+2*rs)
                  temp(m,n)=0;
                  end

           end
      end


end
temp=single(temp);
output = temp;
% figure
% mesh(output)

end