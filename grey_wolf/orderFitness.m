function orderedWolfFitness = orderFitness(subWolfPosition,fobj)

no = size(subWolfPosition,1);


%对每一行进行求和
%第二个参数代表2 对每一行进行适应度函数计算
%存储 适应度值和序号信息
fitnessinfo = zeros(no,2);
%根据某个适应度函数值进行排序
for index=1:no     
   %计算每匹狼的适应度函数
   fitnessindex =  fobj(subWolfPosition(index,:));
   
   %赋值到适应度和序号信息
   fitnessinfo(index,1) =  fitnessindex; 
   fitnessinfo(index,2) = index;  
   
end

 fitnessinfo = sortrows(fitnessinfo,1);
 
 orderedWolfFitness = fitnessinfo;

end