function orderedWolfFitness = orderFitness(subWolfPosition,fobj)

no = size(subWolfPosition,1);


%��ÿһ�н������
%�ڶ�����������2 ��ÿһ�н�����Ӧ�Ⱥ�������
%�洢 ��Ӧ��ֵ�������Ϣ
fitnessinfo = zeros(no,2);
%����ĳ����Ӧ�Ⱥ���ֵ��������
for index=1:no     
   %����ÿƥ�ǵ���Ӧ�Ⱥ���
   fitnessindex =  fobj(subWolfPosition(index,:));
   
   %��ֵ����Ӧ�Ⱥ������Ϣ
   fitnessinfo(index,1) =  fitnessindex; 
   fitnessinfo(index,2) = index;  
   
end

 fitnessinfo = sortrows(fitnessinfo,1);
 
 orderedWolfFitness = fitnessinfo;

end