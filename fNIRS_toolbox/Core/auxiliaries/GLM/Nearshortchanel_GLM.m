function [channelspershortch ] = Nearshortchanel_GLM( lstss,SD )
% lstss is the list of short channels in the probe
% To give the channels closer to  all short channel in the probe in ceils

ml=SD.MeasList;
ml=ml(find(ml(:,4)==1),:);
PosDetshort=zeros(length(lstss),3);
PosCMchannel=zeros(size(ml,1),3);
dist=zeros(size(ml,1),length(lstss));

for i=1:size(ml,1)
    
    if isfield(SD,'DetPos_3d')
        PosCMchannel(i,:)=(SD.SrcPos_3d(ml(i,1),:)+SD.DetPos_3d(ml(i,2),:))/2;
    else
        PosCMchannel(i,:)=(SD.SrcPos(ml(i,1),:)+SD.DetPos(ml(i,2),:))/2; 
    end
end

for i=1:length(lstss)
    
    if isfield(SD,'DetPos_3d')
        PosDetshort(i,:)=SD.DetPos_3d(ml(lstss(i),2),:);
    else
        PosDetshort(i,:)=SD.DetPos(ml(lstss(i),2),:);
    end
    
    for j=1:size(ml,1)
        
        dist(j,i)=sqrt((PosDetshort(i,1)-PosCMchannel(j,1)).^2+(PosDetshort(i,2)-PosCMchannel(j,2)).^2+(PosDetshort(i,3)-PosCMchannel(j,3)).^2);
        
        
    end

end
  for i=1:length(lstss)
          k=0;
          DetsperSC=[];
  for j=1:size(ml,1)
      lst2 =  min(find(dist(j,:)==min(dist(j,:))));
      if i==lst2
          k=k+1;
          DetsperSC(k)=j;
      end
     
   end
     
     
     channelspershortch{i}=sort(DetsperSC);
  end
 
 end
  







