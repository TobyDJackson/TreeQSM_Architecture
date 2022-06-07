function [QSM] = QSM_cell2struct(CylData,BranchData)
    %This function takes the CylData and BranchData that come out of the QSM sometimes and
    %saves them as a struct.
   
    
    Rad=CylData(:,1); Len=CylData(:,2); Sta=CylData(:,3:5); Axe=CylData(:,6:8); CPar=CylData(:,9); CExt=CylData(:,10); BoC=CylData(:,11:13); Added=CylData(:,14); 
    BOrd=BranchData(:,1); BPar=BranchData(:,2); BVol=BranchData(:,3); BLen=BranchData(:,4); BAng=BranchData(:,5); 
    
    cylinder=struct('radius',Rad,'length',Len,'start',Sta,'axis',Axe,'parent',CPar,'extension',CExt,...
         'BranchOrder',BoC(:,2),'branch',BoC(:,1),'PositionInBranch',BoC(:,3),'added',Added);
     
     bases=find(cylinder.PositionInBranch==1);
     BDiam=cylinder.radius(bases)*2;
     
     branch=struct('order',BOrd,'parent',BPar,'volume',BVol,'length',BLen,'angle',BAng,'diameter',BDiam);
     
     QSM=struct('cylinder',cylinder,'branch',branch);
        
end

