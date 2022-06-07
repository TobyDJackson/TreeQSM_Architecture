function [QSM] = QSM_bits2struct(Rad,Len,Sta,Axe,CPar,CExt,BoC,Added,BOrd,BPar,BVol,BLen,BAng)
%This function takes the pieces that come out of the QSM sometimes and
%saves them as a struct.

     cylinder=struct('radius',Rad,'length',Len,'start',Sta,'axis',Axe,'parent',CPar,'extension',CExt,...
         'BranchOrder',BoC(:,2),'branch',BoC(:,1),'PositionInBranch',BoC(:,3),'added',Added);
     
     bases=find(cylinder.PositionInBranch==1);
     BDiam=cylinder.radius(bases)*2;
        branch=struct('order',BOrd,'parent',BPar,'volume',BVol,'length',BLen,'angle',BAng,'diameter',BDiam);
        QSM=struct('cylinder',cylinder,'branch',branch);
        
end

