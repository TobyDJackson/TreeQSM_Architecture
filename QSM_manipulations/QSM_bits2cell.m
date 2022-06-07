function [CylData, BranchData] = QSM_bits2cell(Rad,Len,Sta,Axe,CPar,CExt,BoC,Added,BOrd,BPar,BVol,BLen,BAng)
%This function takes the pieces that come out of the QSM sometimes and
%saves them as a cell
        CylData = single([Rad Len Sta Axe CPar CExt BoC Added]);  
        BranchData = single([BOrd BPar BVol BLen BAng]);
        %ModelData={CylData, BranchData}; %repackages data in cell array for simplify code   

end
