function [CylData,BranchData] = QSM_struct2array(QSM)

           C = QSM.cylinder;
           format long
           CylData = cat(2,C.radius ,C.length, C.start, C.axis, single(C.parent), single(C.extension), single(C.branch) ,...
               single(C.BranchOrder), single(C.PositionInBranch), C.added);
           B = QSM.branch;
           BranchData = cat(2,B.order, B.parent, B.volume, B.length, B.angle);

      %CylData = single([Rad Len Sta Axe CPar CExt BoC Added]);
       % BranchData = single([BOrd BPar BVol BLen BAng]);

end

