function [ CylData_simplified_added ] = add_simplified_area( CylData, CylData_updated )
% This functions is intended to simulate the effect of leaves by adding area to the outer cylinders

[ Arch_Values_original Sail Arch_Names ] = Calculate_architectures_2017( CylData );
Sail_Original=Arch_Values_original(10);
Tot_Original=sum(CylData(:,1).*CylData(:,2));

[ Arch_Values_simplified Sail Arch_Names ] = Calculate_architectures_2017( CylData_updated  );
Sail_simplified=Arch_Values_simplified(10)
Tot_simplified=sum(CylData_updated(:,1).*CylData_updated(:,2));

[length(CylData) Sail_Original Tot_Original length(CylData_updated) Sail_simplified Tot_simplified]
CylData_simplified_added=[length(CylData) Sail_Original Tot_Original length(CylData_updated) Sail_simplified Tot_simplified];

end

