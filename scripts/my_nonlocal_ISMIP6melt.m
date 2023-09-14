
function output = my_nonlocal_ISMIP6melt(X,Y,string),

% switch (oshostname()),
% 	case {'ronne'}
% 		rignotmelt='/home/ModelData/Antarctica/RignotMeltingrate/Ant_MeltingRate.nc';
% 	case {'totten'}
% 		rignotmelt='/totten_1/ModelData/Antarctica/RignotMeltingrate/Ant_MeltingRate.nc';
% 	case {'thwaites','murdo','astrid'}
% 		rignotmelt=['/home/seroussi/Data/Ant_MeltingRate.nc'];
% 	otherwise
% 		error('hostname not supported yet');
% end
rignotmelt='/Users/jbec0008/SAEF/datasets/basalmetrate/melt_rates_NON_LOCAL.nc'

if nargin==2,
	string = 'melt_rate';
end

disp(['   -- ISMIP6 local melt rate: loading ' string]);
xdata = double(ncread(rignotmelt,'x'));
ydata = double(ncread(rignotmelt,'y'));

disp(['   -- ISMIP6 local melt rate: loading' string]);
data  = double(ncread(rignotmelt,string))';

disp(['   -- ISMIP6 local melt rate interpolating' string]);
output = InterpFromGrid(xdata,ydata,data,X(:),Y(:));
output = reshape(output,size(X,1),size(X,2));
