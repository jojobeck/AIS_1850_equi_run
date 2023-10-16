steps = [2];
steps = [12];
% steps =[9];

loadonly = 0;
loadonly = 1;
addpath('scripts');



% org=organizer('repository',['./Models'],'prefix',['ISMIP6Antarctica_'],'steps',steps, 'color', '34;47;2'); 
org=organizer('repository',['./Models'],'prefix',['AIS1850_'],'steps',steps, 'color', '34;47;2'); 
clear steps;
datadir= '/Users/jbec0008/SAEF/datasets/';
data_smb='/Volumes/Crucial X8/SAEF/ISMIP6/data_and_preprocessing/preprocess/SMB_JOHANNA/';
data_2km_xy='/Volumes/Crucial X8/SAEF/ISMIP6/data_and_preprocessing/published_data/ISMIP6/AtmosphericForcing/noresm1-m_rcp8.5/';
data_ocean='/Volumes/Crucial X8/SAEF/ISMIP6/data_and_preprocessing/published_data/ISMIP6/Ocean/';
data_tmp_ocean = '/Users/jbec0008/SAEF/issm_projects/equi_1850/Data/Ocean/';
data_tmp_smb = '/Users/jbec0008/SAEF/issm_projects/equi_1850/Data/Atmosphere/';
% loadonly=1;
%creeate forcing 1 -6
if perform(org,'1850-clim-AtmForcing')% {{{

	% md=loadmodel('./Models/ISMIP6Antarctica__CollapseSSA.math');
	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');
    % datapath='data_smb;

	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;
	
	md.inversion.iscontrol=0;
	md.transient.isthermal=0;
	md.transient.isgroundingline=1;
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);

	% disp('loading TS and SMB climatology data');
	smbclimnc_lonlat           = [data_smb 'UKESM1-0-LL_clim_ssp585_1995-2014_4km.nc'];
	lat                 = double(ncread(smbclimnc_lonlat,'lat'));
	lon                 = double(ncread(smbclimnc_lonlat,'lon'));
	% smb_clim_data       = double(ncread(smbclimnc,'smb_clim'));
	% % ts_clim_data        = double(ncread(smbclimnc,'ts_clim'));

	disp('loading TS and SMB anomoly data');
	smbanomnc           = [data_smb 'MAR-UKESM1-0-LL_asmb_1850-1969_mean_histo_regrid_04000m_EXTENDED.nc'];
	smb_anomaly_data    = double(ncread(smbanomnc,'asmb'));
	tsanomnc           = [data_smb 'tas_Ayr_mean_UKESM1-0-LL_historical_r1i1p1f2_1850_1869.nc'];
	ts_anomaly_data     = double(ncread(tsanomnc,'tas'));

	%Create SMB and TS matrix
	t=[1:size(smb_anomaly_data,3)];
	[x_n y_n]=ll2xy(lat(:,1),lon(:,1),-1);
	y_n = x_n;
	% smb_clim=InterpFromGridToMesh(x_n,y_n,smb_clim_data',md.mesh.x,md.mesh.y,0);
	% ts_clim=InterpFromGridToMesh(x_n,y_n,ts_clim_data',md.mesh.x,md.mesh.y,0);
	temp_matrix_smb = []; temp_matrix_ts = [];
    smb_clim=interpRacmoSMB_ISMIP6(md,datadir);
    ts_clim=interpRacmoTemp_ISMIP6(md,datadir);
	rhoi = md.materials.rho_ice;
	for i = 1:size(smb_anomaly_data,3)
		disp(i);
		%SMB
		temp_smb        = InterpFromGridToMesh(x_n,y_n,smb_anomaly_data(:,:,i)',md.mesh.x,md.mesh.y,0);
		temp_smb        = temp_smb*((31556926/1000)*(1000/rhoi));%m/year
		temp_smb        = temp_smb+smb_clim;
		temp_matrix_smb = [temp_matrix_smb temp_smb];
		%TS
		temp_ts         = InterpFromGridToMesh(x_n,y_n,ts_anomaly_data(:,:,i)',md.mesh.x,md.mesh.y,0);
		% temp_ts         = temp_smb+smb_clim;
		temp_ts         = temp_ts+ts_clim;
		temp_matrix_ts  = [temp_matrix_ts temp_ts];
		clear temp_smb; clear temp_ts;
	end

	%Save Data (1995-2100)
	md.smb.mass_balance       = temp_matrix_smb ;
	miroc_rcp85_smb           = md.smb.mass_balance;
	md.miscellaneous.dummy.ts = temp_matrix_ts;
	miroc_rcp85_ts           = md.miscellaneous.dummy.ts;
	save('Data/Atmosphere/ukesm_histo_clim_smb.mat','miroc_rcp85_smb');
	save('Data/Atmosphere/ukesm_histo_clim_ts.mat','miroc_rcp85_ts');

end %}}}
if perform(org,'1850-cold-clim-AtmForcing')% {{{

	% md=loadmodel('./Models/ISMIP6Antarctica__CollapseSSA.math');
	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');
    % datapath='data_smb;

	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;
	
	md.inversion.iscontrol=0;
	md.transient.isthermal=0;
	md.transient.isgroundingline=1;
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);

	% disp('loading TS and SMB climatology data');
	smbclimnc_lonlat           = [data_smb 'UKESM1-0-LL_clim_ssp585_1995-2014_4km.nc'];
	lat                 = double(ncread(smbclimnc_lonlat,'lat'));
	lon                 = double(ncread(smbclimnc_lonlat,'lon'));
	% smb_clim_data       = double(ncread(smbclimnc,'smb_clim'));
	% % ts_clim_data        = double(ncread(smbclimnc,'ts_clim'));

	disp('loading TS and SMB anomoly data');
	smbanomnc           = [data_tmp_smb 'MAR-UKESM1-0-LL_asmb_1953-1972_mean_histo_regrid_04000m_EXTENDED.nc'];
	smb_anomaly_data    = double(ncread(smbanomnc,'asmb'));
	% tsanomnc           = [data_smb 'tas_Ayr_mean_UKESM1-0-LL_historical_r1i1p1f2_1850_1869.nc'];
	% ts_anomaly_data     = double(ncread(tsanomnc,'tas'));

	%Create SMB and TS matrix
	t=[1:size(smb_anomaly_data,3)];
	[x_n y_n]=ll2xy(lat(:,1),lon(:,1),-1);
	y_n = x_n;
	% smb_clim=InterpFromGridToMesh(x_n,y_n,smb_clim_data',md.mesh.x,md.mesh.y,0);
	% ts_clim=InterpFromGridToMesh(x_n,y_n,ts_clim_data',md.mesh.x,md.mesh.y,0);
	temp_matrix_smb = []; temp_matrix_ts = [];
    smb_clim=interpRacmoSMB_ISMIP6(md,datadir);
    % ts_clim=interpRacmoTemp_ISMIP6(md,datadir);
	rhoi = md.materials.rho_ice;
	for i = 1:size(smb_anomaly_data,3)
		disp(i);
		%SMB
		temp_smb        = InterpFromGridToMesh(x_n,y_n,smb_anomaly_data(:,:,i)',md.mesh.x,md.mesh.y,0);
		temp_smb        = temp_smb*((31556926/1000)*(1000/rhoi));%m/year
		temp_smb        = temp_smb+smb_clim;
		temp_matrix_smb = [temp_matrix_smb temp_smb];
		%%TS
		%temp_ts         = InterpFromGridToMesh(x_n,y_n,ts_anomaly_data(:,:,i)',md.mesh.x,md.mesh.y,0);
		%% temp_ts         = temp_smb+smb_clim;
		%temp_ts         = temp_ts+ts_clim;
		%temp_matrix_ts  = [temp_matrix_ts temp_ts];
		% clear temp_smb; clear temp_ts;
		clear temp_smb; 
	end

	%Save Data (1995-2100)
	md.smb.mass_balance       = temp_matrix_smb ;
	miroc_rcp85_smb           = md.smb.mass_balance;
	% md.miscellaneous.dummy.ts = temp_matrix_ts;
	% miroc_rcp85_ts           = md.miscellaneous.dummy.ts;
	save('Data/Atmosphere/ukesm_cold_clim_smb.mat','miroc_rcp85_smb');
	% save('Data/Atmosphere/ukesm_histo_clim_ts.mat','miroc_rcp85_ts');

end %}}}
if perform(org,'1850-AtmForcing')% {{{

	% md=loadmodel('./Models/ISMIP6Antarctica__CollapseSSA.math');
	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');
    % datapath='data_smb;

	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;
	
	md.inversion.iscontrol=0;
	md.transient.isthermal=0;
	md.transient.isgroundingline=1;
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);

	% disp('loading TS and SMB climatology data');
	smbclimnc_lonlat           = [data_smb 'UKESM1-0-LL_clim_ssp585_1995-2014_4km.nc'];
	lat                 = double(ncread(smbclimnc_lonlat,'lat'));
	lon                 = double(ncread(smbclimnc_lonlat,'lon'));
	% smb_clim_data       = double(ncread(smbclimnc,'smb_clim'));
	% % ts_clim_data        = double(ncread(smbclimnc,'ts_clim'));

	disp('loading TS and SMB anomoly data');
	smbanomnc           = [data_smb 'MAR-UKESM1-0-LL_asmb_1850-1979_histo_regrid_04000m_EXTENDED.nc'];
	smb_anomaly_data    = double(ncread(smbanomnc,'asmb'));
	tsanomnc           = [data_smb 'tas_Ayr_UKESM1-0-LL_historical_r1i1p1f2_1850_2014.nc'];
	ts_anomaly_data     = double(ncread(tsanomnc,'tas'));

	%Create SMB and TS matrix
	t=[1:size(smb_anomaly_data,3)];
	[x_n y_n]=ll2xy(lat(:,1),lon(:,1),-1);
	y_n = x_n;
	% smb_clim=InterpFromGridToMesh(x_n,y_n,smb_clim_data',md.mesh.x,md.mesh.y,0);
	% ts_clim=InterpFromGridToMesh(x_n,y_n,ts_clim_data',md.mesh.x,md.mesh.y,0);
	temp_matrix_smb = []; temp_matrix_ts = [];
    smb_clim=interpRacmoSMB_ISMIP6(md,datadir);
    % ts_clim=interpRacmoTemp_ISMIP6(md,datadir);
	rhoi = md.materials.rho_ice;
	for i = 1:size(smb_anomaly_data,3)
		% disp(i);
		%SMB
		temp_smb        = InterpFromGridToMesh(x_n,y_n,smb_anomaly_data(:,:,i)',md.mesh.x,md.mesh.y,0);
		temp_smb        = temp_smb*((31556926/1000)*(1000/rhoi));%m/year
		temp_smb        = temp_smb+smb_clim;
		temp_matrix_smb = [temp_matrix_smb temp_smb];
		%TS
		% temp_ts         = InterpFromGridToMesh(x_n,y_n,ts_anomaly_data(:,:,i)',md.mesh.x,md.mesh.y,0);
		% % temp_ts         = temp_smb+smb_clim;
		% temp_ts         = temp_ts+ts_clim;
		% temp_matrix_ts  = [temp_matrix_ts temp_ts];
		clear temp_smb; clear temp_ts;
	end

	%Save Data (1995-2100)
	md.smb.mass_balance       = temp_matrix_smb ;
	miroc_rcp85_smb           = md.smb.mass_balance;
	% md.miscellaneous.dummy.ts = temp_matrix_ts;
	% miroc_rcp85_ts           = md.miscellaneous.dummy.ts;
	save('Data/Atmosphere/ukesm_histo_smb.mat','miroc_rcp85_smb');
	% save('Data/Atmosphere/ukesm_histo_clim_ts.mat','miroc_rcp85_ts');

end %}}}
if perform(org,'1850-clim-OcnForcing')% {{{
	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');

	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;
	
	md.inversion.iscontrol=0;
	md.transient.isthermal=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.isgroundingline=1;
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);

	%% Load and remesh basin_id data
    basin_datanc=[data_ocean 'imbie2/imbie2_basin_numbers_8km.nc'];
	x_n            = double(ncread(basin_datanc,'x'));
	y_n            = double(ncread(basin_datanc,'y'));
	basinid_data   = double(ncread(basin_datanc,'basinNumber'));
	basinid        = InterpFromGridToMesh(x_n,y_n,basinid_data',md.mesh.x,md.mesh.y,0)+1;
	basinid        = mean(basinid(md.mesh.elements),2); %per element

	%Remove interpolation errors from basinid data
	for i=1:16
		pos = find(basinid>i-0.5 & basinid < i+0.5);
		basinid(pos) = i;
	end
	save('Data/Ocean/basinid.mat','basinid');

	tfnc           = [data_ocean 'ukesm1-0-ll_ssp585/1850-1994/UKESM1-0-LL_thermal_forcing_8km_x_60m_1850-1896_mean.nc'];
	x_n            = double(ncread(tfnc,'x'));
	y_n            = double(ncread(tfnc,'y'));
	tf_data        = double(ncread(tfnc,'thermal_forcing'));
	z_data         = double(ncread(tfnc,'z'));

	%Build tf cell array
	t = 1:size(tf_data,4);
	tf = cell(1,1,size(tf_data,3));
	for i=1:size(tf_data,3)  %Iterate over depths
      temp_matrix=[];
		for ii=1:size(tf_data,4) %Iterate over time steps
			temp_tfdata=InterpFromGridToMesh(x_n,y_n,tf_data(:,:,i,ii)',md.mesh.x,md.mesh.y,0);
			temp_matrix = [temp_matrix temp_tfdata];
   	end
		temp_matrix = [temp_matrix ; t];
		tf(:,:,i)={temp_matrix};
	end
    obs_clim_tf = tf;
	save('Data/Ocean/ukesm_histo_clim_tf.mat','obs_clim_tf');

    %Set depths
    tf_depths = z_data';
    save('Data/Ocean/tf_depths.mat','tf_depths');
%
%	%Set ISMIP6 basal melt rate parameters
%	md.basalforcings            = basalforcingsismip6(md.basalforcings);
%	md.basalforcings.basin_id   = basinid;
%	md.basalforcings.num_basins = length(unique(basinid));
%	md.basalforcings.delta_t    = delta_t;
%	md.basalforcings.tf_depths  = tf_depths;
%	md.basalforcings.gamma_0    = 14477;
%	md.basalforcings.tf         = tf;
%
%	md.outputdefinition.definitions={};
%	md.timestepping.interp_forcings=1;
%	md.timestepping.final_time=106;
%	md.timestepping.time_step=1/12.;
%	%md.settings.output_frequency=24;
%	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','GroundedArea','FloatingArea','TotalSmb','SmbMassBalance','TotalGroundedBmb','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate',...
%		'IceVolumeScaled','IceVolumeAboveFloatationScaled','GroundedAreaScaled','FloatingAreaScaled','TotalSmbScaled','TotalGroundedBmbScaled','TotalFloatingBmbScaled','BasalforcingsIsmp6TfShelf'};
%
%	%Set melt / friction interpolation schemes
%	md.groundingline.migration = 'SubelementMigration';
%	md.groundingline.friction_interpolation='SubelementFriction1';
%	md.groundingline.melt_interpolation='SubelementMelt1';
%
%	%Solve
%	md.miscellaneous.name='ISMIP6_OCEAN_MELT_TEST';
%	md.verbose=verbose('solution',true,'module',true,'convergence',true);
%	md.cluster=cluster;
%	md.settings.waitonlock=0;
%	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
%
%	savemodel(org,md);
end %}}}
if perform(org,'1850-cold-clim-OcnForcing')% {{{
	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');

	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;
	
	md.inversion.iscontrol=0;
	md.transient.isthermal=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.isgroundingline=1;
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);

	%% Load and remesh basin_id data
    basin_datanc=[data_ocean 'imbie2/imbie2_basin_numbers_8km.nc'];
	x_n            = double(ncread(basin_datanc,'x'));
	y_n            = double(ncread(basin_datanc,'y'));
	basinid_data   = double(ncread(basin_datanc,'basinNumber'));
	basinid        = InterpFromGridToMesh(x_n,y_n,basinid_data',md.mesh.x,md.mesh.y,0)+1;
	basinid        = mean(basinid(md.mesh.elements),2); %per element

	%Remove interpolation errors from basinid data
	for i=1:16
		pos = find(basinid>i-0.5 & basinid < i+0.5);
		basinid(pos) = i;
	end
	save('Data/Ocean/basinid.mat','basinid');

	tfnc           = [data_tmp_ocean 'UKESM1-0-LL_thermal_forcing_8km_x_60m_1953-1972_mean.nc'];
	x_n            = double(ncread(tfnc,'x'));
	y_n            = double(ncread(tfnc,'y'));
	tf_data        = double(ncread(tfnc,'thermal_forcing'));
	z_data         = double(ncread(tfnc,'z'));

	%Build tf cell array
	t = 1:size(tf_data,4);
	tf = cell(1,1,size(tf_data,3));
	for i=1:size(tf_data,3)  %Iterate over depths
      temp_matrix=[];
		for ii=1:size(tf_data,4) %Iterate over time steps
			temp_tfdata=InterpFromGridToMesh(x_n,y_n,tf_data(:,:,i,ii)',md.mesh.x,md.mesh.y,0);
			temp_matrix = [temp_matrix temp_tfdata];
   	end
		temp_matrix = [temp_matrix ; t];
		tf(:,:,i)={temp_matrix};
	end
    obs_clim_tf = tf;
	save('Data/Ocean/ukesm_cold_clim_tf.mat','obs_clim_tf');

    %Set depths
    tf_depths = z_data';
    save('Data/Ocean/tf_depths.mat','tf_depths');
%
%	%Set ISMIP6 basal melt rate parameters
%	md.basalforcings            = basalforcingsismip6(md.basalforcings);
%	md.basalforcings.basin_id   = basinid;
%	md.basalforcings.num_basins = length(unique(basinid));
%	md.basalforcings.delta_t    = delta_t;
%	md.basalforcings.tf_depths  = tf_depths;
%	md.basalforcings.gamma_0    = 14477;
%	md.basalforcings.tf         = tf;
%
%	md.outputdefinition.definitions={};
%	md.timestepping.interp_forcings=1;
%	md.timestepping.final_time=106;
%	md.timestepping.time_step=1/12.;
%	%md.settings.output_frequency=24;
%	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','GroundedArea','FloatingArea','TotalSmb','SmbMassBalance','TotalGroundedBmb','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate',...
%		'IceVolumeScaled','IceVolumeAboveFloatationScaled','GroundedAreaScaled','FloatingAreaScaled','TotalSmbScaled','TotalGroundedBmbScaled','TotalFloatingBmbScaled','BasalforcingsIsmp6TfShelf'};
%
%	%Set melt / friction interpolation schemes
%	md.groundingline.migration = 'SubelementMigration';
%	md.groundingline.friction_interpolation='SubelementFriction1';
%	md.groundingline.melt_interpolation='SubelementMelt1';
%
%	%Solve
%	md.miscellaneous.name='ISMIP6_OCEAN_MELT_TEST';
%	md.verbose=verbose('solution',true,'module',true,'convergence',true);
%	md.cluster=cluster;
%	md.settings.waitonlock=0;
%	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
%
%	savemodel(org,md);
end %}}}
if perform(org,'1850-OcnForcingc')% {{{
	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');

	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;
	
	md.inversion.iscontrol=0;
	md.transient.isthermal=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.isgroundingline=1;
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);

	%% Load and remesh basin_id data
    basin_datanc=[data_ocean 'imbie2/imbie2_basin_numbers_8km.nc'];
	x_n            = double(ncread(basin_datanc,'x'));
	y_n            = double(ncread(basin_datanc,'y'));
	basinid_data   = double(ncread(basin_datanc,'basinNumber'));
	basinid        = InterpFromGridToMesh(x_n,y_n,basinid_data',md.mesh.x,md.mesh.y,0)+1;
	basinid        = mean(basinid(md.mesh.elements),2); %per element

	%Remove interpolation errors from basinid data
	for i=1:16
		pos = find(basinid>i-0.5 & basinid < i+0.5);
		basinid(pos) = i;
	end
	save('Data/Ocean/basinid.mat','basinid');

	tfnc           = [data_ocean 'ukesm1-0-ll_ssp585/1850-1994/UKESM1-0-LL_thermal_forcing_8km_x_60m.nc'];
	x_n            = double(ncread(tfnc,'x'));
	y_n            = double(ncread(tfnc,'y'));
	tf_data        = double(ncread(tfnc,'thermal_forcing'));
	z_data         = double(ncread(tfnc,'z'));

	%Build tf cell array
	t = 1:size(tf_data,4);
	tf = cell(1,1,size(tf_data,3));
	for i=1:size(tf_data,3)  %Iterate over depths
      temp_matrix=[];
		for ii=1:size(tf_data,4) %Iterate over time steps
			temp_tfdata=InterpFromGridToMesh(x_n,y_n,tf_data(:,:,i,ii)',md.mesh.x,md.mesh.y,0);
			temp_matrix = [temp_matrix temp_tfdata];
   	end
		temp_matrix = [temp_matrix ; t];
		tf(:,:,i)={temp_matrix};
	end
    obs_clim_tf = tf;
    save('Data/Ocean/ukesm_histo_tf.mat','obs_clim_tf','-v7.3');

    %Set depths
    tf_depths = z_data';
    save('Data/Ocean/tf_depths.mat','tf_depths');
%
%	%Set ISMIP6 basal melt rate parameters
%	md.basalforcings            = basalforcingsismip6(md.basalforcings);
%	md.basalforcings.basin_id   = basinid;
%	md.basalforcings.num_basins = length(unique(basinid));
%	md.basalforcings.delta_t    = delta_t;
%	md.basalforcings.tf_depths  = tf_depths;
%	md.basalforcings.gamma_0    = 14477;
%	md.basalforcings.tf         = tf;
%
%	md.outputdefinition.definitions={};
%	md.timestepping.interp_forcings=1;
%	md.timestepping.final_time=106;
%	md.timestepping.time_step=1/12.;
%	%md.settings.output_frequency=24;
%	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','GroundedArea','FloatingArea','TotalSmb','SmbMassBalance','TotalGroundedBmb','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate',...
%		'IceVolumeScaled','IceVolumeAboveFloatationScaled','GroundedAreaScaled','FloatingAreaScaled','TotalSmbScaled','TotalGroundedBmbScaled','TotalFloatingBmbScaled','BasalforcingsIsmp6TfShelf'};
%
%	%Set melt / friction interpolation schemes
%	md.groundingline.migration = 'SubelementMigration';
%	md.groundingline.friction_interpolation='SubelementFriction1';
%	md.groundingline.melt_interpolation='SubelementMelt1';
%
%	%Solve
%	md.miscellaneous.name='ISMIP6_OCEAN_MELT_TEST';
%	md.verbose=verbose('solution',true,'module',true,'convergence',true);
%	md.cluster=cluster;
%	md.settings.waitonlock=0;
%	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
%
%	savemodel(org,md);
end %}}}

%no mass changes determine basal meltrate 7-9
if perform(org,'control_1850_clim_nothkchange')% {{{
    % running wiht the non local  parameterization

	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');

    md.smb.mass_balance=interpRacmoSMB_ISMIP6(md,datadir);
	load 'Data/Atmosphere/ukesm_histo_clim_smb.mat';
    md.smb.mass_balance = miroc_rcp85_smb; %already in m/year ice

	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;
	
	md.inversion.iscontrol=0;
	md.transient.isthermal=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.isgroundingline=1;
    md.transient.issmb =0;
	md.masstransport.spcthickness=md.geometry.thickness;

	%Load forcing data
	load 'Data/Ocean/ukesm_histo_clim_tf.mat';
	load 'Data/Ocean/deltat_median.mat';
	load 'Data/Ocean/gamma0_median.mat';
	load 'Data/Ocean/basinid.mat';
	load 'Data/Ocean/tf_depths.mat';

	%Set ISMIP6 basal melt rate parameters
	delta_t                     = deltat_median;
	md.basalforcings            = basalforcingsismip6(md.basalforcings);
	md.basalforcings.basin_id   = basinid;
	md.basalforcings.num_basins = length(unique(basinid));
	md.basalforcings.delta_t    = delta_t;
	md.basalforcings.tf_depths  = tf_depths;
	md.basalforcings.gamma_0    = gamma0_median;
	md.basalforcings.tf         = obs_clim_tf;
    md.basalforcings.islocal = 0;

	%Model specifications
	md.outputdefinition.definitions={};
	md.timestepping.interp_forcing=0;

    md.timestepping.final_time=1/12;
    md.timestepping.time_step=1/12.;
    md.settings.output_frequency=1;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','GroundedArea','FloatingArea','TotalSmb','SmbMassBalance','TotalGroundedBmb','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate',...
		'IceVolumeScaled','IceVolumeAboveFloatationScaled','GroundedAreaScaled','FloatingAreaScaled','TotalSmbScaled','TotalGroundedBmbScaled','TotalFloatingBmbScaled'};

	%Set melt / friction interpolation schemes
	md.groundingline.migration = 'SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='SubelementMelt1';
%
	%Solve
	md.miscellaneous.name='1850_clim_nothkchange';
	md.verbose=verbose('solution',true,'module',true,'convergence',true);
    clustername = 'oshostname';
    cluster = set_cluster(clustername);
    md.cluster =cluster;
	md.settings.waitonlock=Inf;
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
    savemodel(org,md);

end %}}}
if perform(org,'control_1995_clim_nothkchange')% {{{
    % running wiht the non local  parameterization

	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');

    md.smb.mass_balance=interpRacmoSMB_ISMIP6(md,datadir);
	load 'Data/Atmosphere/ukesm_histo_clim_smb.mat';
    md.smb.mass_balance = miroc_rcp85_smb; %already in m/year ice

	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;
	
	md.inversion.iscontrol=0;
	md.transient.isthermal=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.isgroundingline=1;
    md.transient.issmb =0;
	md.masstransport.spcthickness=md.geometry.thickness;

	%Load forcing data
    load './../ISMIP6/Data/Ocean/obs_clim_tf.mat';
	load 'Data/Ocean/deltat_median.mat';
	load 'Data/Ocean/gamma0_median.mat';
	load 'Data/Ocean/basinid.mat';
	load 'Data/Ocean/tf_depths.mat';

	%Set ISMIP6 basal melt rate parameters
	delta_t                     = deltat_median;
	md.basalforcings            = basalforcingsismip6(md.basalforcings);
	md.basalforcings.basin_id   = basinid;
	md.basalforcings.num_basins = length(unique(basinid));
	md.basalforcings.delta_t    = delta_t;
	md.basalforcings.tf_depths  = tf_depths;
	md.basalforcings.gamma_0    = gamma0_median;
	md.basalforcings.tf         = obs_clim_tf;
    md.basalforcings.islocal = 0;

	%Model specifications
	md.outputdefinition.definitions={};
	md.timestepping.interp_forcing=0;

    md.timestepping.final_time=1/12;
    md.timestepping.time_step=1/12.;
    md.settings.output_frequency=1;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','GroundedArea','FloatingArea','TotalSmb','SmbMassBalance','TotalGroundedBmb','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate',...
		'IceVolumeScaled','IceVolumeAboveFloatationScaled','GroundedAreaScaled','FloatingAreaScaled','TotalSmbScaled','TotalGroundedBmbScaled','TotalFloatingBmbScaled'};

	%Set melt / friction interpolation schemes
	md.groundingline.migration = 'SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='SubelementMelt1';
%
	%Solve
	md.miscellaneous.name='1995_clim_nothkchange';
	md.verbose=verbose('solution',true,'module',true,'convergence',true);
    clustername = 'oshostname';
    cluster = set_cluster(clustername);
    md.cluster =cluster;
	md.settings.waitonlock=Inf;
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
    savemodel(org,md);

end %}}}
if perform(org,'control_1850-1995_clim_nothkchange')% {{{

    % running wiht the non local  parameterization

	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');

    md.smb.mass_balance=interpRacmoSMB_ISMIP6(md,datadir);
	load 'Data/Atmosphere/ukesm_histo_clim_smb.mat';
    md.smb.mass_balance = miroc_rcp85_smb; %already in m/year ice

	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;
	
	md.inversion.iscontrol=0;
	md.transient.isthermal=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=0;
	md.transient.isgroundingline=0;
    md.transient.issmb =0;
	md.masstransport.spcthickness=md.geometry.thickness;

	%Load forcing data
	load 'Data/Ocean/ukesm_histo_tf.mat';
	load 'Data/Ocean/deltat_median.mat';
	load 'Data/Ocean/gamma0_median.mat';
	load 'Data/Ocean/basinid.mat';
	load 'Data/Ocean/tf_depths.mat';

	%Set ISMIP6 basal melt rate parameters
	delta_t                     = deltat_median;
	md.basalforcings            = basalforcingsismip6(md.basalforcings);
	md.basalforcings.basin_id   = basinid;
	md.basalforcings.num_basins = length(unique(basinid));
	md.basalforcings.delta_t    = delta_t;
	md.basalforcings.tf_depths  = tf_depths;
	md.basalforcings.gamma_0    = gamma0_median;
	md.basalforcings.tf         = obs_clim_tf;
    md.basalforcings.islocal = 0;

	%Model specifications
	md.outputdefinition.definitions={};
	md.timestepping.interp_forcing=0;

    md.timestepping.final_time=145;
    % md.timestepping.final_time=1/12;
    md.timestepping.time_step=1/12.;
    md.settings.output_frequency=12;
	md.transient.requested_outputs={'TotalGroundedBmb','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate'};

	%Set melt / friction interpolation schemes
	md.groundingline.migration = 'SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='SubelementMelt1';
%
	%Solve
	md.miscellaneous.name='1850-1995_clim_nothkchange';
	md.verbose=verbose('solution',true,'module',true,'convergence',true);
    clustername = 'gadi';
    % clustername = 'oshostname';
    cluster = set_cluster(clustername);
    md.cluster =cluster;
	md.settings.waitonlock=0;
	% md.settings.waitonlock=Inf;
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
    savemodel(org,md);

end %}}}
% clim runs
if perform(org,'control_1850_clim_nonlocal_10ka')% {{{
    % running wiht the non local  parameterization

	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');

    md.smb.mass_balance=interpRacmoSMB_ISMIP6(md,datadir);
	load 'Data/Atmosphere/ukesm_histo_clim_smb.mat';
    md.smb.mass_balance = miroc_rcp85_smb; %already in m/year ice

	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;
	
	md.inversion.iscontrol=0;
	md.transient.isthermal=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.isgroundingline=1;
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);

	%Load forcing data
	load 'Data/Ocean/ukesm_histo_clim_tf.mat';
	load 'Data/Ocean/deltat_median.mat';
	load 'Data/Ocean/gamma0_median.mat';
	load 'Data/Ocean/basinid.mat';
	load 'Data/Ocean/tf_depths.mat';

	%Set ISMIP6 basal melt rate parameters
	delta_t                     = deltat_median;
	md.basalforcings            = basalforcingsismip6(md.basalforcings);
	md.basalforcings.basin_id   = basinid;
	md.basalforcings.num_basins = length(unique(basinid));
	md.basalforcings.delta_t    = delta_t;
	md.basalforcings.tf_depths  = tf_depths;
	md.basalforcings.gamma_0    = gamma0_median;
	md.basalforcings.tf         = obs_clim_tf;
    md.basalforcings.islocal = 0;

	%Model specifications
	md.outputdefinition.definitions={};
	md.timestepping.interp_forcing=0;

    md.timestepping.final_time=10000;
    md.timestepping.time_step=1/12.;
    md.settings.output_frequency=100;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','GroundedArea','FloatingArea','TotalSmb','SmbMassBalance','TotalGroundedBmb','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate',...
		'IceVolumeScaled','IceVolumeAboveFloatationScaled','GroundedAreaScaled','FloatingAreaScaled','TotalSmbScaled','TotalGroundedBmbScaled','TotalFloatingBmbScaled'};

	%Set melt / friction interpolation schemes
	md.groundingline.migration = 'SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='SubelementMelt1';
%
	%Solve
	md.miscellaneous.name='1850_clim_10ka';
	md.verbose=verbose('solution',true,'module',true,'convergence',true);
    clustername = 'gadi';
    cluster = set_cluster(clustername);
	md.settings.waitonlock=0;
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
    if loadonly,
        savemodel(org,md);
    end

end %}}}
if perform(org,'cold_clim_nonlocal_1ka_Cfriction_nn')% {{{
    % running wiht the non local  parameterization

	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');
    x=[min(md.mesh.x):1000:max(md.mesh.x)];   
    y=[min(md.mesh.y):1000:max(md.mesh.y)];   
    c_data = double(ncread("Data/Fields/friction_coefficient_nn.nc","friction_c"));
    friction_c = InterpFromGrid(y,x,c_data,md.mesh.y,md.mesh.x,'nearest');
    pos = find(md.mask.ocean_levelset <0);
    md.friction.coefficient(pos) = friction_c(pos);  


    md.smb.mass_balance=interpRacmoSMB_ISMIP6(md,datadir);
	load 'Data/Atmosphere/ukesm_cold_clim_smb.mat';
    md.smb.mass_balance = miroc_rcp85_smb; %already in m/year ice

	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;
	
	md.inversion.iscontrol=0;
	md.transient.isthermal=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.isgroundingline=1;
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);

	%Load forcing data
	load 'Data/Ocean/ukesm_cold_clim_tf.mat';
	load 'Data/Ocean/deltat_median.mat';
	load 'Data/Ocean/gamma0_median.mat';
	load 'Data/Ocean/basinid.mat';
	load 'Data/Ocean/tf_depths.mat';

	%Set ISMIP6 basal melt rate parameters
	delta_t                     = deltat_median;
	md.basalforcings            = basalforcingsismip6(md.basalforcings);
	md.basalforcings.basin_id   = basinid;
	md.basalforcings.num_basins = length(unique(basinid));
	md.basalforcings.delta_t    = delta_t;
	md.basalforcings.tf_depths  = tf_depths;
	md.basalforcings.gamma_0    = gamma0_median;
	md.basalforcings.tf         = obs_clim_tf;
    md.basalforcings.islocal = 0;

	%Model specifications
	md.outputdefinition.definitions={};
	md.timestepping.interp_forcing=0;

    md.timestepping.final_time=1000;
    md.timestepping.time_step=1/12.;
    md.settings.output_frequency=12*50;
	md.transient.requested_outputs={'default','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate'};

	%Set melt / friction interpolation schemes
	md.groundingline.migration = 'SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='SubelementMelt1';
%
	%Solve
	md.miscellaneous.name='cold_clim_1ka_cnn';
	md.verbose=verbose('solution',true,'module',true,'convergence',true);
    clustername = 'gadi';
    cluster = set_cluster(clustername);
	md.settings.waitonlock=0;
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
    if loadonly,
        savemodel(org,md);
    end

end %}}}
if perform(org,'nobasal_melt_nonlocal_1ka_Cfriction_nn')% {{{
    % running wiht the non local  parameterization

	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');
    x=[min(md.mesh.x):1000:max(md.mesh.x)];   
    y=[min(md.mesh.y):1000:max(md.mesh.y)];   
    c_data = double(ncread("Data/Fields/friction_coefficient_nn.nc","friction_c"));
    friction_c = InterpFromGrid(y,x,c_data,md.mesh.y,md.mesh.x,'nearest');
    pos = find(md.mask.ocean_levelset <0);
    md.friction.coefficient(pos) = friction_c(pos);  


    md.smb.mass_balance=interpRacmoSMB_ISMIP6(md,datadir);
	load 'Data/Atmosphere/ukesm_cold_clim_smb.mat';
    md.smb.mass_balance = miroc_rcp85_smb; %already in m/year ice

	m=((1+sin(71*pi/180))*ones(md.mesh.numberofvertices,1)./(1+sin(abs(md.mesh.lat)*pi/180)));
	md.mesh.scale_factor=(1./m).^2;
	
	md.inversion.iscontrol=0;
	md.transient.isthermal=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.isgroundingline=1;
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);


    pos = md.basalforcings.floatingice_melting_rate >0;
    md.basalforcings.floatingice_melting_rate(pos)=0;
	%Model specifications
	md.outputdefinition.definitions={};
	md.timestepping.interp_forcing=0;

    md.timestepping.final_time=1000;
    md.timestepping.time_step=1/12.;
    md.settings.output_frequency=12*50;
	md.transient.requested_outputs={'default','TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate'};

	%Set melt / friction interpolation schemes
	md.groundingline.migration = 'SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='SubelementMelt1';
%
	%Solve
	md.miscellaneous.name='nobasalmelt_clim_1ka_cnn';
	md.verbose=verbose('solution',true,'module',true,'convergence',true);
    clustername = 'gadi';
    cluster = set_cluster(clustername);
	md.settings.waitonlock=0;
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
    if loadonly,
        savemodel(org,md);
    end

end %}}}
%output input fields
if perform(org,'2d_output_fields')% {{{


	% md=loadmodel('./Models/ISMIP6Antarctica__CollapseSSA.math');
	% md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');
    md =loadmodel('./../ISMIP6/Models/ISMIP6Antarctica_InversionFriction1.mat')
    x=[min(md.mesh.x):1000:max(md.mesh.x)];   
    y=[min(md.mesh.y):1000:max(md.mesh.y)];   
    % datapath='data_smb;
    friction_c = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.friction.coefficient,x,y,NaN);
    nccreate("Data/Fields/friction_coefficient_SSA.nc","friction_c",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    ncwrite("Data/Fields/friction_coefficient_SSA.nc","friction_c",transpose(friction_c));
    
    %expcontourlevelzero(md,md.mask.ice_levelset,-0.0001,'Exp/Level0.exp');mask_gr = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ocean_levelset,x,y,NaN);

    % nccreate("Data/Fields/mask_grounded.nc","mask_grounded",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/mask_grounded.nc","mask_grounded",transpose(mask_gr));

    % bed = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.bed,x,y,NaN);
    % nccreate("Data/Fields/bed.nc","bed",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/bed.nc","bed",transpose(bed));

    % vy = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.initialization.vy,x,y,NaN);
    % nccreate("Data/Fields/vy.nc","vy",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/vy.nc","vy",transpose(vy));
    % vy = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.initialization.vx,x,y,NaN);
    % nccreate("Data/Fields/vx.nc","vx",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/vx.nc","vx",transpose(vy));
    % vx = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.initialization.vx,x,y,NaN);
    % nccreate("Data/Fields/vx.nc","vx",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % % ncwrite("Data/Fields/vx.nc","vx",transpose(vx));
    % vel= InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.initialization.vel,x,y,NaN);
    % nccreate("Data/Fields/vel.nc","vel",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/vel.nc","vel",transpose(vel));
    % velo= InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y, md.inversion.vel_obs,x,y,NaN);
    % nccreate("Data/Fields/velo.nc","velo",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/velo.nc","velo",transpose(velo));

    % thickness= InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.thickness,x,y,NaN);
    % nccreate("Data/Fields/thickness.nc","thickness",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/thickness.nc","thickness",transpose(thickness));
    

    % surface= InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.surface,x,y,NaN);
    % nccreate("Data/Fields/surface.nc","surface",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/surface.nc","surface",transpose(surface));


    % [sx,sy,s]=slope(md); sslope=averaging(md,s,10);
    % disp('      -- Process surface velocity data');
    % vel = md.inversion.vel_obs;
    % flags=(vel==0); pos1=find(flags); pos2=find(~flags);
    % vel(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),vel(pos2),md.mesh.x(pos1),md.mesh.y(pos1));
    % vel=max(vel,0.1);
    % disp('      -- Calculate effective pressure');
    % Neff = md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base;
    % Neff(find(Neff<=0))=1;
    % disp('      -- Deduce friction coefficient');
    % friction_c_calculated=sqrt(md.materials.rho_ice*md.geometry.thickness.*(sslope)./(Neff.*vel/md.constants.yts));
    % friction_c_calculated=min(friction_c_calculated,2000);

    
    
    % friction_c_calculated= InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,friction_c_calculated,x,y,NaN);
    % nccreate("Data/Fields/friction_c_calculated.nc","friction_c_calculated",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/friction_c_calculated.nc","friction_c_calculated",transpose(friction_c_calculated));

    % msk= InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ocean_levelset,x,y,NaN);
    % nccreate("Data/Fields/msk.nc","msk",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/msk.nc","msk",transpose(msk));
    
    % contours=isoline(md, md.mask.ocean_levelset,'value',0,'output','Data/Fields.gl.exp');
    % load 'Data/Ocean/basinid.mat';

    % basin= InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,basinid,x,y,NaN);
    % nccreate("Data/Fields/basin.nc","basin",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/basin.nc","basin",transpose(basin));

    % geothermalflux= InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.basalforcings.geothermalflux,x,y,NaN);
    % nccreate("Data/Fields/geothermalflux.nc","geothermalflux",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/geothermalflux.nc","geothermalflux",transpose(geothermalflux));
    
    % interp field from bedtype
    % p_vel = '/Users/jbec0008/SAEF/datasets/LL-Geo-AntarcticBasins-3a4ea69/';
    % vel_set = "bedtype.nc";
    % bedm = fullfile(p_vel, vel_set);
    % bedtype_mesh = interpBedtypeAntartica(md.mesh.x,md.mesh.y,'bedtype','nearest',bedm);
    % bedtype= InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,bedtype_mesh,x,y,NaN);
    % nccreate("Data/Fields/bedtype.nc","bedtype",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/bedtype.nc","bedtype",transpose(bedtype));
end %}}}
if perform(org,'2d_input_fields')% {{{


	% md=loadmodel('./Models/ISMIP6Antarctica__CollapseSSA.math');
	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');
    c_data = double(ncread("Data/Fields/friction_coefficient_nn.nc","friction_c"));
    x=[min(md.mesh.x):1000:max(md.mesh.x)];   
    y=[min(md.mesh.y):1000:max(md.mesh.y)];   
    friction_c = InterpFromGrid(y,x,c_data,md.mesh.y,md.mesh.x,'nearest');
    m_data = double(ncread("Data/Fields/retreated_diff_mask_grounded.nc","mask_grounded"));
    mask_retreat = InterpFromGrid(y,x,m_data,md.mesh.y,md.mesh.x,'nearest');
    % datapath='data_smb;
    % friction_c = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.friction.coefficient,x,y,NaN);
    % nccreate("Data/Fields/friction_coefficient.nc","friction_c",'Dimensions',{'x' length(x) 'y' length(y)},'Format','classic');
    % ncwrite("Data/Fields/friction_coefficient.nc","friction_c",transpose(friction_c));
    
end %}}}

if perform(org,'InversionFriction1_fromSSAcollapse'),% {{{

	% md=loadmodel(org,'InversionB');
	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');

	%Remove some parts of the domain
	pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/IceAroundShelves.exp',2));
	md.mask.ice_levelset(pos)=-1;

	%Control general
	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.iscontrol=1;
	md.verbose=verbose('solution',false,'control',true);

	%Find all elements that include ice and constrain some of them
	% pos_e = find(min(md.mask.ice_levelset(md.mesh.elements),[],2)<0);
	% flags=zeros(md.mesh.numberofvertices,1); flags(md.mesh.elements(pos_e,:))=1;
	% pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/Constrain.exp',2) & flags);
	% md.stressbalance.spcvx(pos)=md.inversion.vx_obs(pos);
	% md.stressbalance.spcvy(pos)=md.inversion.vy_obs(pos);
	% md.stressbalance.spcvz(pos)=0;

	%Cost functions
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
	md.inversion.cost_functions_coefficients(:,1)=2000;
	md.inversion.cost_functions_coefficients(:,2)=5;
	md.inversion.cost_functions_coefficients(:,3)=.2*50^-3;
	pos=find(md.mask.ice_levelset>0 | md.inversion.vel_obs==0);% no ice or with velocity 0
	md.inversion.cost_functions_coefficients(pos,1:2)=0;

	% PIG correction
    posPig=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/PIG.exp',2));
    md.inversion.cost_functions_coefficients(posPig,1)=0.5;
    md.inversion.cost_functions_coefficients(posPig,2)=2000;
	% md.inversion.cost_functions_coefficients(posPig,3)=.2*50^-2;
    posPig2=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/THW.exp',2));
    md.inversion.cost_functions_coefficients(posPig2,1)=0.5;
    md.inversion.cost_functions_coefficients(posPig2,2)=2000;
	% md.inversion.cost_functions_coefficients(posPig2,3)=.2*50^-2;

	%Controls
	md.groundingline.migration='SubelementMigration'; %VERY important at this resolution
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='SubelementMelt1';
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.maxsteps=300;
	md.inversion.maxiter =300;
	md.inversion.min_parameters=0.05*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=2000*ones(md.mesh.numberofvertices,1);
	md.inversion.control_scaling_factors=1;

    md.miscellaneous.name = 'SSAcollapser_friction_nochange'
    md.settings.waitonlock=0;
	%Go solve/read_clusterjobfile
    % md.miscellaneous.name='SSA_fritction';
    clustername = 'gadi';
    cluster = set_cluster(clustername);
	md.cluster=cluster;
    md=solve(md,'sb','runtimename', false,'loadonly',loadonly);
    if loadonly,
        % rs =read_clusterjobfile('cluster_job_friction');
        % disp(rs);
        % mds=loadresultsfromcluster(md,'runtimename',rs);
        % md.results = mds.results;
        md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;


        savemodel(org,md);
        plot_after_tuning(1,md);
    end
end%}}}
if perform(org,'InversionFriction1_fromSSAcollapse_retreat'),% {{{

	% md=loadmodel(org,'InversionB');
	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');

    c_data = double(ncread("Data/Fields/friction_coefficient_nn.nc","friction_c"));
    x=[min(md.mesh.x):1000:max(md.mesh.x)];   
    y=[min(md.mesh.y):1000:max(md.mesh.y)];   
    friction_c = InterpFromGrid(y,x,c_data,md.mesh.y,md.mesh.x,'nearest');
    m_data = double(ncread("Data/Fields/retreated_diff_mask_grounded.nc","mask_grounded"));
    mask_retreat = InterpFromGrid(y,x,m_data,md.mesh.y,md.mesh.x,'nearest');
    md.friction.coefficient=friction_c;
	%Remove some parts of the domain
	pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/IceAroundShelves.exp',2));
	md.mask.ice_levelset(pos)=-1;

	%Control general
	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.iscontrol=1;
	md.verbose=verbose('solution',false,'control',true);

	%Find all elements that include ice and constrain some of them
	% pos_e = find(min(md.mask.ice_levelset(md.mesh.elements),[],2)<0);
	% flags=zeros(md.mesh.numberofvertices,1); flags(md.mesh.elements(pos_e,:))=1;
	% pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/Constrain.exp',2) & flags);
	% md.stressbalance.spcvx(pos)=md.inversion.vx_obs(pos);
	% md.stressbalance.spcvy(pos)=md.inversion.vy_obs(pos);
	% md.stressbalance.spcvz(pos)=0;

	%Cost functions
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
	md.inversion.cost_functions_coefficients(:,1)=2000;
	md.inversion.cost_functions_coefficients(:,2)=5;
	md.inversion.cost_functions_coefficients(:,3)=.2*50^-3;
	pos=find(md.mask.ice_levelset>0 | md.inversion.vel_obs==0);% no ice or with velocity 0
	md.inversion.cost_functions_coefficients(pos,1:2)=0;
    pos = find( mask_retreat==1);
	md.inversion.cost_functions_coefficients(pos,1:2)=0;

	% PIG correction
    posPig=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/PIG.exp',2));
    md.inversion.cost_functions_coefficients(posPig,1)=0.5;
    md.inversion.cost_functions_coefficients(posPig,2)=2000;
	% md.inversion.cost_functions_coefficients(posPig,3)=.2*50^-2;
    posPig2=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/THW.exp',2));
    md.inversion.cost_functions_coefficients(posPig2,1)=0.5;
    md.inversion.cost_functions_coefficients(posPig2,2)=2000;
	% md.inversion.cost_functions_coefficients(posPig2,3)=.2*50^-2;

	%Controls
	md.groundingline.migration='SubelementMigration'; %VERY important at this resolution
	md.groundingline.friction_interpolation='SubelementFriction1';
	md.groundingline.melt_interpolation='SubelementMelt1';
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.maxsteps=300;
	md.inversion.maxiter =300;
	md.inversion.min_parameters=0.05*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=2000*ones(md.mesh.numberofvertices,1);
	md.inversion.control_scaling_factors=1;

    md.miscellaneous.name = 'SSAcollapser_friction_retreated'
    md.settings.waitonlock=0;
	%Go solve/read_clusterjobfile
    % md.miscellaneous.name='SSA_fritction';
    clustername = 'gadi';
    cluster = set_cluster(clustername);
	md.cluster=cluster;
    md=solve(md,'sb','runtimename', false,'loadonly',loadonly);
    if loadonly,
        % rs =read_clusterjobfile('cluster_job_friction');
        % disp(rs);
        % mds=loadresultsfromcluster(md,'runtimename',rs);
        % md.results = mds.results;
        pos = find( mask_retreat==1);
        md.friction.coefficient(pos)=md.results.StressbalanceSolution.FrictionCoefficient(pos);


        savemodel(org,md);
        plot_after_tuning(13,md);
    end
end%}}}
if perform(org,'Inversion2HO_cretreat'),% {{{
	md=loadmodel('./../ISMIP6/Models/Cycle3_Inversion2B.mat');
    % use steadty state tmperature
    c_data = double(ncread("Data/Fields/friction_coefficient_nn.nc","friction_c"));
    x=[min(md.mesh.x):1000:max(md.mesh.x)];   
    y=[min(md.mesh.y):1000:max(md.mesh.y)];   
    friction_c = InterpFromGrid(y,x,c_data,md.mesh.y,md.mesh.x,'nearest');
    m_data = double(ncread("Data/Fields/retreated_diff_mask_grounded.nc","mask_grounded"));
    mask_retreat = InterpFromGrid(y,x,m_data,md.mesh.y,md.mesh.x,'nearest');
    md.friction.coefficient=friction_c;

	%Control general
	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.iscontrol=1;
	md.verbose=verbose('solution',false,'control',true);

	%Cost functions
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
	md.inversion.cost_functions_coefficients(:,1)=2000*1.5;
	md.inversion.cost_functions_coefficients(:,2)=5;
	md.inversion.cost_functions_coefficients(:,3)=2.5*1e-4;% from L curve analysis
	% pos=find(md.mask.ice_levelset>0 | md.inversion.vel_obs==0);% no ice or with velocity 0
	pos=find(md.mask.ice_levelset>0 | md.inversion.vel_obs<0.1);% no ice or with velocity 0
    % WRONG mask level set ,0 bvel 0.1
    % check for equal weigthing
    %save model if loadonly=1
  
	md.inversion.cost_functions_coefficients(pos,1:2)=0;
    pos = find( mask_retreat==1);
	md.inversion.cost_functions_coefficients(pos,1:2)=0;


	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.iscontrol=1;
	md.verbose=verbose('solution',false,'control',true);


	md.miscellaneous.name='Antarctica_HOfriciton_cretreat';

	%Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.maxsteps=200;
	md.inversion.maxiter =200;
	% md.inversion.maxsteps=20;
	% md.inversion.maxiter =20;
	md.inversion.min_parameters=0.05*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=2000*ones(md.mesh.numberofvertices,1);
	md.inversion.control_scaling_factors=1;

	%Additional parameters
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;

	%Go solve
	md.inversion.iscontrol=1;
    clustername = 'gadi';
    cluster = set_cluster(clustername);
	md.cluster=cluster;
    md.settings.waitonlock=0;
    % md.settings.waitonlock=Inf;
    md=solve(md,'sb','runtimename',false,'loadonly',loadonly);
    if loadonly,
        savemodel(org,md);
    end
end% }}}
%flowline extrapolation
if perform(org,'flowlineextrapolation')% {{{
    modelDir='./Data/';


	% md=loadmodel('./Models/ISMIP6Antarctica__CollapseSSA.math');
	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');
    % Make tmp dircetory for exp files specifically for this function
    tmpDir = fullfile(modelDir, 'tmp');
    mkdir(tmpDir)

    % Extract the XY Coordinates along the GL
    GL = isoline(md, md.mask.ocean_levelset, 'value', 0);
    expwrite(GL, fullfile(tmpDir, 'GL.exp'));

    out = []; % Empty container for results
    % For each segment of the GL; extract each X/Y point; generate a
    % downstream flow line; append the flowline to the output object
    for i=1:length(GL)
        x0 = GL(:,i).x;
        y0 = GL(:,i).y;
        f = flowlines(md.mesh.elements, md.mesh.x, md.mesh.y, md.inversion.vx_obs, md.inversion.vy_obs, x0, y0, 'upstream', 0);
        out = [out; f];
    end

    expwrite(out, fullfile(tmpDir, 'GL_flowlines.exp')); % Write all flowlines to file
    plotmodel(md, 'data', md.mask.ocean_levelset, 'expdisp', fullfile(tmpDir, 'GL_flowlines.exp'), 'edgecolor', 'w');
    % x=[min(md.mesh.x):1000:max(md.mesh.x)];   
    % y=[min(md.mesh.y):1000:max(md.mesh.y)];   
end %}}}
if perform(org,'flowlineextrapolation_downstream')% {{{
    modelDir='./Data/';


	% md=loadmodel('./Models/ISMIP6Antarctica__CollapseSSA.math');
	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');
    % Make tmp dircetory for exp files specifically for this function
    tmpDir = fullfile(modelDir, 'tmp');
    mkdir(tmpDir)

    % Extract the XY Coordinates along the GL
    GL = isoline(md, md.mask.ocean_levelset, 'value', 0);
    expwrite(GL, fullfile(tmpDir, 'GL.exp'));

    out = []; % Empty container for results
    % For each segment of the GL; extract each X/Y point; generate a
    % downstream flow line; append the flowline to the output object
    for i=1:length(GL)
        x0 = GL(:,i).x;
        y0 = GL(:,i).y;
        f = flowlines(md.mesh.elements, md.mesh.x, md.mesh.y, md.inversion.vx_obs, md.inversion.vy_obs, x0, y0, 'downstream', 0);
        out = [out; f];
    end

    expwrite(out, fullfile(tmpDir, 'GL_flowlines_downstream.exp')); % Write all flowlines to file
    plotmodel(md, 'data', md.mask.ocean_levelset, 'expdisp', fullfile(tmpDir, 'GL_flowlines_downstream.exp'), 'edgecolor', 'w');
    % x=[min(md.mesh.x):1000:max(md.mesh.x)];   
    % y=[min(md.mesh.y):1000:max(md.mesh.y)];   
end %}}}
if perform(org,'flowlineextrapolation_velsmae')% {{{
    modelDir='./Data/';


	% md=loadmodel('./Models/ISMIP6Antarctica__CollapseSSA.math');
	md=loadmodel('Models/ISMIP6Antarctica_CollapseSSA.mat');
    % Make tmp dircetory for exp files specifically for this function
    tmpDir = fullfile(modelDir, 'tmp');
    mkdir(tmpDir)

    % Extract the XY Coordinates along the GL
    GL = isoline(md, md.mask.ocean_levelset, 'value', 0);
    expwrite(GL, fullfile(tmpDir, 'GL.exp'));

    out = []; % Empty container for results
    % For each segment of the GL; extract each X/Y point; generate a
    % downstream flow line; append the flowline to the output object
    for i=1:length(GL)
        x0 = GL(:,i).x;
        y0 = GL(:,i).y;
        f = flowlines(md.mesh.elements, md.mesh.x, md.mesh.y, md.inversion.vx_obs, md.inversion.vy_obs, x0, y0, 'upstream', 0);
        out = [out; f];
    end

    expwrite(out, fullfile(tmpDir, 'GL_flowlines.exp')); % Write all flowlines to file
    plotmodel(md, 'data', md.mask.ocean_levelset, 'expdisp', fullfile(tmpDir, 'GL_flowlines.exp'), 'edgecolor', 'w');
    % x=[min(md.mesh.x):1000:max(md.mesh.x)];   
    % y=[min(md.mesh.y):1000:max(md.mesh.y)];   
end %}}}
