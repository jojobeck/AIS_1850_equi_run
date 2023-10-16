function tskin_obs=interpRacmoTemp(md, datadir, varargin)
	disp('   Loading Racmo data from Netcdf'); 
	options=pairoptions(varargin{:});
	outdir=getfieldvalue(options,'outdir','');
	makelocalfile=getfieldvalue(options,'makelocalfile',0);

	% ncdata=[datadir 'racmo/RACMO2.3p2_ERA5_ANT27_Tskin_climatology_1979_01_2021_12.nc'];
	ncdata=[datadir 'racmo/Tskin_1995_2014_clim_annual_mean.nc'];
	%ncdata=[datadir 'racmo/tskin_RACMO2.3_climate_1979_2014_mvw.nc'];
	tskin=ncread(ncdata,'tskin');
	lat=ncread(ncdata,'lat');
	lon=ncread(ncdata,'lon');
	[xtskin ytskin]=ll2xy(lat, lon, -1);

	disp('   Set observed tskin')
	[ny,nx]=size(xtskin);
	x=reshape(xtskin, [ny*nx,1]);
	y=reshape(ytskin, [ny*nx,1]);
	tskinr=reshape(tskin, [ny*nx,1]);
	[xt,yt,tskin_obs]=griddata(x,y,tskinr,md.mesh.x,md.mesh.y);

	if makelocalfile
		cellsize=1000;
		xmin=floor(min(md.mesh.x)/cellsize)*cellsize;
		xmax=ceil(max(md.mesh.x)/cellsize)*cellsize;
		ymin=floor(min(md.mesh.y)/cellsize)*cellsize;
		ymax=ceil(max(md.mesh.y)/cellsize)*cellsize;

		xm=[xmin:cellsize:xmax]';
		ym=[ymin:cellsize:ymax]';
		[XM,YM]=meshgrid(xm,ym);
		[d,d,tdata]=griddata(x,y,tskinr,XM,YM);
		filename=[outdir 'templocal.mat'];
		save(filename,'xm','ym','tdata');
	end
