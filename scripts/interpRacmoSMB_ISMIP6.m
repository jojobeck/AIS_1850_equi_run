function smb_obs=interpRacmoSMB_ISMIP6(md, datadir, varargin)                                                                                                                                                                                                                          
    disp('   Loading Racmo data from Netcdf');
    options=pairoptions(varargin{:});
    outdir=getfieldvalue(options,'outdir','');
    makelocalfile=getfieldvalue(options,'makelocalfile',0);
    method=getfieldvalue(options,'method','natural');
 
    ncdata=[datadir 'racmo/RACMO2.3p2_ERA5_ANT27_smb_monthly_1979_01_2021_12.nc'];
    %ncdata=[datadir 'racmo/SMB_RACMO2.3_monthly_ANT27_1979_2013.nc'];
    smb=ncread(ncdata,'smb');
    smb=squeeze(smb);
    [nx ny nt]=size(smb);
    smbs=sum(reshape(smb,[nx,ny,nt/12,12]),4); %yearly smb
    smbm=mean(smbs(:,:,17:37),3); %mean m/year over period 1995 -2015
    smb=smbm/md.materials.rho_ice; %m/year ice eq.
    lat=ncread(ncdata,'lat');
    lon=ncread(ncdata,'lon');
    [xsmb ysmb]=ll2xy(lat, lon, -1);
 
    disp('   Set observed smb')
    %x=reshape(xsmb, [ny*nx,1]);
    %y=reshape(ysmb, [ny*nx,1]);
    %smbsr=reshape(smbs, [ny*nx,1]);
    [xr,yr,smb_obs]=griddata(xsmb,ysmb,smb,md.mesh.x,md.mesh.y,method);
    %smb_obs=smb_obs*md.materials.rho_water/md.materials.rho_ice;
 
    if makelocalfile
        cellsize=1000;
        xmin=floor(min(md.mesh.x)/cellsize)*cellsize;
        xmax=ceil(max(md.mesh.x)/cellsize)*cellsize;
        ymin=floor(min(md.mesh.y)/cellsize)*cellsize;
        ymax=ceil(max(md.mesh.y)/cellsize)*cellsize;
 
        xm=[xmin:cellsize:xmax]';
        ym=[ymin:cellsize:ymax]';
        [XM,YM]=meshgrid(xm,ym);
        [d,d,smbdata]=griddata(x,y,smbsr,XM,YM);
        filename=[outdir 'smblocal.mat'];
        save(filename,'xm','ym','smbdata');
    end
