function cluster =set_cluster(clustername)
    if strcmpi(clustername,'gadi'),
        disp('--cluster gadi---------------');
        % cluster=gadi('np',30);
        cluster=gadi('np',48);
        % cluster=gadi('np',6);
        cluster.mem=190;
        % cluster.time=60*3;
        cluster.time=60*40;
        cluster.project='bi77';
        % cluster.project='au88';
        cluster.queue='normal';
    else
        cluster=generic('name',oshostname(),'np',6);
    end
