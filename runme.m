steps = [19];
steps = [30];

loadonly = 0;
loadonly = 1;
addpath('scripts');



org=organizer('repository',['./Models'],'prefix',['ISMIP6Antarctica_'],'steps',steps, 'color', '34;47;2'); 
clear steps;
datadir= '/Users/jbec0008/SAEF/datasets/';
% loadonly=1;
