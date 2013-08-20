path(path,'/Users/prudenci/works/personal/svns/sa/stalg/gp_lanl/code_orig/gpmsa_Sourceforge_0609_changed_locally/matlab');
load('gpmsa_pout.mat');

global debugPrudenci
global counterPrudenci
debugPrudenci = 0
counterPrudenci = 0

%%%%%%%%%%%%%%%%%%%%

from = 2000; % ernesto 200 2000
to = length(pout.pvals);
thismany = 500; % ernesto 100 500
ip = round(linspace(from,to,thismany));

PCresponsesurf(pout,ip);

%%%%%%%%%%%%%%%%%%%%

