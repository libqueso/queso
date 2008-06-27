clear all;
path(path,getenv('MCMC_TOOLBOX_PATH_FOR_PECOS_TOOLKIT'));
path(path,pwd);
uqHimmelExOutput;

figure(1); clf
mcmcplot(chainCpp,[],resultsCpp.names,'chainpanel')
subplot(2,2,4)
mcmcplot(sqrt(s2chainCpp),[],[],'dens',2)
title('error std')

chainstats(chainCpp,resultsCpp);

data.ydata = [
            0     0.02090
         4.50     0.01540
         8.67     0.01422
        12.67     0.01335
        17.75     0.01232
        22.67     0.01181
        27.08     0.01139
        32.00     0.01092
        36.00     0.01054
        46.33     0.00978
        57.00     0.009157
        69.00     0.008594
        76.75     0.008395
        90.00     0.007891
       102.00     0.007510
       108.00     0.007370
       147.92     0.006646
       198.00     0.005883
       241.75     0.005322
       270.25     0.004960
       326.25     0.004518
       418.00     0.004075
       501.00     0.003715
    ];

data.y0 = [0.02090;0.02090/3;0;0;0];

figure(2); clf
[t,y] = ode45(@himmelode,linspace(0,600),data.y0,[],mean(chainCpp));
plot(data.ydata(:,1),data.ydata(:,2),'s',t,y,'-')
ylim([0,0.021])
legend('Aobs','A','B','C','D','E',0)
title('Data and fitted model')
