function [S,X,glme_posterior] = getCleanSpikeRates(S,varargin)
%GETCLEANSPIKERATES Return smoothed single-trial fits using Poisson GLME
%
%  [S,X,glme_posterior,glme_prior] = stat.getCleanSpikeRates(S,'Name',value,...);
%
% Inputs
%  S - Table returned using init.spikeData();
%
% Output
%  S - Same table, but with additional variable "lambda" that is parsed by 
%        1) Find mean rate during successful outcomes
%        2) Assign each column a corresponding timestamp: now we will try
%        and fit observations, which have a combination of
%           a) Channel
%              -> Channels are nested in Probe
%           b) Timestep
%           c) Observed spike count (posterior; random effect)
%           d) Grouping of trial outcome
%           e) Duration between reach and grasp
%        3) We will recover lambda, the Poisson rate parameter that is the
%           most likely prior for underlying data given the observations.
%  X - "Expanded" table where each row represents a timestep instead of a
%        single trial.
%  glme_posterior - The model that is used to get the "smoothed output"
%                   (N_fit, CI_fit)
%  glme_prior     - Model used to detrend on a per-channel basis.

pars = struct;
pars.Knots = 1:6:60;
pars.MinSpikeRate = 5; % Spikes/sec
pars.SGOrder = 3;
pars.SGFrameLen = 21;
fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx) == 1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

t = S.Properties.UserData.t;
spikeThresh = pars.MinSpikeRate * (t(end)-t(1)); % Threshold is in terms of total spikes

S.Duration = S.Grasp - S.Reach;
S = movevars(S,'Duration','before','Outcome');

G = findgroups(S(:,{'TrialID','Probe','Channel','Outcome','Duration'}));

% Rearrange data table so that every element has a corresponding time
X = splitapply(@(TrialID,Outcome,Duration,Probe,Channel,Spikes)expandS2X(TrialID,Outcome,Duration,Probe,Channel,Spikes,t),...
      S.TrialID,S.Outcome,S.Duration,...
      S.Probe,S.Channel,S.Spikes,G);
X = vertcat(X{:});
   
% Now, fit using generalized linear mixed-effects model. First step is to
% create a corresponding "detrending" term (average will work well for
% that).
% [Gx,XID] = findgroups(X(:,{'Probe','Channel','t'}));
% XID.Spikes = splitapply(@(Spikes)round(nanmean(Spikes)),X.Spikes,Gx);
% XID.N = splitapply(@(N)round(nanmean(N)),X.N,Gx);

% BinomialSize = XID.N;
% XID.Exclude = BinomialSize < spikeThresh;
% BinomialSize(XID.Exclude) = 1; % Just to make sure fitglme works
% XID.ChannelID = categorical(strcat(XID.Probe,"::",XID.Channel));

% % Get generic trends with time
% tic; fprintf(1,'Fitting GLME for <strong>channel</strong> time trends...');
% glme_prior = fitglme(XID,'Spikes~1+(1+t_ord|ChannelID)',...
%    'BinomialSize',BinomialSize,...
%    'FitMethod','REMPL',...
%    'Distribution','Binomial',...
%    'Link','logit',...
%    'Exclude',XID.Exclude);
% fprintf(1,'complete (%5.2f sec)\n',toc);
% 
% disp(glme_prior);
% disp('<strong>Rsquared:</strong>');
% disp(glme_prior.Rsquared);

% [p,ci] = predict(glme_prior,XID);
% XID.N_prior = p.*XID.N;
% XID.CI_prior = ci.*XID.N;
% XID.N_Spline = splitapply(@(Spikes,t)fitspline(Spikes,t,pars),X.Spikes,X.t,Gx);
% X = outerjoin(X,XID,...
%    'Keys',{'Probe','Channel','t'},...
%    'MergeKeys',true,...
%    'LeftVariables',{'TrialID','Probe','Channel','Outcome','Duration','Spikes','N'},...
%    'RightVariables',{'ChannelID','t','t_2','t_3','N_Spline','Exclude'});
% [p,ci] = predict(glme_prior,X);
% X.N_post = p.*X.N;
% X.CI_post = ci.*X.N;
% X.TrialIDc = categorical(X.TrialID);

X.ChannelID = categorical(strcat(X.Probe,"::",X.Channel));
[Gx,XID] = findgroups(X(:,'ChannelID'));
XID.N_Total_Expect = splitapply(@(N)nanmean(N),X.N,Gx);
XID.Exclude = XID.N_Total_Expect < spikeThresh;
tmp = splitapply(@(Spikes,t,ChannelID)fitspline(Spikes,t,ChannelID,pars),X.Spikes,X.t,X.ChannelID,Gx);
tmp = vertcat(tmp{:});

X = outerjoin(X,XID,...
   'Keys',{'ChannelID'},...
   'MergeKeys',true,...
   'LeftVariables',X.Properties.VariableNames,...
   'RightVariables',{'Exclude','N_Total_Expect'});

X = outerjoin(X,tmp,...
   'Keys',{'ChannelID','t'},...
   'MergeKeys',true,...
   'LeftVariables',X.Properties.VariableNames,...
   'RightVariables',{'N_Spline'});

X.Properties.UserData = S.Properties.UserData;
X.Properties.UserData.NSpikeMinimum = spikeThresh;

% % Now, fit the posterior model
% BinomialSize = X.N;
% BinomialSize(X.N < 1) = 1;

tic; fprintf(1,'Fitting GLME for <strong>trial-specific</strong> trends...');
% glme_posterior = fitglme(X,'Spikes~1+Outcome+(1+N_Spline+Duration*t|TrialID)',...
%    'BinomialSize',BinomialSize,...
%    'FitMethod','REMPL',...
%    'Distribution','Binomial',...
%    'Link','logit',...
%    'Exclude',(X.Exclude) | (X.N < 1));
glme_posterior = fitglme(X,'Spikes~1+Outcome+(1+N_Spline|TrialID)+(1+N+Duration*t|ChannelID)',...
   'FitMethod','REMPL',...
   'Distribution','Poisson',...
   'Link','log',...
   'Exclude',(X.Exclude) | (X.N < 1));
fprintf(1,'complete (%5.2f sec)\n',toc);

disp(glme_posterior);
disp('<strong>Rsquared:</strong>');
disp(glme_posterior.Rsquared);

[p,ci] = predict(glme_posterior,X);
X.N_fit = p;
X.CI_fit = ci;

tic; fprintf(1,'Merging output with original Spikes table...');
G = findgroups(X(:,{'TrialID','Probe','Channel'}));
tmp = splitapply(@(TrialID,Probe,Channel,N_fit,CI_fit,Exclude,N_Spline,N_Total_Expect,t)shrinkX2S(TrialID,Probe,Channel,N_fit,CI_fit,Exclude,N_Spline,N_Total_Expect,t),...
   X.TrialID,X.Probe,X.Channel,X.N_fit,X.CI_fit,X.t,X.Exclude,X.N_Spline,X.N_Total_Expect,G);
tmp = vertcat(tmp{:});

% Last: merge the original table with the output
S = outerjoin(S,tmp,...
   'Keys',{'TrialID','Probe','Channel'},...
   'MergeKeys',true,...
   'LeftVariables',{'TrialID','Probe','Channel','Trial','Reach','Grasp','Duration','Outcome','Forelimb','Spikes'},...
   'RightVariables',{'Exclude','N_Total_Expect','N_Spline','LB_fit','UB_fit','N_fit'});
S.N = sum(S.Spikes,2);
S.ChannelID = categorical(strcat(S.Probe,"::",S.Channel));

fprintf(1,'complete (%5.2f sec)\n',toc);
   function x = expandS2X(TrialID,Outcome,Duration,Probe,Channel,Spikes,t)
      
      N = nansum(Spikes,2); % Now "N" is like the other variables
      Spikes = Spikes'; % Flip so that time runs as rows
      t = t'; % Flip so that time runs as rows    
      nRep = size(Spikes,1); % Number of times to replicate everything

      x = table(repelem(TrialID,nRep,1),...
                 repelem(Outcome,nRep,1),...
                 repelem(Duration,nRep,1),...
                 repelem(Probe,nRep,1),...
                 repelem(Channel,nRep,1),...
                 repmat(t,size(Spikes,2),1),...
                 repelem(N,nRep,1),...
                 Spikes(:),...
            'VariableNames',{'TrialID',...
                             'Outcome',...
                             'Duration',...
                             'Probe',...
                             'Channel',...
                             't',...
                             'N',...
                             'Spikes'});
      x = {x};
      
   end

   function out = fitspline(Spikes,t,ChannelID,pars)
      ChannelID = ChannelID(1);
      [gAll,tAll] = findgroups(t);
      ymu = splitapply(@(n)nanmean(n),Spikes,gAll);
      [tAll,iSort] = sort(tAll,'ascend');
      tAll = reshape(tAll,1,numel(tAll));
      ymu = reshape(ymu(iSort),1,numel(ymu));
      
      ymu_to_knot = sgolayfilt(ymu,pars.SGOrder,pars.SGFrameLen);
      ymu_to_knot = [max(round(ymu_to_knot(1)),0),max(ymu_to_knot(2:(end-1)),0),max(round(ymu_to_knot(end)),0)];
      pars.Knots = union(pars.Knots,[1,numel(ymu_to_knot)]);

      pp = csape(tAll(pars.Knots),ymu_to_knot(pars.Knots),'clamped');
      N_Spline = fnval(pp,tAll)';
      t = tAll';
      ChannelID = repmat(ChannelID,numel(N_Spline),1);
      out = {table(ChannelID,t,N_Spline)};
      
   end

   function s = shrinkX2S(TrialID,Probe,Channel,N_fit,CI_fit,Exclude,N_Spline,N_Total_Expect,t)
      TrialID=TrialID(1);
      Probe=Probe(1);
      Channel=Channel(1);
      N_Spline = N_Spline.';
      N_fit = N_fit.';
      LB_fit = CI_fit(:,1).';
      UB_fit = CI_fit(:,2).';
      Exclude = Exclude(1);
      N_Total_Expect = N_Total_Expect(1);
      t = t.';
      s = table(TrialID,Probe,Channel,N_fit,LB_fit,UB_fit,Exclude,N_Spline,N_Total_Expect,t);
      
      s = {s}; % Return as cell array
   end

end