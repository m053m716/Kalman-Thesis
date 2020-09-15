function [fig,CData] = animateReach(P,trialID,varToAnimate,varargin)
%ANIMATEREACH Animate 3D points of markers for reach
%
%  fig = preview.animateReach(P,trialID,varToAnimate);
%
% Inputs
%  P - Kinematic table after cleaning
%  trialID - Trial to animate
%  varToAnimate - Variable name to animate
%
% Output
%  fig - Figure handle

pars = struct;
pars.C = cfg.gfx('Color_Struct');
pars.Fig = [];
pars.FrameRateApparent = 1/15; % This is only for the (unsaved) animation loop
pars.NFramesPause = 10;
fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

if nargin < 2
   trialID = randi(max(P.TrialID),1);
end

if nargin < 3
   varToAnimate = 'X';
end
P = P(P.TrialID==trialID,:);
X = P.(varToAnimate);
t = P.Properties.UserData.t;
G = findgroups(P.Marker);

if isempty(pars.Fig)
   fig = figure(...
      'Name','Reach Animation',...
      'Units','Normalized',...
      'Color','w',...
      'Position',[0.2 0.2 0.6 0.6]);
   xl = [min(P.(varToAnimate)(1:3:end,:),[],'all'),max(P.(varToAnimate)(1:3:end,:),[],'all')];
   dx = diff(xl); 
   xl = [xl(1) - 0.05*dx, xl(2) + 0.05*dx];
   yl = [min(P.(varToAnimate)(2:3:end,:),[],'all'),max(P.(varToAnimate)(2:3:end,:),[],'all')];
   dy = diff(yl); 
   yl = [yl(1) - 0.05*dy, yl(2) + 0.05*dy];
   zl = [min(P.(varToAnimate)(3:3:end,:),[],'all'),max(P.(varToAnimate)(3:3:end,:),[],'all')];
   dz = diff(zl); 
   zl = [zl(1) - 0.05*dz, zl(2) + 0.05*dz];
   ax = axes(fig,...
      'XColor','k','YColor','k',...
      'View',[-37.5 30],...
      'NextPlot','add',...
      'YDir','reverse','XDir','reverse','ZDir','reverse',...
      'XLim',xl,'YLim',yl,'ZLim',zl);
   box(ax,'on');
   grid(ax,'on');
   uM = unique(P.Marker);
   h = gobjects(numel(uM),1);
   for iH = 1:numel(h)
      h(iH) = line(ax,nan(1,5),nan(1,5),nan(1,5),'LineWidth',2,'Color',pars.C.(uM(iH)));
   end

else
   fig = pars.Fig;
   ax = findobj(fig.Children,'-depth',0,'Type','axes');
   ax = ax(1); % Only want 1
   h = get(ax,'Children');
end

xlabel(ax,'X','FontName','Arial','Color','k');
ylabel(ax,'Y','FontName','Arial','Color','k');
zlabel(ax,'Z','FontName','Arial','Color','k');
set(get(ax,'Title'),'FontName','Arial','Color','k');
if strcmpi(varToAnimate,'X')
   ttxt = 'observed';
else
   ttxt = strrep(varToAnimate(2:end),'_',' ');
end

[~,iPause] = min(abs(t));

if nargout > 1
   CData = getframe(fig);
   CData = repmat(CData,1,1,numel(t)+pars.NFramesPause);
   ii = 0;
   for iT = 1:numel(t)
      data = splitapply(@(X)getTimeSample(X,t(iT),t),X,G);
      arrayfun(@(h,data)cycleData(h,data),h,data);  
      title(ax,sprintf('Grasp (%s) %+5.2f ms',ttxt,t(iT)*1e3));
      pause(pars.FrameRateApparent);
      drawnow;
      CData(:,:,iT+ii) = getframe(fig);
      if iT==iPause
         for ii = 1:pars.NFramesPause
            pause(pars.FrameRateApparent);
            CData(:,:,iT+ii) = getframe(fig);
         end
      end
   end
   
else
   CData = [];
   for iT = 1:numel(t)
      data = splitapply(@(X)getTimeSample(X,t(iT),t),X,G);
      arrayfun(@(h,data)cycleData(h,data),h,data);  
      title(ax,sprintf('Grasp (%s) %+5.2f ms',ttxt,t(iT)*1e3));
      pause(pars.FrameRateApparent);
      drawnow;
      if iT==iPause
         for ii = 1:pars.NFramesPause
            pause(pars.FrameRateApparent);
         end
      end
   end
   
end


   function xt = getTimeSample(X,tSel,tRef)
      xt = {X(:,tRef==tSel)};
      % Returns XYZ as cell array, in order
   end

   function cycleData(h,data)
      h.XData = [data{1}(1), h.XData(1:(end-1))];
      h.YData = [data{1}(2), h.YData(1:(end-1))];
      h.ZData = [data{1}(3), h.ZData(1:(end-1))];
   end

end