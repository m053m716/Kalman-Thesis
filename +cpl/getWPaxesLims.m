function [XL,YL,ZL] = getWPaxesLims(worldPoints)
%% GETWPAXESLIMS  Get axes limits for all world points representing reach trajectories
%
%  [XL,YL,ZL] = GETWPAXESLIMS(worldPoints);
%
% By: Max Murphy  v1.0  11/20/2018  Original version (R2017b)

%%
XL = [inf,-inf]; 
YL = [inf,-inf]; 
ZL = [inf,-inf];
for k = 1:size(worldPoints,1)
   for m = 1:size(worldPoints{k,1},1)
      XL = [min(XL(1),min(worldPoints{k,1}{m,1}(:,1))), ...
            max(XL(2),max(worldPoints{k,1}{m,1}(:,1)))];
      YL = [min(YL(1),min(worldPoints{k,1}{m,1}(:,2))), ...
            max(YL(2),max(worldPoints{k,1}{m,1}(:,2)))];
      ZL = [min(ZL(1),min(worldPoints{k,1}{m,1}(:,3))), ...
            max(ZL(2),max(worldPoints{k,1}{m,1}(:,3)))];
   end
end

end