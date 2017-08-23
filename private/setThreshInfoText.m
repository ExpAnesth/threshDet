function h=setThreshInfoText(baselinePar,thresh,h)

% clear any previous value
if any(ishandle(h))
  delete(h);
end
if all(isfinite(baselinePar))
  % distance of thresh from base line, expressed in terms of noise
  % variability specified in baselinePar
  relativeThresh=(thresh-baselinePar(1))/baselinePar(2);
  h=ultext(['thresh absolute: ' num2str(thresh,3)...
    '; iqhd_{16-84}: ' num2str(relativeThresh,3)],0.005,...
    'color','y','fontsize',10,'fontweight','bold');
end