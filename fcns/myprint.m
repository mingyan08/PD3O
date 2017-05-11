function myprint(name, h)
% This code will produce eps and pdf files out of a MATLAB figure with no
% change. What you see will be what you get.
%
% Example:
%
%   h = figure(1);
%   ...(your plotting code) ...
%   myprint('filename', h);
%
% Your will get filename.eps and filename.pdf.
% -----------------
% Wotao Yin (2011).

if ~ischar(name)
    error('Input 1 must be char');
end

if ~exist('h','var') || isempty(h)
    h = gcf;
end

ppos = get(h,'Position'); ppos(1:2)=0; psize = ppos(3:4);
dpi = get(0,'ScreenPixelsPerInch');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', psize/dpi);
set(gcf, 'PaperPosition', ppos/dpi);

print('-painters','-depsc2','-r0',name);
print('-painters','-dpdf','-r0',name);

end

