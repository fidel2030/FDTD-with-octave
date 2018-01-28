function ret = fdtdplot1d(ER, Ey, Hx, dz)
  ymaxlim = max([ max(Hx) max(Ey)]);
  yminlim = min([ min(Hx) min(Ey)]);
  ylim([yminlim, ymaxlim]);
  za=[0:length(Ey)-1]*dz;
  plot(za, Ey, 'b','linewidth', 2, za, Hx, 'r', 'linewidth', 2);% linewidth=2);
end
