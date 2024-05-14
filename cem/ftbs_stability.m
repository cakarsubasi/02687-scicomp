% Visually check stability of the ftbs scheme.
g_rtl = @(C, thet) abs(1 - C.*(exp(1j*thet) - 1));
g_ltr = @(C, thet) abs(1 - C.*(1 - exp(-1j*thet)));

c = linspace(-1, 1, 100);
thet = linspace(0, 2*pi, 100);
[C, THET] = meshgrid(c, thet);
z_rtl = g_rtl(C, THET);
z_ltr = g_ltr(C, THET);

figure;
contourf(c, thet, z_rtl);
colorbar;
shading flat;
%colormap lines;
xlabel("ak/h");
ylabel("\theta")
exportgraphics(gcf, "ex3_ftbs_stability1.pdf");
figure;
contourf(c, thet, z_ltr);
colorbar;
shading flat;
%colormap lines;
xlabel("ak/h");
ylabel("\theta")
exportgraphics(gcf, "ex3_ftbs_stability2.pdf");