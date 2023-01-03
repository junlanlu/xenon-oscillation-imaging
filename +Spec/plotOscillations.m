function plotOscillations(dyn, BHs, detrend, fitted, save_fig_flag, save_fig_path)

[BHstart, BHend] = findBHs(dyn.t(:,1), BHs);

if length(fitted.freq) == length(dyn.t(BHstart:BHend))
    fitType = 'sine';
else 
    fitType = 'peaks';
end 

% Colors
colors = [0.8500    0.3250    0.0980 
          0.4660    0.6740    0.1880
          0         0.4470    0.7410];

switch fitType
case 'sine'
    figure(2), clf
    set(gcf, 'Position',[875    -200    550    725])
    axFontFs = 13;

    subplot(4,1,1), hold on
    plot(dyn.t(BHstart:BHend),detrend.area*100,'k--')
    plot(dyn.t(BHstart:BHend),fitted.area*100,'Linewidth',2,'color',colors(1,:)); 
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim([dyn.t(BHstart) dyn.t(BHend)])
    xlabel('Time (s)'), ylabel('Amplitude')
    set(gca,'FontSize',axFontFs)
    ylim([-10 10])

    subplot(4,1,2), hold on
    plot(dyn.t(BHstart:BHend),detrend.freq,'k--')
    plot(dyn.t(BHstart:BHend),fitted.freq,'Linewidth',2,'color',colors(1,:)); 
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim([dyn.t(BHstart) dyn.t(BHend)])
    xlabel('Time (s)'), ylabel('Chemical Shift (ppm)')
    set(gca,'FontSize',axFontFs)
    ylim([-.5 .5])

    subplot(4,1,3), hold on
    plot(dyn.t(BHstart:BHend),detrend.fwhm,'k--')
    plot(dyn.t(BHstart:BHend),fitted.fwhm,'Linewidth',2,'color',colors(1,:)); 
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim([dyn.t(BHstart) dyn.t(BHend)])
    xlabel('Time (s)'), ylabel('FWHM (ppm)')
    set(gca,'FontSize',axFontFs)
    ylim([-.5 .5])

    subplot(4,1,4), hold on
    plot(dyn.t(BHstart:BHend),detrend.phase,'k--')
    plot(dyn.t(BHstart:BHend),fitted.phase,'Linewidth',2,'color',colors(1,:)); 
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim([dyn.t(BHstart) dyn.t(BHend)])
    xlabel('Time (s)'), ylabel('Phase (degrees)')
    set(gca,'FontSize',axFontFs)
    ylim([-10 10])
case 'peaks'
        figure(2), clf
    set(gcf, 'Position',[875    -200    550    725])
    axFontFs = 13;

    subplot(4,1,1), hold on
    plot(dyn.t(BHstart:BHend),detrend.area*100,'Linewidth',2,'color',colors(1,:)); 
    plot(fitted.at,fitted.area*100,'.k','MarkerSize',30);
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim([dyn.t(BHstart) dyn.t(BHend)])
    xlabel('Time (s)'), ylabel('Amplitude')
    set(gca,'FontSize',axFontFs)
    ylim([-10 10])

    subplot(4,1,2), hold on
    plot(dyn.t(BHstart:BHend),detrend.freq,'Linewidth',2,'color',colors(1,:)); 
    plot(fitted.ft,fitted.freq,'.k','MarkerSize',30);
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim([dyn.t(BHstart) dyn.t(BHend)])
    xlabel('Time (s)'), ylabel('Chemical Shift (ppm)')
    set(gca,'FontSize',axFontFs)
    ylim([-.5 .5])

    subplot(4,1,3), hold on
    plot(dyn.t(BHstart:BHend),detrend.fwhm,'Linewidth',2,'color',colors(1,:)); 
    plot(fitted.fwt,fitted.fwhm,'.k','MarkerSize',30); 
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim([dyn.t(BHstart) dyn.t(BHend)])
    xlabel('Time (s)'), ylabel('FWHM (ppm)')
    set(gca,'FontSize',axFontFs)
    ylim([-.5 .5])

    subplot(4,1,4), hold on
    plot(dyn.t(BHstart:BHend),detrend.phase,'Linewidth',2,'color',colors(1,:)); 
    plot(fitted.pt,fitted.phase,'.k','MarkerSize',30);
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim([dyn.t(BHstart) dyn.t(BHend)])
    xlabel('Time (s)'), ylabel('Phase (degrees)')
    set(gca,'FontSize',axFontFs)
    ylim([-10 10])
end 

if save_fig_flag == 1
    options.Format = 'tiff';
    hgexport(gcf,save_fig_path,options)
end  
