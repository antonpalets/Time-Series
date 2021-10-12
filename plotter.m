function plotter(data, m, maOrder)
    if nargin<2
        m=25;
    end
    if nargin<3
        maOrder=0;
    end
    
    figure
    set(gcf, 'Position',  [100, 400, 1500, 400])
    subplot(131)
    %stem(0:m,  acf(data, m))
    acf(data, m, 0.05, 1, maOrder, 1);
    
    subplot(132)
    %stem(0:m, pacf(data, m))
    pacf(data, m, 0.05, 1, 1);
    
    subplot(133)
    normplot(data)
    