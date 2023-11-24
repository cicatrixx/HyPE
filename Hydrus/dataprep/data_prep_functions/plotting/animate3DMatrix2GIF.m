%% animate layers in a 3D matrix
function animate3DMatrix2GIF(data_in,plotprefix, ptitle, gif_filename)
% filename = fullfile(figpath,'Q5kmVs500m_logval.gif'); % Specify the output file name
    
    nslides=size(data_in, 3);    
    %% Set up figs for making GIF
    data_max =  max(data_in,[],'all')
    data_min =  min(data_in,[],'all')
    fig=figure;%('Position', get(0, 'Screensize'));
    axes('Units', 'normalized', 'Position', [0 0 1 1])
    for m=1:nslides
        imagescnan(data_in(:,:,m));
        caxis([data_min, data_max])
        axis image
        title([plotprefix, ptitle{m}], Interpreter="none")
        colorbar
    
        drawnow
        frame = getframe(fig);
        im{m} = frame2im(frame);
    end
    
    %% animate plots
    delaytsec=1.25;
    for nn = 1:nslides
        [A,map] = rgb2ind(im{nn},256);
        if nn == 1
            imwrite(A,map,gif_filename,'gif','LoopCount',Inf,'DelayTime',delaytsec);
        else
            imwrite(A,map,gif_filename,'gif','WriteMode','append','DelayTime',delaytsec);
        end
    end