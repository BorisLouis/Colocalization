function [] = AlinePlot(selPos, IM, IM_DAPI, IM_marker, IM_nuc)

% Function to plot overlap between diffraction limited virus and DNA marker


figure
set(gcf, 'color', 'w');
pause
for i = 1:length (selPos.frame)
    fr = selPos.frame(i);
    %get contour nucleus from DAPI
    mask = bwboundaries(IM_nuc(:,:, fr));
    radi = 20;
    cod = [selPos.x(i) selPos.y(i)];
    
    %get image from Virus channel
    imVirus =  double(IM(:,:,fr));
    % reduce valus for better visualization
    imVirus(imVirus > max(imVirus(:))./5) = max(imVirus(:))./5;
    imVirus = (imVirus-min(imVirus(:)))./(max(max(imVirus-min(imVirus(:))))); %normalize
    
    %get image from Marker channel
    imMarker = double(IM_marker(:,:,fr));
    % reduce valus for better visualization
    imMarker(imMarker > max(imMarker(:))./3) = max(imMarker(:))./3;
    imMarker = (imMarker-min(imMarker(:)))./(max(max(imMarker-min(imMarker(:))))); %normalize
    %imMarker = H_norimage(imMarker);
    
    
    %plot virus image alone - red
    im = zeros([size(imVirus) 3]);
    im(:,:,2) = imVirus;
    subplot(2,3,1); 
    image(im); axis equal tight
    hold on
    plot(selPos.x(i), selPos.y(i), 'oy');
    for j = 1:length(mask)
        plot(mask{j}(:,2), mask{j}(:,1), 'c', 'Linewidth', 2);
    end
    hold off
    title('Virus')
    
    subplot(2,3,4); 
    image(im); axis equal tight
    axis([cod(1)-radi cod(1)+radi cod(2)-radi cod(2)+radi])
    title(['Virus #' num2str(selPos.particle(i))])
    
    %plot virus image alone - red
    im = zeros([size(imVirus) 3]);
    im(:,:,1) = imMarker;
    subplot(2,3,2); 
    image(im); axis equal tight
    hold on
    for j = 1:length(mask)
        plot(mask{j}(:,2), mask{j}(:,1), 'c', 'Linewidth', 2);
    end
    hold off
    title('DNA marker')
    
    subplot(2,3,5); 
    image(im); axis equal tight
    axis([cod(1)-radi cod(1)+radi cod(2)-radi cod(2)+radi])
    
    %plot overlay
    im = zeros([size(imVirus) 3]);
    im(:,:,1) = imMarker;
    im(:,:,2) = imVirus;
    subplot(2,3,6); 
    image(im); axis equal tight
    axis([cod(1)-radi cod(1)+radi cod(2)-radi cod(2)+radi])
    title(['pval(mean) = ' num2str(selPos.pval_mean(i))])
    
    %plot DAPI
    imDAPI = double(IM_DAPI(:,:,fr));
    subplot(2,3,3);
    imagesc(imDAPI);axis equal tight
    colormap gray
    hold on
    for j = 1:length(mask)
        plot(mask{j}(:,2), mask{j}(:,1), 'c', 'Linewidth', 2);
    end
    hold off
    title('DAPI')
    
    pause
    
end
    close(gcf)