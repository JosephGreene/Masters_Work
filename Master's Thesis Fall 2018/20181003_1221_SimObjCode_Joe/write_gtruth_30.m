clear
load('f_obj30.mat');
fn = 'gtruth30';
mkdir(fn);
h = figure(1);
for i = 1:20
    subplot(4,5,i);
    imagesc(f_obj30(:,:,i));
    axis image;
    title(sprintf('slice%d',i));
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
end

filename = sprintf('%s/slice%d.png',fn,i);
saveas(h,'gtruth30.png','png');