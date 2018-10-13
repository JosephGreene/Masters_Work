clear
load('f_obj5.mat');
fn = 'gtruth5';
mkdir(fn);
for i = 1:5
    h = imagesc(f_obj5(:,:,i));
    axis image;
    filename = sprintf('%s/slice%d.png',fn,i);
    saveas(h,filename,'png');
end