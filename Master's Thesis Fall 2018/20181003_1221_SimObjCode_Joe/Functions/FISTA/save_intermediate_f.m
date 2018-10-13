% make output directories for results
born_order = forwardObj.bornorder;
f = xhat;
nz = forwardObj.sigSize(3);
iters_time = toc;

f_outdir = sprintf('%s/born%d/iters_%d',results_dir,born_order,iIter);
mkdir(f_outdir);
mkdir([f_outdir,'/f_abs']);
mkdir([f_outdir,'/f_real']);
mkdir([f_outdir,'/f_imag']);

%%% f-reconstruction result - saving

% Save f as mat file

fn =[f_outdir,'/f_born',num2str(born_order),'.mat'];

save(fn,'f','cost_all','cost_fid','cost_reg','iters_time');

% Save abs(f) as image files
fn =[f_outdir,'/f_abs'];
f_min= min(abs(f(:)));
f_max = max(abs(f(:)));
for m =1 : nz
    slice = mat2gray(abs(f(:,:,m)),[f_min,f_max]);
    filename = [fn,'/TV_born',num2str(born_order),'_abs_', num2str(m),'.bmp'];
    imwrite(slice,filename,'bmp')
end

% Save real(f) as image files
fn =[f_outdir,'/f_real'];
f_min= min(real(f(:)));
f_max = max(real(f(:)));
for m =1 : nz
    slice = mat2gray(real(f(:,:,m)),[f_min,f_max]);
    filename = [fn,'/TV_born',num2str(born_order),'_real_', num2str(m),'.bmp'];
    imwrite(slice,filename,'bmp')
end

% Save imag(f) as image files
fn =[f_outdir,'/f_imag'];
f_min= min(imag(f(:)));
f_max = max(imag(f(:)));
for m =1 : nz
    slice = mat2gray(imag(f(:,:,m)),[f_min,f_max]);
    filename = [fn,'/TV_born',num2str(born_order),'_imag_', num2str(m),'.bmp'];
    imwrite(slice,filename,'bmp')
end