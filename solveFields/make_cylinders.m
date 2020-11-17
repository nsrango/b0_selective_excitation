prefix = '/home/narango/Dropbox (MIT)/ultimate_shim/basis';

deg_ords = floor(linspace(5,33,(35-5)/2));
    if ~exist(prefix, 'dir')
      mkdir(prefix)
    end
for i = 1:length(deg_ords)


    nx = 111;
    ny = 111;
    nz = 40;
    coordsx = linspace(-.55,.55,nx);
    coordsy = linspace(-.55,.55,ny);
    coordsz = linspace(-.4,.4,nz);

    [X,Y,Z] = meshgrid(coordsx,coordsy,coordsz);

    rhoReference = sqrt(2)*.55*1.1;
    degreeMax = deg_ords(i);
    orderMax = deg_ords(i);
    rk = gpuArray(createTargetPointCylinder(rhoReference,4*rhoReference,degreeMax,orderMax));
    figure(1);
    plot3(rk(:,1),rk(:,2),rk(:,3),'*');


    bfull = surface_basis(rk,X,Y,Z);
    save(fullfile(prefix,sprintf('cylinder_o%d_rsmall.mat',deg_ords(i))),'bfull','-v7.3','-nocompression');
end



