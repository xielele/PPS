function dmoo()    

global T T0 PROBLEM DIM 

DIM = 20;

hold off;

if isunix() > 0
    libs = sprintf('%s/include',pwd);
    loadlibrary(libs, 'setdmoo.h');
else
    loadlibrary('dset', 'include/setdmoo.h');
end
PROBLEM = calllib('dset', 'Problem');
TINTERV = calllib('dset', 'TInterval');
T       = calllib('dset', 'TInit');
TSTEP   = calllib('dset', 'TStep');
unloadlibrary('dset');

T0      = T;

if isunix() > 0
    libs = sprintf('%s/include',pwd);
    loadlibrary(libs, 'librun.h');
else
    loadlibrary('ddemo', 'include/librun.h');
end
s       = calllib('ddemo', 'Init', 'dmoo.set');
gen     = 0;
try
    while s >= 0 
        s = calllib('ddemo','Step',s);
                    
        UpdatePlot();

        tit = sprintf('Gen=%d,TW=%d',gen,T);
        title(tit);
        
        gen = gen + 1;
        if gen > 1 && mod(gen, TINTERV) == 0
            T = T + TSTEP;
        end        
        
        pause(0.01);
    end
catch
%    'there is an error!'
%     unloadlibrary('solverdmeda');
%     return;
    disp('program error!');
end

unloadlibrary('ddemo');
end

%% ========================================================================
function UpdatePlot()

% DRAWF = 0;

[pf, ps, xra, yra, zra, fra1, fra2, fra3] = Front();

fid = fopen( 'dmoo.pop' );
num = fscanf( fid, '%i', 4);
ds  = fscanf( fid, '%f', [num(3)+num(4), num(1)+num(2)] );
fclose(fid);
dim = num(3);
% if DRAWF >0
    %%% draw PF
    subplot(1,2,1);
    hold off;
    df  = ds(1:num(3),:);
    if dim == 2
        PPlot( pf, dim, 'b-',1 ); hold on;
    else
        PPlot( pf, dim, 'b.',1 ); hold on;
    end
    PPlot( df(:,1:num(1)), dim, 'rs',2 ); hold on;
    PPlot( df(:,num(1)+1:num(1)+num(2)), dim, 'co',2 );  hold on;
    xlabel('F1'); ylabel('F2'); xlim(fra1);  ylim(fra2);
    if dim == 3
        zlabel('F3'); zlim(fra3);
    end
% else
    %%% draw PS
    subplot(1,2,2); 
    hold off;
    if dim == 2
        PPlot( ps, dim, 'b-',1 ); hold on;
    else
        PPlot( ps, dim, 'b.',1 ); hold on;
    end    
    dx  = ds(num(3)+1:num(3)+num(4),:); 
    PPlot( dx(:,1:num(1)), dim, 'rs',2 ); hold on;
    PPlot( dx(:,num(1)+1:num(1)+num(2)), dim, 'cs',2 );  hold on; 
    xlabel('x1'); ylabel('x2'); xlim(xra);  ylim(yra);
    if dim == 3
        zlabel('x3'); zlim(zra);
    end
% end
end

%%
function [x1,x2] = mc(t, index)
global T0
% centre moving functions
switch index
    case 8
        x1 = 0;
        x2 = sin(0.5*t*pi);    
    case 9
        x1 = 0.5*t*pi;
        x2 = 0.5*t*pi;
    case {10,13,14}
        x1 = 0.5*t*pi;
        x2 = sin(0.5*t*pi);
    case 11
        r  = 0.5;
        x1 = r*cos(1*t*pi)+2*r*cos(0.5*t*pi);
        x2 = r*sin(1*t*pi)+2*r*sin(0.5*t*pi);    
    case 12
        x1 = 0.5*t*pi;
        x2 = sin(0.5*t*pi);
        if x2>=0
            x2 = x2-1.0;
        else
            x2 = x2+1.0;
        end        
end
end

%%
function [pf, ps, xrange, yrange, zrange, frange1, frange2, frange3] = Front()

global T PROBLEM DIM
    
frange1 = [0 1]; frange2 = [0 1]; frange3 = [0 1];
switch PROBLEM
    case 1 % FDA1
        G       = sin(0.5*pi*T);
        ps      = ones(2,50);
        ps(1,:) = linspace(0,1,50);
        ps(2,:) = ones(1,50)*G;
        pf      = zeros(2,50);
        pf(1,:) = linspace(0,1,50);
        pf(2,:) = 1.0-pf(1,:).^0.5;
        xrange  = [0,1];
        yrange  = [-1,1];
        zrange  = [0 0];
     case 2 % FDA2
        H       = 0.75 + 0.7*sin(0.5*pi*T);
        ps      = ones(2,50);
        ps(1,:) = linspace(0,1,50);
        ps(2,:) = zeros(1,50);
        pf      = zeros(2,50);
        pf(1,:) = linspace(0,1,50);
        pf(2,:) = 1.0-pf(1,:).^(1.0/H);
        xrange  = [0,1];
        yrange  = [-1,1];
        zrange  = [0 0];
      case 3 % FDA3
        G       = abs(sin(0.5*pi*T));
        ps      = ones(2,50);
        ps(1,:) = linspace(0,1,50);
        ps(2,:) = ones(1,50)*G;
        pf      = zeros(2,50);
        pf(1,:) = linspace(0,1,50);
        pf(2,:) = 1+G-(pf(1,:)*(1+G)).^0.5;
        xrange  = [0,1];
        yrange  = [-1,1];
        zrange  = [0 0];
        frange1 = [0 1]; frange2 = [0 2]; frange3 = [0 1];
      case 4 % FDA4
        G       = abs(sin(0.5*pi*T));
        [s,t]   = meshgrid(0:0.1:1,0:0.1:1);
        ps      = ones(3,121);
        ps(1,:) = reshape(s,[1,121]);
        ps(2,:) = reshape(t,[1,121]);
        ps(3,:) = ones(1,121)*G;
        pf      = zeros(3,121);
        pf(1,:) = cos(0.5*pi*ps(1,:)).*cos(0.5*pi*ps(2,:));
        pf(2,:) = cos(0.5*pi*ps(1,:)).*sin(0.5*pi*ps(2,:));        
        pf(3,:) = sin(0.5*pi*ps(1,:));
        xrange  = [0,1];
        yrange  = [0,1];
        zrange  = [0 1];  
    case 5 % DMOP1
        H       = 1.25+0.75*sin(0.5*pi*T);
        ps      = ones(2,50);
        ps(1,:) = linspace(0,1,50);
        ps(2,:) = zeros(1,50);
        pf      = zeros(2,50);
        pf(1,:) = linspace(0,1,50);
        pf(2,:) = 1.0-pf(1,:).^H;
        xrange  = [0,1];
        yrange  = [-1,1];
        zrange  = [0 0]; 
     case 6 % DMOP2
        H       = 1.25+0.75*sin(0.5*pi*T);
        G       = sin(0.5*pi*T);
        ps      = ones(2,50);
        ps(1,:) = linspace(0,1,50);
        ps(2,:) = ones(1,50)*G;
        pf      = zeros(2,50);
        pf(1,:) = linspace(0,1,50);
        pf(2,:) = 1.0-pf(1,:).^H;
        xrange  = [0,1];
        yrange  = [-1,1];
        zrange  = [0 0]; 
    case 7 % DMOP3
        G       = sin(0.5*pi*T);
        ps      = ones(2,50);
        ps(1,:) = linspace(0,1,50);
        ps(2,:) = ones(1,50)*G;
        pf      = zeros(2,50);
        pf(1,:) = linspace(0,1,50);
        pf(2,:) = 1.0-pf(1,:).^0.5;
        xrange  = [0,1];
        yrange  = [-1,1];
        zrange  = [0 0]; 
    case {8,9,10,11,12}
        H       = 1.25+0.75*sin(0.5*pi*T);
        [x1,y1] = mc(T,PROBLEM);
        ps      = ones(2,50);
        ps(1,:) = linspace(0,1,50);
        ps(2,:) = 1-ps(1,:).^H;
        ps(1,:) = ps(1,:)+x1;
        ps(2,:) = ps(2,:)+y1;
        pf      = zeros(2,50);
        pf(1,:) = linspace(0,1,50);
        pf(2,:) = 1.0-pf(1,:).^H;   
        xrange  = [x1,x1+1];
        yrange  = [y1,y1+1];
        zrange  = [0 0];
    case 13
        H       = 1.25+0.75*sin(0.5*pi*T);
        [x1,y1] = mc(T,PROBLEM);
        [s,t]   = meshgrid(0:0.1:1,0:0.1:1);
        ps      = ones(3,121);
        ps(1,:) = reshape(s,[1,121]);
        ps(2,:) = reshape(t,[1,121]);
        ps(3,:) = 1-ps(1,:).^H;
        pf      = zeros(3,121);
        pf(1,:) = cos(0.5*pi*ps(1,:)).*cos(0.5*pi*ps(2,:));
        pf(2,:) = cos(0.5*pi*ps(1,:)).*sin(0.5*pi*ps(2,:));        
        pf(3,:) = sin(0.5*pi*ps(1,:));
        ps(1,:) = ps(1,:) + x1;
        ps(2,:) = ps(2,:) + y1;
        ps(3,:) = ps(3,:) + y1;  
        xrange  = [x1,x1+1];
        yrange  = [y1,y1+1];
        zrange  = [y1,y1+1]; 
    case 14
        H       = 1.25+0.75*sin(0.5*pi*T);
        [x1,y1] = mc(T,PROBLEM);
        ps      = ones(2,50);
        ps(1,:) = linspace(0,1,50);
        if mod(floor(T*10.0),2)==0 
            ps(2,:) = 1-ps(1,:).^H;
        else
            ps(2,:) = ps(1,:).^H;
        end
        ps(1,:) = ps(1,:)+x1;
        ps(2,:) = ps(2,:)+y1;
        pf      = zeros(2,50);
        pf(1,:) = linspace(0,1,50);
        pf(2,:) = 1.0-pf(1,:).^H;   
        xrange  = [x1,x1+1];
        yrange  = [y1,y1+1];
        zrange  = [0 0];        
end
end

%%
function PPlot( data, dim, par, size )

switch dim
    case 2
        plot( data(1,:),data(2,:),par, 'MarkerSize',size);
    case 3
        plot3( data(1,:),data(2,:),data(3,:),par, 'MarkerSize',size);
end
end
