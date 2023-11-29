function test()    

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

gen     = 0;
for s = 1:250 

    UpdatePlot(s);

    tit = sprintf('Gen=%d,TW=%d',gen,T);
    title(tit);

    gen = gen + 1;
    if gen > 1 && mod(gen, TINTERV) == 0
        T = T + TSTEP;
    end        

    pause(0.01);
end
end

%% ========================================================================
function UpdatePlot(s)

% DRAWF = 0;

[pf, ps, xlow, xupp, ylow, yupp] = Front();

snam= sprintf('data\\%d.pop',s);
fid = fopen( snam );
num = fscanf( fid, '%i', 4);
ds  = fscanf( fid, '%f', [num(3)+num(4), num(1)+num(2)] );
fclose(fid);
% if DRAWF >0
    %%% draw PF
    subplot(1,2,1);
    hold off;
    df  = ds(1:num(3),:);
    PPlot( pf, 2, 'b-',1 ); hold on;
    PPlot( df(:,1:num(1)), 2, 'rs',2 ); hold on;
    PPlot( df(:,num(1)+1:num(1)+num(2)), 2, 'co',2 );  hold on;
    xlabel('F1'); ylabel('F2');
    xlim([0 1]);  ylim([0 1]);
% else
    %%% draw PS
    subplot(1,2,2); 
    hold off;
    PPlot( ps, 2, 'b-',1 ); hold on;
    dx  = ds(num(3)+1:num(3)+num(4),:); 
    PPlot( dx(:,1:num(1)), 2, 'rs',2 ); hold on;
    PPlot( dx(:,num(1)+1:num(1)+num(2)), 2, 'cs',2 );  hold on; 
    xlabel('x1'); ylabel('x2');
    xlim([xlow xupp]);  ylim([ylow yupp]);
% end
end

%%
function [x1,x2] = mc(t, index)
global T0
% centre moving functions
switch index
    case 2
        x1 = 0;
        x2 = sin(0.5*t*pi);    
    case 3
        x1 = 0.5*t*pi;
        x2 = 0.5*t*pi;
    case 4
        x1 = 0.5*t*pi;
        x2 = sin(0.5*t*pi);
    case 5
        x1 = 0.5*t*pi;
        x2 = sin(0.5*t*pi);
%         x2 = x2 + (x2>=0)*1;
%         x2 = x2 - (x2< 0)*1;
        if mod(0.25*(t-T0),2)<1
            x2 = x2+1.0;
        else
            x2 = x2-1.0;
        end
    case 6
        r  = 0.5;
        x1 = r*cos(1*t*pi)+2*r*cos(0.5*t*pi);
        x2 = r*sin(1*t*pi)+2*r*sin(0.5*t*pi);      
end
end

%%
function x2 = ms(x1, t)
% shape moving functions

global DIM

h  = 1.5+sin(0.5*t*pi);% - 1.5 + 1 + 2.0*2.0/DIM;
x2 = x1.^h;
end

%%
function [pf, ps, x1, x2, y1, y2] = Front()

global T PROBLEM
    
switch PROBLEM
    case 1 % F1
        G       = sin(0.5*pi*T);
        ps      = ones(2,50);
        ps(1,:) = linspace(0,1,50);
        ps(2,:) = ones(1,50)*G;
        pf      = zeros(2,50);
        pf(1,:) = linspace(0,1,50);
        pf(2,:) = 1.0-pf(1,:).^0.5;
        x1      = 0;
        x2      = 1;
        y1      = -1;
        y2      = 1;
%     case 2 % F2
%         G       = sin(0.5*pi*T);
%         H       = 1.5+G;
%         ps      = ones(2,50);
%         ps(1,:) = linspace(0,1,50);
%         ps(2,:) = ps(1,:).^H - G;     
%         pf      = zeros(2,50);
%         pf(1,:) = linspace(0,1,50);
%         pf(2,:) = 1.0-pf(1,:).^H; 
%         x1      = 0;
%         x2      = 1;
%         y1      = -1;
%         y2      = 2;          
    otherwise %F3-F6
        H       = 1.5+sin(0.5*pi*T);
        [x1,y1] = mc(T,PROBLEM);
        ps      = ones(2,50);
        ps(1,:) = linspace(0,1,50);
        ps(2,:) = 1-ms(ps(1,:),T);
        ps(1,:) = ps(1,:)+x1;
        ps(2,:) = ps(2,:)+y1;
        pf      = zeros(2,50);
        pf(1,:) = linspace(0,1,50);
        pf(2,:) = 1.0-pf(1,:).^H;   
        x2      = x1 + 1;
        y2      = y1 + 1;
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
