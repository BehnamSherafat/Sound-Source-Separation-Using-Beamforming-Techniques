function [finalpos, azimuth_s, elevation_s, distance_s, finalsrp,finalfe, SearchVolume]=srplems(s, mic_loc, fs, lsb, usb, c)
%% This function uses SRP-PHAT with Stochastic Region Contraction (SRC) global-optimization algorithm to locate a single source using a frame of data and M microphones.
%%% Hoang Do - LEMS, Division of Engineering, Brown University.
%%% Paper on this algorithm can be found at:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H. Do, H. F. Silverman, and Y. Yu, “A real-time srp-phat source location implementation using stochastic region                  %
% contraction(src) on a large-aperture microphone array,” in Proc. IEEE Int. Conf. Acoust. Speech, Signal Process., Honolulu,      %
% HI, April 2007, vol. 1, pp. 121–124.                                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Note: This code is not fully optimized. Suggestions to improve the code is welcome.
%%% Please cite us if you use this work. Thank you.
%% Inputs:
%%% 1) s is "A FRAME" of data (L x M), L should be a power of 2
%%% 2) mic_loc is the microphone 3D-locations (M x 3) ( in meters)
%%% 3) fs: sampling rate (Hz)
%%% 4) lsb: a row-vector of the lower rectangular search boundary, e.g., [-2 -1 0] (meters)
%%% 5) usb: row-vector of the upper rectangular search boundary, e.g., [2 0 6] (m)
%%% It also calls other 2 subroutines: src and fe.
%% Outputs:
%%% 1) finalpos: estimated location of the source
%%% 2) finalsrp: srp-phat value of the point source
%%% 3) finalfe: number of fe's evaluated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Usage example:
%%% [finalpos,finalsrp,finalfe]=srplems(s, mic_loc, fs, lsb, usb)
%%% In our work, we use:
%%% 24 microphones ---> mic_loc is a (24x3) matrix.
%%% A 2048-sample framelength ---> s is a (2048 x 24) matrix
%%% 20 KHz sampling rate -> fs = 20,000.
%%% A search volume of 4m x 1m x 6m with the rectangular boundaries:
%%% lsb = [-2 -1 0] (m)
%%% usb = [2 0 6] (m)
%%% If the values of fs, lsb, usb are not specified in the inputs,
%%% fs=20,000Hz, lsb = [-2 -1 0], usb = [2 0 6] will be used.
%% Initialize variables:

warning off all


if nargin < 5, usb=[2 0 6]; end
if nargin < 4, lsb=[-2 -1 0]; end
if nargin < 3, fs=20000; end


L = size(s,1); %%% determine frame-length
M = size(mic_loc,1); %%% number of microphones (also number of channels of the input X).
np=M*(M-1)/2; %%%number of independent pairs

%steplength = L/4;  %%% non-overlapped length of a frame (75% overlap)
dftsize = 2^nextpow2(L);       %%% dft size equals to frame size
magiconst=10*fs/c;

%% Determine the maximum end-fire length (in samples) of a microphone pairs:
mdist=pdist(mic_loc);
efmax=max(mdist);
efs=2*(fix(efmax*fs/c)); %%% End-fire is symmetric about the 0th sample so multiplying by 2.
%efs=801;
hefs=round(efs/2); %%%half of the end-fire region
%% Doing the GCC-PHAT:
sf=fft(s,dftsize);                    %%%FFT of the original signals
csf=conj(sf);
yv=zeros(np,efs);                     %%%%Initialize yv to store the SRP-PHAT samples

p=1;
for i=1:M-1
    su1mic=sf(:,i)*ones(1,M);      %%%Create M copies of each signal's FFT
    prodall=su1mic.*csf;           %%%%Calculate the cross-power spectrum: = fft(x1).*conj(fft(x2))
    ss=prodall(:,i+1:M);           %%%% ss will be the cross-power spectra of microphone pairs (i,i+1), (i,i+2)...(i,M)
    ss=ss./(abs(ss)+eps);          %%%% PHAT weighting
    
    ssifft=real(ifft(ss,dftsize)); %%%% Get the GCC-PHAT ssifft, which is the IFFT of the cross-power spectra
    %newssifft=[ssifft(end-hefs:end,:); ssifft(1:efs-hefs-1,:)]; %%% Only select 'efs' samples (the beginning+end portions)
    newssifft=[ssifft(end-hefs+1:end,:); ssifft(1:efs-hefs,:)];
    newssifftr = newssifft';          %%%% Transpose it
    yv(p:p+M-i-1,:)=newssifftr;    %%%%Store in yv
    p=p+M-i;                       %%%% Update the current index of yv
    clear su1mic prodall ss ssifft newssifft newssifftr
end

%% Doing cubic-splines interpolation (factor of 10) on the GCC-PHAT:

xx=[1:.1:efs];
x=[1:efs];
yintp=spline(x,yv,xx);
yintpt=yintp';

%% Initialize to do SRC:

efsintp=length(xx)/2;

row1=([0:np-1]*2*efsintp)';

% randpts=150000;  %%% J0 in SRC. Depending on the size of the search volume, choose an appropriate value here (Here, 3000 is for a V_{search}= 4m x 1m x 6m)
% npoints=5000;   %%% Best N points. Again, choose an appropriate number according to your problem.
randpts=100000;
npoints=3000;

%%Doing SRC:
[bestpos, azimuth_s, elevation_s, distance_s, bestsrp,nofes, SearchVolume]=src(magiconst,lsb,usb,randpts,npoints,mic_loc,yintpt,row1,efsintp);

% Outputs:
finalpos=bestpos;  %%%Final source location estimate
finalsrp=bestsrp/np;  %%%Final source's SRP-PHAT value (normalized by number of pairs)
finalfe=nofes+randpts; %%% Number of fe's used.
        
%% SRC algorithm:


    function [bestpos, azimuth_s, elevation_s, distance_s, bestsrp,nofes, SearchVolume] = src(magiconst,bstart,bend,randpts,npoints,mic_loc,yintpt,row1,efsintp)
        %%%This function does SRC algorithm.
        %%%It calls subroutine "fe".
        
        if bstart==bend               %%%check if the lower search-boundary is the same as the upper one, i.e., V_{search}=0
            fprintf('Search volume is 0! Please expand it');
        end
        
        rng;
        %%Throw "randpts" random points in the defined boundary:
        for i=1:randpts
            [yval(i),position(i,:)] = fe(bstart,bend,magiconst,mic_loc,yintpt,row1,efsintp);
        end
        
        figure(11);
        view(43,24);
        for j = 1:M
            scatter3(mic_loc(j,1),mic_loc(j,2),mic_loc(j,3),'filled','MarkerFaceColor','r','MarkerEdgeColor','r');
            hold on
        end
        hold on
        ptCloud = pointCloud(position,'Intensity',yval');
        pcshow(ptCloud)
        colorbar
        
        newbstart=[min(position(:,1)) min(position(:,2)) min(position(:,3))];
        newbend=[max(position(:,1)) max(position(:,2)) max(position(:,3))];
        SearchVolume(1)=(newbend(1,1)-newbstart(1,1))*(newbend(1,2)-newbstart(1,2))*(newbend(1,3)-newbstart(1,3));
        
      
%         threshold = prctile(yval,perc_pts);
        idx_vec = find(yval > prctile(yval,98));
%         threshold = mean(yval);
%         idx_vec = find(yval >= threshold);
        srp_vec = yval(idx_vec);
        bestvec = position(idx_vec,:);
        
        ptCloud_in = pointCloud(bestvec,'Intensity',srp_vec');
        [ptCloudOut,~,~] = pcdenoise(ptCloud_in, 'NumNeighbors', 30, 'Threshold', 0.02);
        
        srp_vec = ptCloudOut.Intensity';
        bestvec = ptCloudOut.Location;
        
        position0=position;
        bestvec0=bestvec;
        
        figure(2);
        view(43,24);
        for j = 1:M
            scatter3(mic_loc(j,1),mic_loc(j,2),mic_loc(j,3),'filled','MarkerFaceColor','r','MarkerEdgeColor','r');
            hold on
        end
        hold on
        ptCloud2 = pointCloud(bestvec,'Intensity',srp_vec');
        pcshow(ptCloud2)
        colorbar
        hold on
        
        % finding best fit 3D line (direction of arrival)
        average_points = sum(srp_vec'.*bestvec,1)/sum(srp_vec',1);
        v1=[0,0,0];
        v2=average_points;
        v=[v2;v1];
        plot3(v(:,1),v(:,2),v(:,3),'b')
        
        % find azimuth and elevation of best fit line (direction of arrival)
        [azimuth,elevation,r] = cart2sph(average_points(1), average_points(2), average_points(3));
        azimuth_s(1) = azimuth*(180/pi); % azimuth is the counterclockwise angle in the x-y plane measured in radians from the positive x-axis. The value of the angle is in the range [-pi pi].
        elevation_s(1) = elevation*(180/pi); % elevation is the elevation angle in radians from the x-y plane. The value of the angle is in the range [-pi/2, pi/2].
        distance_s(1) = r;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%Define the new boundary vectors and search Volume:
        newbstart=[min(bestvec(:,1)) min(bestvec(:,2)) min(bestvec(:,3))];
        newbend=[max(bestvec(:,1)) max(bestvec(:,2)) max(bestvec(:,3))];
        SearchVolume(2)=(newbend(1,1)-newbstart(1,1))*(newbend(1,2)-newbstart(1,2))*(newbend(1,3)-newbstart(1,3));
        
        
        minVoxel = 1;  %%%the threshold (subject to change, depending on your data) of the sub-volume,i.e., a (cube-root of 100)-cm cube
%         maxFEs=3000000;      %%% threshold for the maximum number of functional evaluations (fe's) allowed to evaluate
        maxFEs = 100000; 
        
        TotalFEs=randpts;
        iter=1;
        countfe=0;
        kk=1;
        searchVol = SearchVolume(2);
        while searchVol>minVoxel && TotalFEs <=maxFEs
            
            threshold = mean(srp_vec);
            gidx = find(srp_vec > threshold);
            a5=srp_vec(gidx);
            newbestvec=bestvec(gidx,:); 
            la5=numel(a5);
            nonewpts(kk)=npoints-la5;   %%% the number of new points needed to evaluate
 
            ptCloud = pointCloud(newbestvec,'Intensity',a5');
            pcshow(ptCloud)
            colorbar
            
            figure(kk+2);
            view(43,24);
            for j = 1:M
                scatter3(mic_loc(j,1),mic_loc(j,2),mic_loc(j,3),'filled','MarkerFaceColor','r','MarkerEdgeColor','r');
                hold on
            end
            hold on
            ptCloud = pointCloud(newbestvec,'Intensity',a5');
            pcshow(ptCloud)
            colorbar
            hold on
                       
            % finding best fit 3D line (direction of arrival)
            average_points = sum(a5'.*newbestvec,1)/sum(a5',1);
            v1=[0,0,0];
            v2=average_points;
            v=[v2;v1];
            plot3(v(:,1),v(:,2),v(:,3),'b')

            % find azimuth and elevation of best fit line (direction of arrival)
            [azimuth,elevation,r] = cart2sph(average_points(1), average_points(2), average_points(3));
            azimuth_s(kk+1) = azimuth*(180/pi); % azimuth is the counterclockwise angle in the x-y plane measured in radians from the positive x-axis. The value of the angle is in the range [-pi pi].
            elevation_s(kk+1) = elevation*(180/pi); % elevation is the elevation angle in radians from the x-y plane. The value of the angle is in the range [-pi/2, pi/2].
            distance_s(kk+1) = r;
            
            newbstart=[min(newbestvec(:,1)) min(newbestvec(:,2)) min(newbestvec(:,3))];  %%% Contract search-volume to the volume containing only the high points
            newbend=[max(newbestvec(:,1)) max(newbestvec(:,2)) max(newbestvec(:,3))];

            if nonewpts(kk)>=1
                for i=1: nonewpts(kk)
                    a4=0;
                    while a4 < threshold && countfe <=maxFEs   %%%keep evaluating new random points until we fill up 'nonewpts' points that have SRP values above the mean
                        [a4,goodpos] = fe(newbstart,newbend,magiconst,mic_loc,yintpt,row1,efsintp);
                        countfe=countfe+1;
                        
                    end
                    
                    if a4~=0
                        a5=[a5 a4];
                        newbestvec=[newbestvec; goodpos];
                    end
                    
                end
            end
            
            srp_vec=a5;
            clear bestvec
            bestvec=newbestvec;
            clear a5 newbestvec a4 goodpos newbstart newbend;
            
            searchVol=(max(bestvec(:,1))-min(bestvec(:,1)))*(max(bestvec(:,2))-min(bestvec(:,2)))*(max(bestvec(:,3))-min(bestvec(:,3))); %%% The new search-volume containing only the best N points
            SearchVolume(kk+2) = searchVol;
            TotalFEs=randpts+countfe;
            kk=kk+1;
        end
        
        [v,i]=max(srp_vec);
        bestsrp=v;
        bestpos=bestvec(i,:);
        nofes=TotalFEs;
        clear bestvec
        clear a5 newbestvec a4 goodpos newbstart newbend;
        clear searchVol srp_vec
    end

%% Functional Evaluation sub-routine (calculate SRP-PHAT value for a point in the search space):

    function [yval1,position1] = fe(bstart,bend,magiconst,mic_loc,yintpt,row1,efsintp)
        %%%This function evaluates an 'fe' in the search space, gives back the 'fe' position and its SRP-PHAT value.
        
        
        M=size(mic_loc,1); %%%number of mics
        
        %%%generate a random point in the search space:
        y=rand(1,3);
        x=(bstart.*(1-y)+bend.*(y));   %This vector x defines the coordinates (x,y,z) of the rand point x
        
        
        
        %%%Find the distances from M-microphones to the point:
        a1=ones(M,1);
        xx1=a1*x;
        xdiff1=xx1-mic_loc;
        dists=sqrt(sum(xdiff1.*xdiff1,2));
        
        
        %%%%Differences in distances:
        ddiffs_ones=ones(M,1)*dists';
        ddm=ddiffs_ones-ddiffs_ones';
        
        %%% Calculate the TDOA index:
        %  v=ddm(gidM);
        v=nonzeros(tril(ddm,0));
        
        %ddiffsi32=int32(round(magiconst*v+efsintp));
        %ddiffsi32=floor(magiconst*v+efsintp)+1;
        ddiffsi32=round(magiconst*v+efsintp+1.5);
        row=row1+ddiffsi32;
        
        %%%Pull out the GCC-PHAT value corresponding to that TDOA:
        v1=yintpt(row);
        %%% SRP-PHAT value of the point:
        yval1=sum(v1);
        position1=x;
    end
end
