function [traj,rad_ang] = radial_trajectory(kdim,goldenangle,varargin)
% Analytical description of the k-space trajectory for radial sampling.
% Either calculated with BART or custom code. Image dimensions are also
% encoded in the trajectory.
%
% Note:
% - BART doesnt work with tiny golden angle yet
% - Always assumes similar partitions for stack-of-stars

% Check input
kdim=c12d(kdim);

% BART switch
if ~isempty(varargin)
    isbart=1;
else
    isbart=0;
end

if isbart
    traj_call=['traj -r -x ',num2str(kdim(1)),' -y ',num2str(kdim(2))];
    
    if goldenangle > 0
        traj_call=[traj_call,' -G'];
    end
    
    traj=bart(traj_call);
else
    % Pre-allocate trajectory matrix
    traj=zeros(3,kdim(1),kdim(2),prod(kdim([3 5:end])));
    
    % Get radial angles for uniform (rev) or golden angle
    if goldenangle > 0
        d_ang=(pi/(((1+sqrt(5))/2)+goldenangle-1));
    else
        d_ang=pi/(kdim(2));
    end
    rad_ang=0:d_ang:d_ang*(kdim(2)-1);
    
    % Line reversal for uniform
    if goldenangle == 0
        rad_ang(2:2:end)=rad_ang(2:2:end)+pi;
        rad_ang=mod(rad_ang,2*pi);
    end
    
    % Calculate samples for single spoke
    %kx=linspace(0,2*kdim(1)-1,kdim(1)+1)'-(2*kdim(1)-1)/2;kx(end)=[];
    kx=linspace(0,kdim(1)-1,kdim(1)+1)'-(kdim(1)-1)/2;kx(end)=[];
    
    % Modulate successive spokes
    for ky=1:kdim(2)
        traj(1,:,ky,:)=repmat(kx*exp(1j*rad_ang(ky)),[1 1 size(traj,4)]);
    end
    
    % Reshape to image size
    traj=reshape(traj,[3 kdim([1:3 5:end])]);
    
    % Split in channels
    traj(2,:)=real(traj(1,:));
    traj(1,:)=imag(traj(1,:));
    
    % Simple z for stack-of-stars (linear increment)
    if kdim(3) > 1
        if mod(kdim(3),2)==0 % is_odd
            kz=linspace(-kdim(3)/2,kdim(3)/2,kdim(3)+1);kz(end)=[];
        else
            kz=linspace(-kdim(3)/2,kdim(3)/2,kdim(3));
        end
        for z=1:size(kdim(3))
            traj(3,:,:,:,:,:,:,:,:,:)=repmat(permute(kz,[1 3 4 2]),[1 kdim(1:2) 1 kdim(5:end)]);
        end
    end
end

disp('+Radial analytical trajectory is calculated.')
% END
end