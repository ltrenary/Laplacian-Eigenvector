function [lp] = laplacianevec_field(lon,lat, mask ,neofs,field)
%%% PURPOSE: 
%This code finds the eigenvectors of the Laplacian operator and the component time series for user provided data. 
%
%%%% INPUT
% lon = vector of longitude (SX x 1), where SX = number of longitude
% lat = vector of latitude (SY x 1), where SY = number of latitude
% mask = logical matrix containing 1 and 0's that isolates the regions
% of interst. Desired region is listed using 1's.
% neofs = number of eigenvectors to retain.
% field (optional) = data to be projected onto the eigenvectors of the Laplacian
% operator to yield time series. This is an optional input.
%
%%% OUTPUT
% lp.eof = the eigenvectors of the Laplacian operator
% lp.pc (optional) = the time series recovered from projecting user specified data
% onto the eigenvectors.
% lp.fexpvar = the fraction of explained variance per eigenvector.
% Code written by L. Trenary (5/2017).

if (nargin <5)
    field =0.0;
end


% Just making sure data dimensions are consistent
mask = permute(mask,[2 1]);
[dim1, dim2] = size(mask);
[dim1f,dim2f,dim3f]=size(field);

if (dim1 == dim2f)
    field = permute(field,[2 1 3]);
end
%% DEFINE GRID
 nlon=length(lon);
  nlat=length(lat);

  x=double(lon)/360*pi*2;
  y=double(lat)/360*pi*2;
  dth = abs(y(2)-y(1));
  dph = abs(x(2)-x(1));

  
  idum=0; 
  clear theta phi
  
  for jj=1:nlat
    for ii=1:nlon
      if (mask(ii,jj));
        idum=idum+1;
        theta(idum) = y(jj);
        phi(idum)  =  x(ii);
      end
    end
  end
  n=length(theta);

  if nnz(cos(theta)<0)
    disp('abs(theta) > pi/2!')
    return
  end

%% Defining the laplace operator

K = zeros(n,n);
B = zeros(n,1);

%Calculate the greens function. See eqns. (13), (15), (20), and (22) in DelSole and Tippet (2015;
%J.CLIM).

 for i1=1:n
    i2=1:i1;
    term1 = 2*sin( (theta(i1)-theta(i2))/2 ).^2;
    term2 = 2*cos(theta(i1))*cos(theta(i2));
    term3 = sin( (phi(i1)-phi(i2))/2 ).^2;
    G = -1/(4*pi)*log(term1+term2.*term3);
    K(i1,i2) = G.*sqrt(cos(theta(i1))*cos(theta(i2)))*dth*dph;
    K(i2,i1) = K(i1,i2);
    r = sqrt(dth*dph*cos(theta(i1))/pi);
    K(i1,i1) = r^2*(1-2*log(r/sqrt(2)))/4;
    B(i1) = sqrt(cos(theta(i1)));
  end
  beta = 1./sqrt(sum(B.^2));
  B = beta*B;
%% COMPUTE ORTHOGONAL COMPLEMENT TO CONSTANT VECTOR. Refer to eqn. (24) in DelSole and Tippett (2015, J. CLIM).
u = B.*ones(n,1);
u = u/norm(u);


K = (K - u*(u'*K) - (K*u)*u' + u*(u'*K*u)*u');

%% FIND THE EIGENVECTORS OF LAPLACIAN
%[Z,d] = eigs(K,neofs);
[Z,d,v]=svd(K);

Bi = 1./B;
YY = NaN(nlat*nlon, neofs+1);

% Eigenfunctions 
  YY(mask(:),1) = Bi.*u;
  YY(mask(:),2:(neofs+1)) = Bi(:,ones(1,neofs)).*Z(:,[1:neofs]);

% Psuedo Inverse (need to scale fields by the cos of latitud).
YYi = NaN(nlat*nlon, neofs+1);
  YYi(mask(:),1) = B.*u;
  YYi(mask(:),2:(neofs+1)) = B(:,ones(1,neofs)).*Z(:,[1:neofs]);


%% Project onto data
area = B.^2;


lp.eofi = reshape(YYi,nlon,nlat,neofs+1);
lp.eof = reshape(YY,nlon,nlat,neofs+1);

if (nargin == 5)
[ny,nx,ntime] = size(field);
fld = reshape(field,[ny*nx, ntime]);
lp.pc =  fld(mask(:),:)'*YYi(mask(:),:);
% Refer to (34) in DelSole and Tippett (2015, J. CLIM).
u = B.*ones(n,1);
for np =1:neofs
    at = nansum(lp.pc(:,np).^2);
    bt = sum((fld(mask(:),:,:).^2)'*area);
lp.fexpvar(np,1) =(at/bt)*100.;
end


else
    return
end


end
