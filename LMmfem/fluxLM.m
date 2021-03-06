% FluxLM.m file for the approximate flux for LMmfem in the paper
% "Three Matlab Implementations of the Lowest-Order Raviart-Thomas
% MFEM with a Posteriori Error Control" by
% C.Bahriawati and C.Carstensen

function p=fluxLM(element,coordinate,u)
MidPoint=reshape(sum(reshape(coordinate(element',:),3,2*size(element,1)))...
                 ,size(element,1),2)/3;
P=reshape(u(1:3*size(element,1)),3,size(element,1));
p=reshape([ones(3,1)*P(1,:)+(ones(3,1)*P(3,:)).*...
         (reshape(coordinate(element',1),3,size(element,1))-...
          ones(3,1)*MidPoint(:,1)'),ones(3,1)*P(2,:)+(ones(3,1)*P(3,:)).*...
         (reshape(coordinate(element',2),3,size(element,1))-ones(3,1)*...
          MidPoint(:,2)')],3*size(element,1),2);

