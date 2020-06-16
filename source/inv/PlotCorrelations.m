function [Ax] = PlotCorrelations (m, mNames, mbnds)
% plots 1D and 2D marginal distributions of models in m
% 
% INPUTS
% m         matrix of models (Niter x Nvar)
% mNames    names of model parameters (1 x Nvar)
% bnds      bounds on parameters (Nvar x 2, or empty)

figure;

[Niter, Nvar] = size(m);
Ax = tight_subplot(Nvar,Nvar,0);
Ax = reshape(Ax, Nvar, Nvar)';

if isempty(mbnds)
    mbnds = [min(m)', max(m)'];
end

for icol = 1:Nvar
    for irow = 1:Nvar
        axes(Ax(irow,icol))
        if (icol == irow)
            histogram(m(:,icol), 20);
            xlim(mbnds(icol,:));
            set(gca,'YTickLabel',[]);
        elseif (icol<irow)
            Ax(irow,icol).Visible = 'off';
            continue;
        else
            dscatter(m(:,icol), m(:,irow), 'plottype', 'contour');
            axis([mbnds(icol,:), mbnds(irow,:)]);
            set(gca,'XTick',[],'XTickLabel',[]);
            
            if icol==Nvar
                set(gca,'YAxisLocation','right');
            else
                set(gca,'YTick',[],'YTickLabel',[]);
            end
        end
        
        if irow==1, title(mNames{icol}); end
    end
    
end




end