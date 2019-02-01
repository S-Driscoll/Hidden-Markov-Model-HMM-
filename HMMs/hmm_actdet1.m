function [ STATES1,STATES,STATESp,STATESn,sig_lev_rms ] = hmm_actdet1( sig_rms,nBins,sti,f,winN1)
    %%%% Seq derived from RMS
    clear sig_freq sig_bins sig_lev_rms
    flag =0;
    [sig_freq, sig_bins] = hist(sig_rms, nBins);
    % obtain the bin distance
    bin_dist = (sig_bins(2)-sig_bins(1));
    quant_error = 0;
    % indexing each float value to the nearest bin
    for bb = 1:length(sig_rms)
          for cc = 1:nBins
              dist = abs(sig_bins - sig_rms(bb));
              [minval, indx] = min(dist);
              quant_error = quant_error + (sig_rms(bb) - sig_bins(indx));
              sig_lev_rms(bb) = indx;
          end
    end
    clear TRGUESS EMITGUESS ESTTR ESTEMIT STATES  
    [TRGUESS,EMITGUESS] = hmmestimate(sig_lev_rms,sti); 
    [ESTTR,ESTEMIT] = hmmtrain(sig_lev_rms,TRGUESS,EMITGUESS);
    STATES = hmmviterbi(sig_lev_rms,ESTTR,ESTEMIT);
    STATES = STATES(:)-1;
    STATESp = find(diff(STATES)==1); STATESn = find(diff(STATES)==-1);
    %%%winN1 = round(1*f);
    if length(STATESp)~=length(STATESn)
        %figure;plot([sig_rms,1e-4*(sti),1e-4*(STATES+1.5)]);
        %STATESp,STATESn
        %dbstop in hmm_actdet1.m at 29
        flag=1;
    end
    if length(STATESp)<length(STATESn)
        STATESn(1)=[];
    elseif length(STATESp)>length(STATESn)
        STATESp(end)=[];
    end
    if length(STATESp)>1 && length(STATESn)>1
        ind1 = find((STATESp(2:end)-STATESn(1:end-1))<floor(0.15*f)); %%% Remove Spurious 0's in states
        STATESp(ind1+1)=[];STATESn(ind1)=[];
        flag=1;
    end
    
    ind = find((STATESn-STATESp)>winN1);
    
    if any(ind) %%% Remove Spurious 1's in states
        STATESp = STATESp(ind);
        STATESn = STATESn(ind);
    end
    
    if length(STATESp)>1 || length(STATESn)>1
        %[STATESp,STATESn]
        %dbstop in hmm_actdet1.m at 46
        flag =1;
        STATESp = STATESp(1);
        STATESn = STATESn(end);
        
    end
    STATES1=zeros(length(sig_rms),1);
    STATES1(STATESp:STATESn)=1;
%     if flag
%         figure;plot([sig_rms,1e-4*(sti),1.8e-4*(STATES1),1e-4*(STATES+1.5)]);legend('rms','sti','resti1','resti')
%         dbstop in hmm_actdet1.m at 62
%         
%     end
    
end

