function out_struct = MetaResponseAnal_VmResponseWin(expt,input_struct)
%input_struct for this needs
%-sigdata_filt which is the average across trials for a given stimulus
%-sigon/sigoff
%-basetimes
out_struct.response_vm_confint = [];
out_struct.spont_var = [];
out_struct.up_win = [];
out_struct.low_win = [];
out_struct.allVmNotVec = [];
out_struct.allVmNotVar = [];
out_struct.vm_mean_up = [];
out_struct.vm_max_up = [];
out_struct.vm_var_up = [];
out_struct.vm_mean_low = [];
out_struct.vm_min_low = [];
out_struct.vm_var_low = [];
out_struct.allVmRespVec_low = [];
out_struct.allVmRespVec_up = [];

allfields = fieldnames(input_struct);
for ifield = 1:size(allfields,1)
    s = [allfields{ifield} ' = input_struct.' allfields{ifield} ';'];
    eval(s)
end

conf_p = 98;
windowsize = 200;
binsize = 10;
bin_p = 0.85;

out_struct.response_vm_confint = getCDFconf (mean(sigdata_filt(:,basetimes(1):sigon)),conf_p);
out_struct.spont_var = mean(var(sigdata_filt(:,basetimes(1):sigon)));

[up_inds, low_inds] = WindowResponse(mean(sigdata_filt(:,sigon:sigoff)), out_struct.response_vm_confint, windowsize, binsize, bin_p);
out_struct.up_win = getWindowEdges (up_inds, 1, 1)+sigon;
out_struct.low_win = getWindowEdges (low_inds, 1, 1)+sigon;

allresp = union(up_inds,low_inds);
allind = [1:size(sigdata_filt,2)];
allind(allresp) = 0;
notind = find(allind ~= 0);
out_struct.allVmNotVec = mean(sigdata_filt(:,notind));
out_struct.allVmNotVar = var(sigdata_filt(:,notind));

allVmRespVec_up = [];
allVmRespVec_low=[];
respind = 1;
for iresp = 1:size(out_struct.up_win,2)
    %     spkinds = [];
    thisresp = mean(sigdata_filt(:,out_struct.up_win(1,iresp):out_struct.up_win(2,iresp)));
    out_struct.vm_mean_up (respind) = mean(thisresp);
    out_struct.vm_max_up (respind) = max(thisresp);
    out_struct.vm_var_up (respind) = var(thisresp);
    allVmRespVec_up = [allVmRespVec_up, thisresp];
    %     for itrial=1:size(spiketimes,2)
    %         these_spk = spiketimes{itrial};
    %         spkinds = [spkinds, intersect(these_spk, [up_win(1,iresp):up_win(2,iresp)])];
    %     end
    %     out_struct.response_spk (respind) = (size(spkinds,2)/size(spiketimes,2))/(diff(up_win(:,iresp))*expt.wc.dt);
    respind = respind + 1;
end

respind = 1;
for iresp = 1:size(out_struct.low_win,2)
    %     spkinds = [];
    thisresp = mean(sigdata_filt(:,out_struct.low_win(1,iresp):out_struct.low_win(2,iresp)));
    
    out_struct.vm_mean_low (respind) = mean(thisresp);
    out_struct.vm_min_low (respind) = min(thisresp);
    out_struct.vm_var_low (respind) = var(thisresp);
    allVmRespVec_low = [allVmRespVec_low, thisresp];
    respind = respind + 1;
end
out_struct.allVmRespVec_low = allVmRespVec_low;
out_struct.allVmRespVec_up = allVmRespVec_up;