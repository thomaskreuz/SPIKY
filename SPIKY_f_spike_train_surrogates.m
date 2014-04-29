% This function generates various kinds of spike train surrogates (which differ in the
% selection of properties that are maintained). Currently the properties maintained 
% are either (1) the individual spike numbers, or (2) the individual interspike interval
% distribution (in addition) or (3) the pooled spike train or (4) the peri-stimulus time histogram (PSTH).

function surro_spikes=SPIKY_f_spike_train_surrogates(spikes,para)

num_trains=length(spikes);
num_spikes=cellfun('length',spikes);
surro_spikes=cell(1,num_trains);
if para.choice==1 % keep number of spikes for each spike train constant
   for trac=1:num_trains
       surro_spikes{trac}=sort(para.tmin+rand(1,num_spikes(trac))*(para.tmax-para.tmin));
       surro_spikes{trac}=single(unique(round(surro_spikes{trac}/para.dts)*para.dts));
   end
elseif para.choice==2 % keep interspike interval distribution for each spike train constant
   for trac=1:num_trains
       isi=diff([para.tmin spikes{trac} para.tmax]);
       surro_spikes{trac}=cumsum(isi(randperm(length(isi))));
       surro_spikes{trac}=single(unique(round(surro_spikes{trac}(1:end-1)/para.dts)*para.dts));
   end
elseif para.choice==3 % keep pooled spike train
   all_spikes=[spikes{:}];
   num_all_spikes=length(all_spikes);
   ori_labels=zeros(1,num_all_spikes);
   sc=0;
   for trac=1:num_trains
       ori_labels(sc+(1:num_spikes(trac)))=trac;
       sc=sc+num_spikes(trac);
   end
   surro_labels=ori_labels(randperm(num_all_spikes));
   for trac=1:num_trains
       surro_spikes{trac}=all_spikes(surro_labels==trac);
   end
elseif para.choice==4 % keep PSTH
   % PSTH
   num_bins = 100;
   bins = min([spikes{:}]):(max([spikes{:}])-min([spikes{:}]))/num_bins:max([spikes{:}]);
   psth = histc([spikes{:}], bins);
   
   % ICDF
   lambda = 2;
   ipsth = icdf(psth, 2);
   
   %Generating (still working on)
end
