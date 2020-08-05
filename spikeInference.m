oasis_setup;
spikes_OASIS = zeros(size(DFoF));
denoised_Ca_OASIS = zeros(size(DFoF));
clear options_OASIS;
temp_denoised_Ca = {};
temp_spikes = {};
for i=1:size(sequences,1)
	parfor k=1:size(DFoF,1) % if using parfor, subscript should be i+offset, not sequences(i,1).
        [temp_denoised_Ca{k,i}, temp_spikes{k,i}, options_OASIS(k,i)] = ...
		deconvolveCa(DFoF(k,sequences(i,1):sequences(i,2)), 'constrained', 'ar1', 'optimize_pars', 'optimize_b', 'tau_range', [10,60]);
	end
	disp(['Temporal segment ',num2str(i),' finished.']);
end
for i=1:size(sequences,1)
    for k=1:size(DFoF,1)
        denoised_Ca_OASIS(k,sequences(i,1):sequences(i,2)) = temp_denoised_Ca{k,i};
        spikes_OASIS(k,sequences(i,1):sequences(i,2)) = temp_spikes{k,i};
    end
end
clear temp_denoised_Ca;
clear temp_spikes;


%% calculate SNR to filter out too noisy traces
RSS = zeros(1,size(DFoF,1));
signal = zeros(1,size(DFoF,1));
for k=1:size(DFoF,1)
	temp_b = zeros(1,size(DFoF,2));
	for ii=1:size(sequences,1)
		temp_b(sequences(ii,1):sequences(ii,2)) = options_OASIS(k,ii).b*ones(1,sequences(ii,2) - sequences(ii,1) + 1);
	end
	RSS(k) = std(temp_b + denoised_Ca_OASIS(k,:) - DFoF(k,:));
	temp_spikes = spikes_OASIS(k,:);
	temp_denoised_Ca = denoised_Ca_OASIS(k,:);
	temp_denoised_Ca(temp_spikes<=0) = [];
	signal(k) = mean(temp_denoised_Ca);
	if isempty(temp_denoised_Ca)
		signal(k) = 0;
	end
end
SNR = signal./RSS;
clear temp_b;
clear temp_spikes;
clear temp_denoised_Ca;


%% display
temp = randperm(size(DFoF,1));
figure('Name','deconvolution by OASIS');
for i=1:12
	subplot(4,3,i);
	k = temp(i);
	plot(DFoF(k,:));
	hold on;
	plot(denoised_Ca_OASIS(k,:));
	plot(spikes_OASIS(k,:));
	temp_b = zeros(1,size(DFoF,2));
	temp_sn = zeros(1,size(DFoF,2));
	for ii=1:size(sequences,1)
		temp_b(sequences(ii,1):sequences(ii,2)) = options_OASIS(k,ii).b*ones(1,sequences(ii,2) - sequences(ii,1) + 1);
		temp_sn(sequences(ii,1):sequences(ii,2)) = options_OASIS(k,ii).sn*ones(1,sequences(ii,2) - sequences(ii,1) + 1);
	end
	plot(temp_b);
	plot(temp_sn);
	volume = length(find(A3(:,k)>0));
	title([int2str(k),'(volume: ',num2str(volume),') (SNR: ',num2str(SNR(k)),')']);
end
legend('DFoF','denoised Ca','spikes','b','sn');
clear temp_b;
clear temp_sn;

figure('Name','deconvolution by OASIS');
for i=1:6
	subplot(3,2,i);
	k = temp1(i);
	plot(DFoF(k,:));
	hold on;
	plot(denoised_Ca_OASIS(k,:));
	plot(spikes_OASIS(k,:));
	temp_b = zeros(1,size(DFoF,2));
	temp_sn = zeros(1,size(DFoF,2));
	for ii=1:size(sequences,1)
		temp_b(sequences(ii,1):sequences(ii,2)) = options_OASIS(k,ii).b*ones(1,sequences(ii,2) - sequences(ii,1) + 1);
		temp_sn(sequences(ii,1):sequences(ii,2)) = options_OASIS(k,ii).sn*ones(1,sequences(ii,2) - sequences(ii,1) + 1);
	end
	plot(temp_b);
	plot(temp_sn);
	volume = length(find(A3(:,k)>0));
	title([int2str(k),'(volume: ',num2str(volume),') (SNR: ',num2str(SNR(k)),')']);
end
legend('DFoF','denoised Ca','spikes','b','sn');
clear temp_b;
clear temp_sn;

figure('Name','deconvolution by OASIS');
for i=1:6
	subplot(3,2,i);
	k = temp3(i);
	plot(DFoF(k,:));
	hold on;
	plot(denoised_Ca_OASIS(k,:));
	plot(spikes_OASIS(k,:));
	temp_b = zeros(1,size(DFoF,2));
	temp_sn = zeros(1,size(DFoF,2));
	for ii=1:size(sequences,1)
		temp_b(sequences(ii,1):sequences(ii,2)) = options_OASIS(k,ii).b*ones(1,sequences(ii,2) - sequences(ii,1) + 1);
		temp_sn(sequences(ii,1):sequences(ii,2)) = options_OASIS(k,ii).sn*ones(1,sequences(ii,2) - sequences(ii,1) + 1);
	end
	plot(temp_b);
	plot(temp_sn);
	volume = length(find(A3(:,k)>0));
	title([int2str(k),'(volume: ',num2str(volume),') (SNR: ',num2str(SNR(k)),')']);
end
legend('DFoF','denoised Ca','spikes','b','sn');
clear temp_b;
clear temp_sn;