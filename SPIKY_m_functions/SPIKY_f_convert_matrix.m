% If needed, this function converts the input spikes into the input format used by SPIKY,
% namely a cell with the spike times of the different spike trains as different cell arrays.
% It works both for padded matrices with spike times as rows as well as
% with matrices of zeros and ones with the ones indicating the spike times

function spikes=SPIKY_f_convert_matrix(matrix,dt,tmin)

if iscell(matrix)
    if length(matrix)>1
        num_trains=length(matrix);
        spikes=cell(1,num_trains);
        for trac=1:num_trains
            spikes{trac}=single(matrix{trac});
        end
    else
        num_trains=size(matrix{1},1);
        spikes=cell(1,num_trains);
        if isempty(setdiff(unique(reshape(matrix{1},1,numel(matrix{1}))),[0 1]))   % matrix with zeros and ones
            if nargin<3
                tmin=0;
            end
            if nargin<2
                dt=1;
            end
            for trac=1:num_trains
                spikes{trac}=single(tmin+find(matrix{1}(trac,:)>0)*dt);
            end
        else
            for trac=1:num_trains
                spikes{trac}=single(matrix{1}(trac,1:find(matrix{1}(trac,:),1,'last')));
            end
        end
    end
else
    num_trains=size(matrix,1);
    spikes=cell(1,num_trains);
    if isempty(setdiff(unique(reshape(matrix,1,numel(matrix))),[0 1]))   % matrix with zeros and ones
        if nargin<3
            tmin=0;
        end
        if nargin<2
            dt=1;
        end
        for trac=1:num_trains
            spikes{trac}=single(tmin+find(matrix(trac,:)>0)*dt);
        end
    else
        for trac=1:num_trains
            spikes{trac}=single(matrix(trac,1:find(matrix(trac,:),1,'last')));
        end        
    end
end

