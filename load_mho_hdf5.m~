function [mho_data, mho_metadata] = load_mho_hdf5(file_path)
%% load_mho_hdf5 load the hdf5 file of the Millstone Hill IS radar measurements
% Created by Lei Cai on 13.06.2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%     fp:       file full path
% Output:
%     mho_data:  A 1-D struct array with the following fields:
%       .antenna
%       .pulse_code
%       .pulse_length
%       .variables:     a MATLAB struct. The field names indicate the different variables
%           Most frequently used variables are listed below, see also var_name_dict:
%             .AZ1:     Azimuth angle
%             .EL1:     Elevation angle
%             .n_e:     Electron density
%             .n_e_err: Electron density error
%             .v_i_los: lign-of-sight velocity (pos=away)
%             .v_i_los_err: lign-of-sight velocity error
%             .status: goodness of fit, NOTE: the defination is different from EISCAT status, see the note in mho_metadata.experiment_notes.
%     mho_metadata: A struct for the metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file_path = 'mlh160314g.005.hdf5';

var_name_dict = struct( ...
    'AZ1', 'az1', ...
    'AZ2', 'az2', ...
    'EL1', 'el1', ...
    'EL2', 'el2', ...
    'PULSE_LENGTH', 'pl', ...
    'T_SYS', 'systmp', ...
    'POWER_NORM', 'pnrmd', ...
    'P_Tx', 'power', ...
    'MODE_TYPE', 'mdtyp', ...
    'POWER_LENGTH_F', 'pulf', ...
    'LAG_SPACING', 'dtau', ...
    'IPP', 'ipp', ...
    'f_Tx', 'tfreq', ...
    'v_PHASE_Tx', 'vtx', ...
    'v_PHASE_Tx_err', 'dvtx', ...
    'SCAN_TYPE', 'scntyp', ...
    'CYCN', 'cycn', ...
    'POSN', 'posn', ...
    'RANGE_RES', 'mresl', ...
    'RANGE', 'range', ...
    'SNR', 'snp3', ...
    'RESIDUAL', 'wchsq', ...
    'STATUS', 'gfit', ...
    'FIT_TYPE', 'fit_type', ...
    'FPI_QUALITY', 'fpi_dataqual', ...
    'ACF_NORM', 'fa', ...
    'ACF_NORM_ERR', 'dfa', ...
    'n_pp', 'popl', ...
    'n_pp_err', 'dpopl', ...
    'n_e', 'ne', ...
    'n_e_err', 'dne', ...
    'T_i', 'ti', ...
    'T_i_err', 'dti', ...
    'T_r', 'tr', ...
    'T_r_err', 'dtr', ...
    'nu_i', 'co', ...
    'nu_i_err', 'dco', ...
    'v_i_los', 'vo', ...
    'v_i_los_err', 'dvo', ...
    'comp_H_p', 'ph+', ...
    'comp_H_p_err', 'dph+', ...
    'comp_mix', 'pm', ...
    'comp_mix_err', 'dpm', ...
    'v_DOP_los', 'vdopp', ...
    'v_DOP_los_err', 'dvdopp', ...
    'HEIGHT', 'gdalt' ...
);

field_names = fieldnames(var_name_dict);
for i=1 : length(field_names)
  queried_h5_var_names{i} = var_name_dict.(field_names{i});
end


% h5 = h5read(file_path);
data_groups = h5info(file_path, '/Data/Array Layout').Groups;
num_groups = length(data_groups);

mho_data = [];
for i=1 : num_groups
  group_name = data_groups(i).Name;
  
  rm = regexp(group_name, 'kinst=([+-]?\d+\.?\d*)', 'tokens');
  antenna_id = str2num(rm{1}{1});

  rm = regexp(group_name, 'mdtyp=([+-]?\d+\.?\d*)', 'tokens');
  pulse_code_id = str2num(rm{1}{1});

  rm = regexp(group_name, 'pl=([+-]?\d+\.?\d*)', 'tokens');
  pulse_code_length = str2num(rm{1}{1})*1e6; % in microseconds

  h5_vars = h5info(file_path, [group_name '/2D Parameters/']).Datasets;
  
  for j=1: length(h5_vars)
    h5_var_name = h5_vars(j).Name;

    h5_var_data = h5read(file_path, [group_name '/2D Parameters/' h5_var_name]);
    
    ind_var = find(ismember(queried_h5_var_names, h5_var_name));
    
    if isempty(ind_var)
      continue
    end


    var_name = field_names{ind_var};

    mho_data_tmp.variables.(var_name) = h5_var_data';
      
  end

  h5_vars = h5info(file_path, [group_name '/1D Parameters/']).Datasets;
  
  for j=1: length(h5_vars)
    h5_var_name = h5_vars(j).Name;

    h5_var_data = h5read(file_path, [group_name '/1D Parameters/' h5_var_name]);
    
    ind_var = find(ismember(queried_h5_var_names, h5_var_name));
    
    if isempty(ind_var)
      continue
    end


    var_name = field_names{ind_var};

    mho_data_tmp.variables.(var_name) = h5_var_data';
      
  end
  [nrow, ncol] = size(mho_data_tmp.variables.AZ1);
  mho_data_tmp.variables.RANGE = repmat(h5read(file_path, [group_name '/range']), 1, ncol);

  timestamps = h5read(file_path, [group_name '/timestamps']);
  dtl = datetime(timestamps,'ConvertFrom','epochtime','TicksPerSecond',1,'Format','dd-MMM-yyyy HH:mm:ss');
  mho_data_tmp.variables.timestamps = datenum(dtl)';

  switch antenna_id
    case 31
      antenna = 'misa';
    case 32
      antenna = 'zenith';
  end
  mho_data_tmp.antenna = antenna;
  
  switch pulse_code_id
    case 97
      pulse_code = 'alternating';
    case 98
      pulse_code = 'barker';
    case 115
      pulse_code = 'signle pulse';
  end
  mho_data_tmp.pulse_code = pulse_code;
  mho_data_tmp.pulse_length = pulse_code_length;

  mho_data = [mho_data mho_data_tmp];
end

mho_metadata.experiment_notes = h5read(file_path, '/Metadata/Experiment Notes').FileNotes';
mho_metadata.experiment_parameters.name = h5read(file_path, '/Metadata/Experiment Parameters').name';
mho_metadata.experiment_parameters.value = h5read(file_path, '/Metadata/Experiment Parameters').value';
mho_metadata.data_variables.name = h5read(file_path, '/Metadata/Data Parameters').mnemonic';
mho_metadata.data_variables.desc = h5read(file_path, '/Metadata/Data Parameters').description';
mho_metadata.data_variables.isError = h5read(file_path, '/Metadata/Data Parameters').isError;
mho_metadata.data_variables.unit = h5read(file_path, '/Metadata/Data Parameters').units';
mho_metadata.data_variables.kind = h5read(file_path, '/Metadata/Data Parameters').category';
end
