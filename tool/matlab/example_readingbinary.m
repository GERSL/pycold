% example for reading COLD result
path = '/Users/coloury/Dropbox/Documents/QTProjects/S-CCD/test/spectral_336_3980_obs_ccd.dat';
m = memmapfile(path, 'Format', {'int32',  [1 1] 't_start'; 'int32',  [1 1] 't_end'; 'int32', [1 1] 't_break';
    'int32',  [1 1] 'pos'; 'int32', [1 1] 'num_obs'; 'int32', [1 1] 'category'; 'int32', [1 1] 't_confirmed';
    'int32',  [1 1] 'change_prob'; 'double', [8 7]  'coefs'; 'double', [1 7]  'rmse';
    'double', [1 7]  'magnitude'});
rec_cg_c = transpose(m.Data)

% example for reading S-CCD result
path = '/Users/coloury/Dropbox/Documents/QTProjects/S-CCD/test/spectral_336_3980_obs_ccd.dat';
m = memmapfile(path, 'Format', {'int32',  [1 1] 't_start'; 'int32',  [1 1] 't_end'; 'int32', [1 1] 't_break';
    'int32',  [1 1] 'pos'; 'int32', [1 1] 'num_obs'; 'int16', [1 1] 'category'; 'int16', [1 1] 'land_type';
    'int32', [1 1] 't_confirmed'; 'int32',  [1 1] 'change_prob'; 'double', [6 7] 'coefs'; 'double', [1 7] 'rmse';
    'double', [1 7] 'magnitude'});
rec_cg_c_s = transpose(m.Data)
