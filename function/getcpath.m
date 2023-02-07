function folder = getcpath()

folder = which('NEURALTBX_C_Folder');
folder = folder(1:end-20);