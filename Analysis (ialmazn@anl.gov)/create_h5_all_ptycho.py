##convert all files to h5 files
##use anaconda for python
#/APSshare/anaconda/x86_64/bin/python /home/beams25/USER26ID/savedata/2018R1/20180306/Analysis/create_h5_all.py
#/APSshare/anaconda/x86_64/bin/python /home/beams25/USER26ID/savedata/2018R1/20180306/Analysis/create_h5_all_ptycho.py

##transfert to my account
#scp -r -P 5022 429h mrichard@firewall.esrf.fr:/data/id01/inhouse/richard/aps/data/

import os
import PIL.Image as Image
import numpy as np
import h5py
import numpy


scans = numpy.arange(13,17,1)

for jj in range(len(scans)):
    path = '/home/beams25/USER26ID/savedata/2018R1/20180312/Images/'+str(scans[jj])+'/'
    l_file = os.listdir(path)
    if (os.path.isdir('/home/beams25/USER26ID/savedata/2018R1/20180312/Images/'+str(scans[jj])+'h/')):
        print('directory exists')
    else:  
        os.mkdir('/home/beams25/USER26ID/savedata/2018R1/20180312/Images/'+str(scans[jj])+'h/')
        for ii in range(len(l_file)):
            raster = Image.open(path+l_file[ii])
            imarray=np.array(raster)
            data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180312/Images/'+str(scans[jj])+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
            dims = np.shape(imarray)
            data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
            data_h5['/raw_data'][:,:] = imarray
            data_h5.close()
