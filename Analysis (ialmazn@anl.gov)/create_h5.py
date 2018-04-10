##use anaconda for python
#/APSshare/anaconda/x86_64/bin/python /home/beams25/USER26ID/savedata/2018R1/20180306/Analysis/create_h5.py
##transfert to my account
#scp -r -P 5022 429h mrichard@firewall.esrf.fr:/data/id01/inhouse/richard/aps/sample3/

import os
import PIL.Image as Image
import numpy as np
import h5py
import numpy


scan = 598
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)
os.mkdir('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/')

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    dims = np.shape(imarray)
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

'''
scan = 426
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    dims = np.shape(imarray)
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()
 
'''
'''
scan = 444
path = '/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    dims = shape(imarray)
    data_h5 = h5py.File('/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

scan = 446
path = '/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    dims = shape(imarray)
    data_h5 = h5py.File('/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

scan = 448
path = '/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    dims = shape(imarray)
    data_h5 = h5py.File('/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

scan = 450
path = '/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    dims = shape(imarray)
    data_h5 = h5py.File('/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

scan = 452
path = '/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    dims = shape(imarray)
    data_h5 = h5py.File('/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

scan = 454
path = '/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    dims = shape(imarray)
    data_h5 = h5py.File('/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

scan = 456
path = '/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    dims = shape(imarray)
    data_h5 = h5py.File('/mntdirect/_data_id01_inhouse/richard/aps/sample2/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()


'''

'''

scan = 444
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    dims = np.shape(imarray)
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

scan = 446
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    dims = np.shape(imarray)
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

scan = 448
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    dims = np.shape(imarray)
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

scan = 450
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    dims = np.shape(imarray)
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()


scan = 452
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    dims = np.shape(imarray)
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

scan = 454
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    dims = np.shape(imarray)
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

scan = 456
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    dims = np.shape(imarray)
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

'''
'''

scan = 458
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    dims = np.shape(imarray)
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()


scan = 460
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    dims = np.shape(imarray)
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

scan = 462
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    dims = np.shape(imarray)
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

scan = 464
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(len(l_file)):
    raster = Image.open(path+l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'+os.path.splitext(l_file[ii])[0]+'.h5','w')
    dims = np.shape(imarray)
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()

'''



