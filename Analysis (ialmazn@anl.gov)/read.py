from readMDA import *

filename = '26idbSOFT_0077.mda'
d = readMDA(filename)
x = d[1].p[0].data 
y = d[2].p[0].data

#open fluorescence image
Zn_fluo = array(d[2].d[6].data)

figure()
imshow(Zn_fluo[::-1,:],extent=[min(min(y)),max(max(y)),min(x),max(x)])
colorbar()
axis('tight')
title(str(filename))

#open diffraction
diff = array(d[2].d[7].data)

figure()
imshow(log10(diff[::-1,:]),extent=[min(x),max(x),min(min(y)),max(max(y))])
colorbar()
axis('tight')
title(str(filename))

#find values
gamma = d[2].d[51].data
theta = d[2].d[52].data

#images of the night
scans = [75,77,83,85,87,91,93]

index = 0
figure(figsize=(15,25))
for ii in range(len(scans)):
    index += 1
    subplot(7,2,index)
    filename = '26idbSOFT_%04d.mda'%scans[ii]
    d = readMDA(filename)
    x = d[1].p[0].data 
    y = d[2].p[0].data
    theta = d[2].d[52].data
    imshow(Zn_fluo[::-1,:],extent=[min(min(y)),max(max(y)),min(x),max(x)])
    colorbar()
    axis('tight')
    index += 1
    title('Zn fluo: scan ='+str(scans[ii])+' - theta =%02.1f'%theta[0][0])
    subplot(7,2,index)
    diff = array(d[2].d[7].data)
    imshow(log10(diff[::-1,:]),extent=[min(min(y)),max(max(y)),min(x),max(x)])
    colorbar()
    axis('tight')
    title('Diff: scan ='+str(scans[ii])+' - theta =%02.1f'%theta[0][0])   

imgr = d[2].d[19].data

#convert tiff into h5 files
import matplotlib.pyplot as plt
tiff_file = '/mntdirect/_data_id01_inhouse/richard/aps/sample1/75/scan_75_img_028964.tif'
I = plt.imread(tiff_file)

import Image
raster = Image.open(tiff_file)
print raster.format, raster.size, raster.mode, raster.info
imarray=np.array(raster)

figure()
imshow(log10(imarray))

#convert into h5 file
import h5py
hf = h5py.File('data.h5', 'w')
hf.create_dataset('data', data=imarray)
hf.close()
dims = shape(imarray)



cd savedata/2018R1/20180306/Images/75h

data_h5 = h5py.File('test.h5','w')
data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
data_h5['/raw_data'][:,:] = imarray
data_h5.close()


#read h5 file
hf = h5py.File('data.h5', 'r')
figure()
imshow(log10(hf['data']))


hf = h5py.File('test.h5', 'r')
figure()
imshow(log10(hf['raw_data']))

cd savedata/2018R1/20180306/Images/75h

import os
import PIL.Image as Image
scan = 176
path = '/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'/'
l_file = os.listdir(path)

for ii in range(l_file):
    raster = Image.open(l_file[ii])
    imarray=np.array(raster)
    data_h5 = h5py.File('/home/beams25/USER26ID/savedata/2018R1/20180306/Images/'+str(scan)+'h/'os.path.splitext(l_file[ii])[0]+'.h5','w')
    data_h5.create_dataset('/raw_data',(dims[0],dims[1]),dtype='int16',maxshape=(dims[0],dims[1]),compression="gzip", compression_opts=9)
    data_h5['/raw_data'][:,:] = imarray
    data_h5.close()
 

#images of the night
scans = [213,215,217,219,221,223,225,227]

index = 0
fig=figure(figsize=(15,25))
for ii in range(len(scans)):
    index += 1
    subplot(8,2,index)
    filename = '26idbSOFT_%04d.mda'%scans[ii]
    d = readMDA(filename)
    x = d[1].p[0].data 
    y = d[2].p[0].data
    theta = d[2].d[52].data
    Zn_fluo = array(d[2].d[6].data)
    imshow(Zn_fluo[::-1,:],extent=[min(min(y)),max(max(y)),min(x),max(x)])
    hot()
    colorbar()
    axis('tight')
    index += 1
    title('Zn fluo: scan ='+str(scans[ii])+' - theta =%02.2f'%theta[0][0])
    subplot(8,2,index)
    diff = array(d[2].d[7].data)
    imshow(log10(diff[::-1,:]),extent=[min(min(y)),max(max(y)),min(x),max(x)])
    colorbar()
    axis('tight')
    hot()
    title('Diff: scan ='+str(scans[ii])+' - theta =%02.2f'%theta[0][0])   


#images of the night
scans = [176,178,180,182,184,186,188,190]

index = 0
fig=figure(figsize=(15,25))
for ii in range(len(scans)):
    index += 1
    subplot(8,2,index)
    filename = '26idbSOFT_%04d.mda'%scans[ii]
    d = readMDA(filename)
    x = d[1].p[0].data 
    y = d[2].p[0].data
    theta = d[2].d[52].data
    Zn_fluo = array(d[2].d[6].data)
    imshow(Zn_fluo[::-1,:],extent=[min(min(y)),max(max(y)),min(x),max(x)])
    hot()
    colorbar()
    axis('tight')
    index += 1
    title('Zn fluo: scan ='+str(scans[ii])+' - theta =%02.2f'%theta[0][0])
    subplot(8,2,index)
    diff = array(d[2].d[7].data)
    imshow(log10(diff[::-1,:]),extent=[min(min(y)),max(max(y)),min(x),max(x)])
    colorbar()
    axis('tight')
    hot()
    title('Diff: scan ='+str(scans[ii])+' - theta =%02.2f'%theta[0][0])   





