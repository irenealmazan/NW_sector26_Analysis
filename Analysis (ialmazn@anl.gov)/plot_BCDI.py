import sys, os, time, re
import scipy, numpy, pylab
from scipy.optimize import leastsq
import scipy.signal
import h5py
sys.path.append('/mntdirect/_users/mrichard/test/Jan/macro_Dom/lib/python2.7/site-packages')
import xray as xu
import cPickle as pickle
import multiprocessing
import scipy.misc
os.environ["PYOPENCL_COMPILER_OUTPUT"] = "0" #set to 1 to see the compilation going on
from silx.image import sift
from numpy import median
try:
    from PyMca5.PyMca import ArraySave
    from PyMca5.PyMca import EdfFile
except:
    from PyMca import ArraySave
    from PyMca import EdfFile
import PyMca5

# not sorted never import *
from pylab import *
from scipy import *
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib import rc, font_manager
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show
from readMDA import *



cdict = {'red':  ((0.0, 1.0, 1.0),
                  (0.11, 0.0, 0.0),
                  (0.36, 0.0, 0.0),
                  (0.62, 1.0, 1.0),
                  (0.87, 1.0, 1.0),
                  (1.0, 0.0, 0.0)),
          'green': ((0.0, 1.0, 1.0),
                  (0.11, 0.0, 0.0),
                  (0.36, 1.0, 1.0),
                  (0.62, 1.0, 1.0),
                  (0.87, 0.0, 0.0),
                  (1.0, 0.0, 0.0)),
          'blue': ((0.0, 1.0, 1.0),
                  (0.11, 1.0, 1.0),
                  (0.36, 1.0, 1.0),
                  (0.62, 0.0, 0.0),
                  (0.87, 0.0, 0.0),
                  (1.0, 0.0, 0.0))}
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

filename = '26idbSOFT_0423.mda'
d = readMDA(filename)

ccdn = array(d[2].d[19].data)
fluo = array(d[2].d[6].data)

figure()
imshow(fluo)
axis('tight')

dims = shape(fluo)
x = arange(0,dims[1],1)
max_fluo = zeros((dims[0]))
for ii in range(dims[0]):
    max_fluo[ii] = find(fluo[ii,:]==fluo[ii,:].max())[0]

imgr = zeros((dims[0]))
for ii in range(dims[0]):
    imgr[ii] = ccdn[ii,max_fluo[ii]] 

scan_no = 423
datadir = '/mntdirect/_data_id01_inhouse/richard/aps/sample3'
for idx in range(len(imgr)):
    i = int(imgr[idx]) # read ccd image from EDF file
    img_name = 'scan_%03d_img_%05d.h5'%(scan_no,i)
    hdf5file = os.path.join(datadir,'%03d'%(scan_no)+'h/',img_name)
    if not os.path.isfile(hdf5file):
        print "no file"
        break
    else:
        h = h5py.File(hdf5file)
        print hdf5file
        ccd = h['raw_data'][:,:]
        CCD = xu.blockAverage2D(ccd,3,3)#, nav[0],nav[1], roi=roi)
        if idx==0:
            intensity = numpy.zeros( (len(ccdn),) + CCD.shape )
        intensity[idx,:,:] = scipy.signal.medfilt(CCD,[3,3])  


figure()
imshow(xu.maplog(intensity.sum(axis=0),5,0),cmap=my_cmap)

figure()
imshow(xu.maplog(intensity.sum(axis=1),5,0),cmap=my_cmap)

figure()
imshow(xu.maplog(intensity.sum(axis=2),5,0),cmap=my_cmap)

###second scan

filename = '26idbSOFT_0426.mda'
d = readMDA(filename)

ccdn = array(d[2].d[19].data)
fluo = array(d[2].d[6].data)

figure()
imshow(fluo)
axis('tight')

dims = shape(fluo)
x = arange(0,dims[1],1)
max_fluo = zeros((dims[0]))
for ii in range(dims[0]):
    max_fluo[ii] = find(fluo[ii,:]==fluo[ii,:].max())[0]

imgr = zeros((dims[0]))
for ii in range(dims[0]):
    imgr[ii] = ccdn[ii,max_fluo[ii]] 

scan_no = 426
datadir = '/mntdirect/_data_id01_inhouse/richard/aps/sample3'
for idx in range(len(imgr)):
    i = int(imgr[idx]) # read ccd image from EDF file
    img_name = 'scan_%03d_img_%05d.h5'%(scan_no,i)
    hdf5file = os.path.join(datadir,'%03d'%(scan_no)+'h/',img_name)
    if not os.path.isfile(hdf5file):
        print "no file"
        break
    else:
        h = h5py.File(hdf5file)
        print hdf5file
        ccd = h['raw_data'][:,:]
        CCD = xu.blockAverage2D(ccd,3,3)#, nav[0],nav[1], roi=roi)
        if idx==0:
            intensity = numpy.zeros( (len(ccdn),) + CCD.shape )
        intensity[idx,:,:] = scipy.signal.medfilt(CCD,[3,3])  



figure()
imshow(xu.maplog(intensity.sum(axis=0),5,0),cmap=my_cmap)

figure()
imshow(xu.maplog(intensity.sum(axis=1),5,0),cmap=my_cmap)

figure()
imshow(xu.maplog(intensity.sum(axis=2),5,0),cmap=my_cmap)



#create in 3D
nx,ny,nz = shape(intensity)
QX,QY,QZ = np.mgrid[0:nx:1j * nx,
                      0:ny:1j * ny,
                      0:nz:1j*nz]
INT = xu.maplog(intensity,4,0)

from mayavi import mlab
mlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0))
mlab.contour3d(QX, QY, QZ, INT, contours=8, opacity=0.5)
mlab.colorbar(title="log(int)", orientation="vertical")



#/APSshare/anaconda/x86_64/bin/python /home/beams25/USER26ID/savedata/2018R1/20180306/Images/create_vtk.py


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





