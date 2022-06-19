# Ventricle Volume developed by Arash Rabbani, rabarash@yahoo.com
# This py file is responsible for semantic segmentation of the cine cardiac MRI data into left ventricle, right ventricle, and myocardium
# After running this file, please run Main.m to generate the final outputs
import numpy as np
from tensorflow import keras
import cv2
import scipy.io as sio
from skimage.measure import label
from skimage.measure import regionprops

def flatout(A,prc=None): 
    # replaces outliers with .01 and 0.99 percentiles
    if prc==None:
        prc=.01
    p1=np.quantile(A,prc)
    p2=np.quantile(A,1-prc)
    A[A>p2]=p2
    A[A<p1]=p1
    return A
def normal(A):
    # normalizes an array
    A_min = np.min(A)
    return (A-A_min)/(np.max(A)-A_min)
def prep(A,ss=128):
    # preprocessing of images to size 128 x 128
    A=np.squeeze(A)
    S=A.shape;
    if S[0]<S[1]:
        m=round((S[1]-S[0])/2); x1=1; x2=S[0]; y1=m+1; y2=m+S[0];
    if S[0]>=S[1]:
        m=round((S[0]-S[1])/2); x1=m+1; x2=m+S[1]; y1=1; y2=S[1];
    A=A[x1:x2,y1:y2]
    S=A.shape
    m=round(S[0]*.1);
    A=A[m:-m,m:-m];
    A=flatout(A,.01)
    Rat=A.shape[0]/ss/.9
    A=cv2.resize(A, (ss,ss),interpolation=cv2.INTER_LINEAR)
    A=normal(A)
    return A,Rat
def remisol(AG,bkvalue,typ='SA'):
    # removes isolated objects
    for J in range(AG.shape[0]):
        AGNew=AG[J,...]
        # remove general islands
        ag=(AG[J,...]!=bkvalue)
        L=label(ag)
        N=len(np.unique(L))-1
        # print(N)
        if N>1:
            RE=regionprops(L)
            MAX=0; 
            for I in np.arange(0,N):
                if RE[I]['area']>MAX:
                    MAX=RE[I]['area']
                    ID=I
            MASK=ag*0; MASK=np.uint8(L==ID+1)
            AGNew[(ag-MASK)==1]=bkvalue
         # remove label based islands
        for K in range(3):
            ag=(AG[J,...]==K)
            L=label(ag)
            N=len(np.unique(L))-1
            if N>1:
                RE=regionprops(L)
                MAX=0; 
                for I in np.arange(0,N):
                    if RE[I]['area']>MAX:
                        MAX=RE[I]['area']
                        ID=I
                MASK=ag*0; MASK=np.uint8(L==ID+1)
                AGNew[(ag-MASK)==1]=bkvalue
                ddd=1
        AG[J,...]=AGNew   
    return AG
model1=keras.models.load_model('Model1.h5') # loading deep learning pre-trained model 1
model2=keras.models.load_model('Model2.h5') # loading deep learning pre-trained model 2
f=sio.loadmat('Raw.mat') # loading cine CMR input data as 4D array with [x,y,n slices, n frames]
A1=f['M']
SLNum=A1.shape[2];
S=[A1.shape[0], A1.shape[1]]
Frames=A1.shape[3]
AAG=np.zeros((SLNum,128,128,Frames))
for J in range(Frames):
    im=np.zeros((S[0],S[1],SLNum))
    try:
        for I in range(SLNum):
            im[:,:,I]=A1[:,:,I,J]
    except:
       continue 
    X,Rat=prep(im)
    X=np.rollaxis(X,2); 
    X=np.reshape(X,[X.shape[0],X.shape[1],X.shape[2],1])
    Y2=model1.predict(X) #prediction stage 1
    X2=np.concatenate([X,Y2],axis=-1)
    Y3=model2.predict(X2) # prediction stage 2
    AG=np.argmax(Y3,axis=3);
    AG=remisol(AG,3)
    AAG[:,:,:,J]=AG
    print('Frame '+str(J))
    AAG=np.uint8(AAG)
sio.savemat('seg.mat', {'AAG':AAG}) # saving the output segmented array as seg.mat readable by MATALB 
V=np.zeros((Frames,SLNum,3))
