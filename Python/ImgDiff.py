import cv2
import numpy as np
import matplotlib.pyplot as plt

file = "Z:\\Jordan\\Yackle\\Videos\\3813168_20160509-1.avi"
cam = cv2.VideoCapture(file)
points = []
count = 0
grabbed = True

#==============================================
def process(frame):
    kernel = np.ones((2,2),np.uint8)
    smoothing = np.ones((5,5),np.uint8)

    # series of erosions and dilations to eliminate background
    frame = cv2.morphologyEx(frame,cv2.MORPH_ERODE,smoothing,iterations=1)
    frame = cv2.morphologyEx(frame,cv2.MORPH_DILATE,smoothing,iterations=2)
    frame = cv2.morphologyEx(frame,cv2.MORPH_OPEN,smoothing,iterations=2)
    frame = cv2.morphologyEx(frame,cv2.MORPH_CLOSE,smoothing,iterations=3)

    return frame
#=================================================

while grabbed == True:
    grabbed,frame = cam.read()
    frame = cv2.blur(frame,(2,2)) # smooth the image
    thresh = cv2.cvtColor(frame,cv2.COLOR_BGR2GRAY) # grayscale
    _,thresh = cv2.threshold(thresh,60,255,cv2.THRESH_BINARY_INV) # binary threshold

    thresh = process(thresh) # create binary threshold

    if count > 200:
        dframe = thresh - preframe
        dframe = process(dframe)
        M = cv2.moments(dframe)
        #x = np.uint(['m01']/M['m00'])
        #y = np.uint(M['m10']/M['m00'])

        #cv2.circle(dframe,(x,y),2,(0, 0, 255),thickness=2)
        #cv2.imshow("dframe",dframe)
        plt.subplot(2,2,1)
        plt.imshow(dframe,cmap='gray')
        plt.subplot(2,2,2)
        plt.imshow(frame,cmap='gray')
        plt.subplot(2,2,3)
        plt.imshow(thresh,cmap='gray')
        plt.subplot(2,2,4)
        plt.imshow(preframe,cmap='gray')
        plt.show()

        #points.append((x,y))
    else:
        preframe = thresh

    count += 1
