
"""
track.py

Use for object tracking in video frame...specifically for tracking the black headpiece of the mouse during free movement.
This will convert the video/image frame to Mono8 and threshold the BW figure and invert, as it assumes the headpiece is black.
If the mouse is also black, will likely track the entire headpiece-mouse complex, if not, should track the headpiece nicely.

Note this will not work well if the background is also black!

Assumes that numpy/scipy and openCV (as cv2) are already installed on your computer
"""

# imports
import numpy as np
import scipy as sp
import argparse
import cv2
import csv
import os
import time

# initialize the current frame of the video, along with the list of
# ROI points along with whether or not this is input mode
frame = None
roiPts = []
center = []
inputMode = False

print
#outdir = 'Z:\Jordan\Yackle\Videos'
#outfile = 'test'
outdir = raw_input("output directory: ")
outfile = raw_input("ouput file name: ")
if not os.path.isdir(outdir):
	os.mkdir(outdir)
os.chdir(outdir)

#=============================
def selectROI(event, x, y, flags, param):
    # grab the reference to the current frame, list of ROI
    # points and whether or not it is ROI selection mode
    global frame, roiPts, inputMode

    # if we are in ROI selection mode, the mouse was clicked,
    # and we do not already have four points, then update the
    # list of ROI points with the (x, y) location of the click
    # and draw the circle
    if inputMode and event == cv2.EVENT_LBUTTONDOWN and len(roiPts) < 4:
        roiPts.append((x, y))
        cv2.circle(frame, (x, y), 4, (0, 255, 0), 2)
        cv2.imshow("frame", frame)

#=============================
def preProcess(img):
    # globals
    kernel = np.ones((2,2),np.uint8) # for denoising binary thresholded images
	threshold = 70

    # convert to gray scale
    img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # blur and threshold
    img = cv2.blur(img, (3,3))
    _,img = cv2.threshold(img,threshold,255,cv2.THRESH_BINARY_INV) # assuming headpiece is black

    # smooth the thresholded to remove noise
    img = cv2.morphologyEx(img, cv2.MORPH_CLOSE, kernel, iterations=1) # remove holes in the image
    img = cv2.morphologyEx(img, cv2.MORPH_OPEN, kernel, iterations=2) # remove speckled noise

    return img

#=============================
def boundingBox():
    # global variables
    global frame, roiPts, inputMode, threshold, kernel

    # indicate that we are in input mode and clone the
    # frame
    inputMode = True
    orig = frame.copy()
    roiPts = [] # overwrite previous roi points

    # keep looping until 4 reference ROI points have
    # been selected; press any key to exit ROI selction
    # mode once 4 points have been selected
    while len(roiPts) < 4:
        cv2.imshow("frame", frame)
        cv2.waitKey(0)

    # determine the top-left and bottom-right points
    roiPts = np.array(roiPts)
    s = roiPts.sum(axis = 1)
    tl = roiPts[np.argmin(s)]
    br = roiPts[np.argmax(s)]

    # grab the ROI for the bounding box and convert it
    # to the GRAY colorspace, then perform an inverse binary threshold
    # and remove speckled noise
    #orig = preProcess(orig)

    # compute a histogram for the ROI and store the
    # bounding box as well as the center of the box
    roi = orig[tl[1]:br[1], tl[0]:br[0]]
    roiHist = cv2.calcHist([roi], [0], None, [16], [0, 180])
    roiHist = cv2.normalize(roiHist, roiHist, 0, 255, cv2.NORM_MINMAX)
    roiBox = (tl[0], tl[1], br[0], br[1])

    return roiBox, roiHist

#=============================
def main():
    # construct the argument parse and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("-v", "--video",
        help = "path to the (optional) video file")
    args = vars(ap.parse_args())

    # grab the reference to the current frame, list of ROI
    # points and whether or not it is ROI selection mode
    global frame, roiPts, center, inputMode, threshold, kernel, outdir, outfile

    # if the video path was not supplied, grab the reference to the
    # camera...else load the reference video
    if not args.get("video", False):
        camera = cv2.VideoCapture(0) # change this to readimage!
    else:
        camera = cv2.VideoCapture(args["video"])

    # setup the mouse callback
    cv2.namedWindow("frame")
    cv2.setMouseCallback("frame", selectROI)

    # initialize the termination criteria for cam shift, indicating
    # a maximum of ten iterations or movement by a least one pixel
    # along with the bounding box of the ROI
    termination = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 10, 1)
    roiBox = None
    grabbed = True
    # ==========================================================================
    # Perform the action over each image frame
    # keep looping over the frames
    timestamp = []
	framecount = []
	count = 0
    while grabbed == True:
        # grab the current frame
        grabbed, frame = camera.read()

        # check to see if we have reached the end
        if not grabbed:
            break

        # if the ROI has been computed
        if roiBox is not None:
			count += 1
            framecount.append(count)
			timestamp.append(time.time())

            # preprocess the frame and calculate the back projection
            #frame2 = preProcess(frame)
            backProj = cv2.calcBackProject([frame], [0], roiHist, [0, 100], 1)

            # apply cam shift to the back projection, convert the
            # points to a bounding box, and then draw them
            r, roiBox = cv2.CamShift(backProj, roiBox, termination)
            pts = np.int0(cv2.cv.BoxPoints(r))
            cv2.polylines(frame, [pts], True, (0, 255, 0), 2)

            # append the center of the object position to "centers"
            temp = np.int32(r[0])
            c = (temp[0],temp[1]) # converted to integers
            cv2.circle(frame, c, 1, (0, 0, 255), thickness=2) # draw the center circle
            center.append(c) # store the center information
        else:
            print('Please draw the bounding box roi points')
            roiBox,roiHist = boundingBox() # call function to draw ROI

        # show the frame and record if the user presses a key
        cv2.imshow("frame", frame)
        key = cv2.waitKey(1) & 0xFF

        # handle if the 'i' key is pressed, then go into ROI
        # selection mode...will allow us to redraw the ROI points
        if key == ord("i"):
            roiBox, roiHist = boundingBox()

        # if the 'q' key is pressed, stop the loop
        elif key == ord("q"):
            break

    # cleanup the camera and close any open windows
    camera.release()
    cv2.destroyAllWindows()

    # save the center positions as a csv file and get cumulative distance travelled
    print("saving the data...")
    os.chdir(outdir)
    filename = outfile + ".csv"
    np.savetxt(filename, np.c_[timestamp,framecount,center], fmt='%10.2f %06d %i %i', delimiter = ',')

if __name__ == "__main__":
    main()
