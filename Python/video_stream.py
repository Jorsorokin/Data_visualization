"""
video_stream.py

Uses openCV to acquire video frames from a connected webcam/CV-compatible camera.
The camera images are streamed as .jpgs to the folder specified by the user.
The software will create a new directory called "video" where it will save the images
as well as a .csv file containing timestamps of each image and meta-data including
the date-time created, the frame number, the average frame rate.

* Note: if you type "1" when the software asks "wait for subfolder?", it will
continuously monitor the root directory and only stream if it detects a new folder.
This is a hack in order to pseudo-synchronize the camera with hardware like the
TDT or OpenEphys while recording data such as photometry or ephys *

Dependencies:
- openCV (as cv2)
- openCV-compatible camera
- numpy
- concurrent

By Jordan Sorokin, 3/28/2017
"""

# imports
import cv2, numpy, csv, time, os
from concurrent import futures # to speed up saving
executor = futures.ThreadPoolExecutor(max_workers = 3)

# Functions
#-------------------------
def waitToCapture(directory):
	# wait for a subfolder creation in the rootfolder "directory"
	prechange = dict( [(files,None) for files in os.listdir( directory )] )
	while True:
		postchange = dict( [(files,None) for files in os.listdir( directory )] )
		change = [files for files in postchange if not files in prechange]
		if change:
			break
	print change
	return change

#-------------------------
def saveImg(img,counter):
    name = "{0:0=5d}".format(counter)
    cv2.imwrite(name + ".jpg", img)

#------------------------
def changeDirectory(directory):
	if not os.path.isdir(directory):
		os.mkdir(directory)
	os.chdir(directory)

#-------------------------
def main(outdir,executor,waiting=0):

	# presets/looping params
	grabbed = True # will break loop if camera fails for some reason
	count = 0
	timestamp = []
	framecount = []

	# initialize the camera
	camera = cv2.VideoCapture(0) # if only 1 camera plugged in, will be the first index (0)
	print('Camera initialized!')

	# optimize function calls
	readFrame = camera.read
	showFrame = cv2.imshow
	saveFrame = executor.submit
	closeFrame = cv2.waitKey

	# check for changes in root directory
	if waiting == 1:
		print('Waiting for action')
		newdir = waitToCapture(outdir)
		outdir = outdir + '\\' + newdir[0] # concatenate the newly detected subfolder
		timestamp.append( time.time() ) # append the time of the folder initialization...for synching the camera
		framecount.append(count) # frame 0 will be the start of the folder...no frames collected yet

	# concatenate a "video" subfolder for storage
	outdir = outdir + '\\' + 'video' # make a new subfolder to store the images and .csv file
	changeDirectory(outdir)

	# begin infinite loop and write the image to the specified folder
	print
	print('Press "q" to end video capture')
	while grabbed == True:

		# read an image & update the timestamp/framecount
		grabbed, frame = readFrame()
		count += 1
		timestamp.append( time.time() )
		framecount.append(count)

		# show the frame
		showFrame('Video',frame)

		# save the image
		saveFrame(saveImg,frame,count)

		# check for user input
		key = closeFrame(1) & 0xFF
		if key == ord('q'):
			break

		"""
		Still need to figure out how best to automatically stop acquisition if waiting for
		hardware. Maybe checking the file-sizes of the folders, and seeing if any file size
		stops increasing compared to last loop. However, this assumes that the file sizes
		change continuously, whereas they may change via data-block dumping at some predefined
		time interval. Also assumes that file sizes would be changing anyway...

		Could try to access more of the low-level features of the system to see if
		windows has changed its representation of the initialized folder. Of course, if no
		changes are expexcted to occur beyond the initial change when delaying camera acquisition,
		then one can manually stop acquisition by pressing "q".
			- JMS
		"""

	# close the camera
	camera.release()
	cv2.destroyAllWindows()

	# save the timestamp .CSV file
	print
	print("Frames displayed: %i" % count)
	print("Average framerate: " + str(count / (timestamp[-1] - timestamp[0])))
	print
	print("saving timestamps...")
	os.chdir(outdir)
	filename = 'timestamps.csv'
	numpy.savetxt(filename, numpy.c_[timestamp,framecount], fmt='%10.2f %i', delimiter = ',')


# prompt the user for the root directory to save files, and whether or not to
# wait for any subfolders
outdir = raw_input("output directory: ")
waitForSubfolder = input("wait for subfolder? [1/0] ")
changeDirectory(outdir)

# call "main" function
main(outdir,executor,waiting=waitForSubfolder)
