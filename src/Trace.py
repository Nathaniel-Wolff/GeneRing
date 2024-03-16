import scipy
#NW Edits: changed import of numpy to np
import numpy as np
import scipy.spatial
import itertools
import Chromatin

__version__="01.00.01"
__original_author__ ="Robert Shelansky"
__author__ = "Nathaniel Wolff"


class Dummy:
	"""Class is called if molecule's midpoint cannot be approximated. So, in the call in moleculify"""
	def __init__(self):
		pass
	def getBubbles(self):
		return np.array([])

	def getBubblesInWindow(self, windowStart = 0, windowEnd = -1):
		return np.array([])

class Trace:
	"""
	Trace represents an individual trace of a molecule from an image file. It can be instantiated by giving it
	a reference to a list of tuples of coordinates from starting at the end of the molecule (the fork) and leading
	to the beginning of the molecule. It contains all the methods required to do analysis of an individual trace.
	It also contains helper functions which do both the analysis and manipulate an individual "realization of the fit" --
	meaning, that for each coordinate/point, the function determines whether it's part of a linker or not.
	"""
	def __init__(self, trace):
		"""
		Takes a list of coordinate tuples and computes metrics required for realizing a specific bubble linker path.
		usable metrics are as follows.
		_trace:
			#array of x,y coordinates of one single _trace
		_ld:
			#distance between successive points linked diff (ld) and the distance between all points as a matrix (d)
			#index 0 refers to distance between 0,1 in _trace
			#index -1 refers to distance between -2,-1 in _trace
		_cld:
			#cumulative distance between coordinates starting at 0,1
			#there is no index 0 
			#index i refers to the distance traveled to get to index i+1 in _trace	
		_ll:
			#length of the whole molecule in the coordinate system
		_d:
			#distance between every point and every other
		"""
		#print("trace", trace[:10])
		self._trace      =scipy.array(trace)

		self._ld         =scipy.array([scipy.spatial.distance.euclidean(i,j) for i,j in zip(self._trace[:-1],self._trace[1:])])


		self._cld        =scipy.concatenate(([0],scipy.cumsum(self._ld)))

		self._ll         =scipy.sum(self._ld)
		#
		self._dMatrix          =scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(self._trace,'euclidean'))
		#NW Edit: the below originally called on scipy to get the .ma attribute but this attribute is actually in numpy (imported above)
		#NW Edit: there were two _._d attributes and the first one was bound to a different name for clarity

		#this, the actual d, masks self distances
		self._d          =np.ma.masked_where(self._dMatrix==0,self._dMatrix) #mask self distances

		#compute the length and save to the object
		self.__len__()

	def __getitem__(self, arg):
		"""
		returns the coordinate at index arg. (   [x,y]  )
		"""
		return (self._trace[arg])

	def __len__(self):
		"""
		returns the length in coordinates of the trace. The number of coordinate entries.
		"""
		#NW Edits: Made __len__ a class attribute; see Moleculify for explanation
		self.length = len(self._trace)
		return(len(self._trace))

	#below is a block of helper methods: partition, smooth, and segment.
	# Each of these is used in regionify, which receives a guess of the midpoint to start with
	#this is used in the initial call of partition

	def partition(self, midpoint):
		"""
		returns two lists forward and reverse who have length len(self)= len(forward) + len(reverse)
		each index contains the distance of the closest coordinate of the other strand of the molecule as defined by midpoint
		"""

		forward=scipy.amin(self._d[0:midpoint+1,midpoint+1:len(self)],axis=1)
		reverse=scipy.amin(self._d[0:midpoint+1,midpoint+1:len(self)],axis=0)
		return forward, reverse

	def smooth(self, mol, smooth):
		"""
		takes a list where every element is some value and performs a sliding window average around that coordinate
		creating a new list whose length is the same as inputed.
		"""
		mol         =scipy.convolve(mol,scipy.ones(smooth)/smooth,mode='same')
		return(mol)

	def segment(self, mol, threshold, end):
		"""
		takes a list whose values are some score of whether that specific coordinate is a linker.
		this is typically used on the lists returned from partition but does not have to be.
		All stretches of elements whose threshold is above or below some value are then determined.
		and a list of regions and there corresponding label are returned. It also ignores the last "end"
		coordinates of mol.

		regions = [(0,end),(start,end),...,(start,end),(start,len(mol))]
		labels   = [1,0,...1,0]  1 for greater then threshold 0 for below.

		"""
		#NW Debugging:
		#fixed mask
		#mask previously returned an array of Boolean. scipy.ediff1d cannot take non-float 64 or int 64 arguments
		#now mask returns an array of the integer equivalents of booleans as intended
		mask        = np.array(mol > threshold, dtype = np.int64)

		#NW Debugging:
		#fixed breaks: Changed scipy.concatenate to numpy.concatenate
		breaks      =np.concatenate(( [0],scipy.argwhere(scipy.ediff1d(mask[:-end]))[:,0]+1, [len(mol)] ))

		regions     =scipy.array([(s,e) for s,e in zip(breaks[:-1],breaks[1:])])
		labels      =scipy.array([int(mask[s]) for s,e in regions])
		return(regions,labels)
	def regionify(self, guess, threshold=4, smooth=10, end=5):
		"""
		returns a pair of tuples that is a best guess at regions given a specific
		midpoint, threshold, end, and smooth factor. A convenience function which calls other
		trace methods.
		"""

		#NW Debugging:
		#regionify calls the partition method to generate a,b.


		#using midpoint guess, forms forward and reverse strand lists
		a,b  = self.partition(guess)
		#NW Edits: changed regionify tuples to same variable names as in regionify (is this ok to do??)

		#using the a smoothed forward list and determines which regions are linkers or not
		#see documentation for segement

		#ff = forward regions, fl = forward region labels (1 = bubble, 0 = linker)
		ff,fl = self.segment  (self.smooth(a      , smooth), threshold=threshold, end= end)

		# rf = reverse regions, fl = reverse region labels (1 = bubble, 0 = linker)
		rf,rl = self.segment  (self.smooth(b[::-1], smooth), threshold=threshold, end= end)

		rf    = len(self)-rf[:,[1,0]]
		rl   = rl


		return((ff,fl),(rf,rl))

	def _midpoint(self, left, right):
		"""
		given two coordinates _midpoint returns the middle coordinate in terms of actual distance.
		"""
		return(scipy.searchsorted(self._cld, self._cld[left]+(self._cld[right]-self._cld[left]) / 2 ))

	#Not finished reviewing...
	def midpoint(self, guess = 0, threshold=4, smooth=10, end=5):
		#NW Edits: Decided to call this method at initialization such that the midpoint can be a class attribute.
		#This allows other methods that call on the midpoint to easily access it without a call to the midpoint method
		###Going to set guess' default to 0 for now.
		#NW Edit: Moved this method up for readability.
		"""
		Takes some parameters and returns the best guess of the midpoint.
		It does this iteratively. It attempts to define regions based on guess.
		Then, looks at the second to last set of regions in the shorter strand and attempts to recalculate the midpoint
		given that region and its corresponding region in the longer strand. does this until convergence.
		or (until it attempts number of regions found in one long molecule)--upper bound not super important.
		helps in absence of convergence.
		"""

		"""NW Edit: This now saves the midpoint parameters to the class. 
		This makes calls to functions that use these parameters easier to do, just referring to the instance of the object.
		These functions are solve_molecule, label, and scale. """
		#NW Debugging:
		#call to regionify didn't had guess as posititional argument, changed to a keyword argument
		#changed regionify tuples to same variable names as in regionify (is this ok to do??)

		#either user supplies the guess or scipy guesses the midpoint using 1/2 of cumulative distance
		guess           = guess or scipy.searchsorted(self._cld, (self._cld[-1]/2))
		#print("cld", self._cld[-1]/2)

		#regionify call; returns regions and bubble/linker indices
		(ff,fl),(rf,rl) = self.regionify(guess = guess, threshold = threshold, smooth=smooth, end=end)
		for i in range(min(len(ff),len(rf))):
			i         = min(len(ff),len(rf))-2
			new_guess = self._midpoint(ff[i][1],rf[i][0])
			if guess==new_guess:
				#NW Debugging:
				#Before guess was saved as a class attribute after exiting the for loop, but this is after the function call ends
				#Now it is saved before return
				self.guess = guess
				return(guess)

			guess = new_guess
			(ff,fl),(rf,rl) = self.regionify(guess = guess, threshold = threshold, smooth=smooth, end=end)

		#NW Edits: Made midpoint parameters class attributes; see moleculify for explanation
		#To do this, midpoint is called at instantiation and is saved as an attribute, used to define a and b
		#print("done")


		self.ff = ff
		self.fl = fl
		self.rf = rf
		self.rl = rl
	def zip(self, fregions, flabels, rregions, rlabels):
		"""
		takes the definitions of forward regions and reverse regions and returns the concatenated version.
		This is simply the region definitions for the molecule in its entirety
		"""
		regions = scipy.concatenate((fregions ,rregions))
		labels  = scipy.concatenate((flabels,rlabels))
		return(regions,labels)

	def msd_of_linker_end_points(self,fregions,flabels,rregions,rlabels):
		"""
		takes definitions of regions and labels for forward and reverse
		strands assuming each region in the forward direction corresponds to a region
		in the reverse direction and returns the mean squaared distance
		of the start and end points of each region. 
		"""
		if len(fregions) != len(rregions):
			return(float('Inf'))
		flinks=fregions[scipy.logical_not(flabels)]
		rlinks=rregions[scipy.logical_not(rlabels)]
		s=sum([ self._d[t1,t2] for t1,t2 in zip(flinks[:,0]  ,rlinks[:,1] -1)  if self._d[t1,t2]])   /len(flinks)
		f=sum([ self._d[t1,t2] for t1,t2 in zip(flinks[:,1]-1,rlinks[:,0])     if self._d[t1,t2]])   /len(flinks)
		return((s+f)/2)

	def msd_of_region_sizes(self,fregions,flabels,rregions,rlabels):
		"""
		takes definitions of regions and labels for forward and reverse
		strands assuming each region in the forward direction coorosponds to a region
		in the reverse direction and returns the mean squared distance between the sizes of 
		each region.
		"""
		if len(fregions) != len(rregions):
			return(float('Inf'))
		flen = self._cld[fregions[-1,-1]-1] - self._cld[fregions[0,0]]
		rlen = self._cld[rregions[0,-1] -1] - self._cld[rregions[-1,0]]
		dif  = sum([((self._cld[ff-1]-self._cld[fs])/flen-(self._cld[rf-1]-self._cld[rs])/rlen) **2            for (fs,ff),(rs,rf) in  zip(fregions,rregions)])
		dif  = dif*(flen+rlen)/2
		return(dif)

	def sd_segment_size(self,fregions,flabels,rregions,rlabels):
		"""
		takes definitions of regions and labels for forward and reverse
		strands assuming each region in the forward direction corresponds to a region
		in the reverse direction and returns how similar in size each fragment is.
		"""
		if len(fregions) != len(rregions):
			return(float('Inf'))
		flen = self._cld[fregions[-1,-1]-1] - self._cld[fregions[0,0]]
		rlen = self._cld[rregions[0,-1] -1] - self._cld[rregions[-1,0]]
		return((flen - rlen)**2)

	# def find_midpoint(self, guess=None ,threshold=4, sensitivity=10, end=5):
	# 	guess           = guess or scipy.searchsorted(self._cld, (self._cld[-1]/2))
	# 	(ff,fl),(rf,rl) = regionify(self, guess, threshold= threshold, sensitivity=sensitivity, end=end)
	# 	i=min(len(ff),len(rf))-2
	# 	guess = _midpoint(self, ff[i][1],rf[i][0])
	# 	return(guess)

	def edgebuffer(self, threshold, smooth):
		"""
		Calculates how many coordinates to ignore on the end by determining
		the ceiling of the minimum number of coordinates to meet threshold
		"""
		return(int(scipy.ceil(threshold/min(self._ld))))

	def solve_molecule(self, midpoint = 0, threshold=4, smooth=10, end=5):
		"""
		given specific settings produces a list of objects which represent a realization of a trace.
		It is again a convenience function like regionify
		"""
		#NW Edits: default for midpoint is "currently" set to 0. This will be adjusted later
		molecule       =self.smooth    (scipy.concatenate(self.partition(midpoint)),smooth=smooth)
		(fr,fl),(rr,rl)=self.regionify (midpoint, threshold=threshold , smooth=smooth,end=end)
		regions, labels = self.zip       (fr,fl,rr[::-1],rl[::-1])
		return (midpoint,molecule,(fr,fl),(rr,rl),(regions,labels))



	def returnLocalBubblesOrLinkers(self, window = (0,0), threshold=4, smooth=10, end=5):
		#omitting midpoint
		"""
		1) Solves the molecule and renames the forward and reverse arrays, bubble/linker arrays
		2) Uses the window supplied by the user to take only part of forward and reverse arrays that belongs to window
		3) Save all of the bubble sizes to an array and return them
		"""

		bubbleSizes = []
		# either user supplies the guess or scipy guesses the midpoint using 1/2 of cumulative distance
		midpointGuess = scipy.searchsorted(self._cld, (self._cld[-1]/2))

		solvedMolecule = self.smooth    (scipy.concatenate(self.partition(midpointGuess)),smooth=smooth)
		(forwardList,forwardStatus),(reverseList,reverseStatus)=self.regionify (midpointGuess, threshold=threshold , smooth=smooth,end=end)

		#scale the window
		scaledWindowTop = round(window[0]/(self.scale(self.length)))
		scaledWindowBottom = round(window[1]/(self.scale(self.length)))

		#searching the forward list for the index of the first region that belongs to the window
		#equivalent to bisect_left

		firstRegionIndex = scipy.searchsorted(
			np.array([region[0] for region in forwardList]), scaledWindowTop, side = "left")

		#same for last region
		lastRegionIndex = scipy.searchsorted(np.array([region[1] for region in forwardList]),
											 scaledWindowBottom, side = "left")

		#iterating through regions spanning the window and collecting the bubble sizes
		for region,identity in zip(forwardList[firstRegionIndex:lastRegionIndex+1],forwardStatus[firstRegionIndex:lastRegionIndex+1]):
			if identity == 0:
				bubble = region
				#rescaling bubble sizes in terms of base pairs and not coordinates
				bubbleSizes.append(   (region[1] - region[0])*  (self.scale(self.length)   ))
		return bubbleSizes

	def label(self,fl,rl):
		"""
		given two lists of regions or labels: it returns a list of of length len(trace) whos
		values are which region it belongs to if you attempt to zip of the molecule from the end.
		"""
		return(scipy.concatenate((list(range(len(fl))),list(reversed(range(len(rl)))))))

	def scale(self,length):
		"""
		calculates the number of basepairs per coordinate distance.
		"""
		return(length / (self._ll / 2))

	def moleculify(self,fr,fl,rr,rl,length):
		"""
		takes a representation of a trace for region definitions and labels and a length in base pairs
		of a molecule and returns a Chromatin.molecule version.
		"""

		#mean length is now a class attribute
		#NW Edits: changed to remove last condition FOR NOW
		#instead of returning none, there should be some kind of dummy instance of the molecule object. how to construct?

		#NW Edit: Counter that returns number of failed traces is now included in moleculify
		failedTraces = 0
		if len(fl)!=len(rl)> 0:
			failedTraces+=1
			thisDummy = Dummy()
			return(thisDummy)
		else:
			region_lengths       = scipy.array([sum((self._cld[r1[1]-2] - self._cld[r1[0]], self._cld[r2[1]-2] - self._cld[r2[0]]))/2  for r1,r2 in zip(fr,rr)])
			exclusive_end_pts    = scipy.ceil(self.length * scipy.ndarray.round(scipy.cumsum(region_lengths)/sum(region_lengths),decimals=3))

			inclusive_start_pts  = scipy.concatenate(([0],exclusive_end_pts[:-1]))

			regions  = scipy.array([(s,e) for s,e in zip(inclusive_start_pts,exclusive_end_pts)])
			molecule=Chromatin.Molecule([ Chromatin.Region(l,self.length-e,self.length-s,e-s) for (s,e),l in reversed(list(zip(regions,self.fl)))])
			return(molecule)

