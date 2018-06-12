"""Misc. useful things"""

class Tracker(object):
	"""Tracks progress of task"""
	def __init__(self, name, size, currPercentage=0, count=0):
		self.name = name	#name of task
		self.size = size	#total size of task (e.g. number of lines in file)
		self.currPercentage = currPercentage	#current percentage of task complete
		self.count = count		#absolute amount of task complete (e.g. number of lines of file read)
	def increment(self):
		if self.size !=0 and self.size is not None:
			self.count += 1
			newPercentage = self.currPercentage + 1
			if float(self.count)/self.size >= float(newPercentage)/100:	#if at least X% of the file has been read, print percentage
				self.currPercentage = newPercentage
				print "{} {}% complete".format(self.name, self.currPercentage)	

def args_are_valid(args, names, intervals):	
	valid_args = True
	for (arg, name, interval) in zip(args, names, intervals):
		lower_bound = interval[0]
		upper_bound = interval[1]
		if lower_bound is not None:
			if arg <= float(lower_bound):
				print "Error. {} must be > {}.".format(name, lower_bound)
				valid_args = False
		if upper_bound is not None:
			if arg >= float(upper_bound):
				print "Error. {} must be < {}.".format(name, upper_bound)
				valid_args = False
	return valid_args

def get_res_string(res):
	"""Converts resolution in bp to string (e.g. 10kb)"""
	res_kb = res/1000
	if res_kb < 1000:
		return str(res_kb) + "kb"
	else:
		return str(res_kb/1000) + "mb"	 
