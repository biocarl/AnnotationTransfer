from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Data.CodonTable import TranslationError
import sys
import re
from exceptions import IndexError
import sys, getopt


#print to stderr (for piping)
def eprint(string) :
	print >> sys.stderr, string

#main
def main(argv):

	################################################
	#This is a  comparsion function to determine intersection of two feature positions 
	def compare ( feature1 , feature2 ) : 

		if isContained (feature1, feature2) :
			return 0

		if isLeft (feature1, feature2) :
			return -1

		if isRight (feature1, feature2) :
			return 1
	#bool if there is an interference of the features locations
	def isContained( feature1, feature2 ):

		bool = False

		if (feature1.location.start in feature2 or feature1.location.end in feature2 or feature2.location.start in feature1 or feature2.location.end in feature1):
			bool = True
		
		#This is for the case that start and end of the features locations are lying exactly in the gaps of the Compound-Location
		if isinstance(feature1.location, CompoundLocation) or isinstance(feature2.location, CompoundLocation) :
			return not isLeft(feature1, feature2) and not isRight(feature1, feature2)

		return bool

	# feature2 is left of feature1 (no interference)
	def isLeft (feature1, feature2 ) : 
		return ( feature2.location.end <= feature1.location.start )

	# feature2 is right of feature1 (no interference)
	def isRight ( feature1, feature2) :
		return ( feature2.location.start >= feature1.location.end )



	#This is classic binarySearch, which also holds a global variable with the last used interval if the search has no match
	def binSearch (list, element ):
		def binarySearch_recursive(list, element , start, end) :
			global feature_ratt_left,feature_ratt_right

			if end < start:
				return False
			mid = (start + end)/2

			if isContained(list[mid],element):
				return True
			elif isRight(list[mid],element):
				feature_ratt_left = mid +1
				feature_ratt_right = end
				return binarySearch_recursive(list, element, mid+1, end)
			else:
				feature_ratt_left = start
				feature_ratt_right = mid -1
				return binarySearch_recursive(list, element, start, mid-1)

		global feature_ratt_left,feature_ratt_right
		return binarySearch_recursive( list, element, 0, len (list)-1)


	#This is a simple comparison function to sort features according to their position
	def compare_featureSort(feature1, feature2):
		if feature1.location.start < feature2.location.start:
			return -1
		elif feature1.location.start > feature2.location.start:
			return 1
		else:
			if feature1.location.end < feature2.location.end :
				return -1
			elif feature1.location.end > feature2.location.end :
				return 1
			else : 
				return 0


	#Parameters
	ratt_file = ''
	rast_file = ''
	#name of output-file
	outputfile = 'output'
	#prints output to STOUT
	bool_OutputPipe = False

	#last positions before end of BinarySearch
	global feature_ratt_left,feature_ratt_right
	feature_ratt_left = 0
	feature_ratt_right = 0

	#

	try:
		opts, args = getopt.getopt(argv,":l:fph:o",["rast=","ratt=","output=","pipe"])

		for opt, arg in opts:
			if opt == '-h':
				eprint ("For options see main.py")
				sys.exit()
	
			elif opt in ("l","--rast"):
				rast_file = arg

			elif opt in ("f","--ratt"):
				ratt_file = arg

			elif opt in ("-p","--pipe"):
				bool_OutputPipe = True

			elif opt in ("-o","--output"):
				outputfile = arg

		if (ratt_file == '' or rast_file == '') : 
			raise getopt.GetoptError ("Provide Files")

	except ( getopt.GetoptError, IndexError ) :
		eprint ('Wrong CLI-Arguments, see script.py -h')
		sys.exit(2)
	
	#CODE
	###################################


	ratt_record = SeqIO.read(open(ratt_file,"r"), "embl")

	#here the unique features of RAST are stored
	uniq_features = list()

	for rast_record in SeqIO.parse(open(rast_file,"r"), "embl") :
		count = 0
#		#Iteration RAST
		for feature_rast in rast_record.features :	

		#Feature locations of RAST which are not interfering with feature locations transferred by RATT
			if (not binSearch(ratt_record.features , feature_rast)) :

				#only CDS-type is transfered  of RAST-featrues, if you want all (e.g. tRNA, gene,..) then delete condition
				#if feature_rast.type == 'CDS' :
				uniq_features.append(feature_rast)
				#here you can tag the new rast-features
				#feature_rast.qualifiers.update({'rastAnnotated':''})

					
	#Inserting the new features on the right spot, by appending and then sorting the whole list (the way above is faster but doesn't work so far)
	ratt_record.features.extend(uniq_features)
	ratt_record.features.sort(compare_featureSort)

	################
	#Save everything to a new File

	if ( bool_OutputPipe ) :
		#Print to STOUT
		print(ratt_record.format("embl"))
	else : 
		output_handle = open(outputfile+".embl", "w")
		SeqIO.write(ratt_record, output_handle, "embl")
		output_handle.close()

if __name__ == "__main__":

	main(sys.argv[1:])

