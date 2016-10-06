from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Data.CodonTable import TranslationError
from exceptions import IndexError
import sys, getopt
import re

########
#TODOS
#1) Pipe doesn't work if local ID is to long, typically done by RAST, shorten Locus-Tag? not feasable for NCBI-Sub anyway : only warning...
#2)


#print to stderr (for piping)
def eprint(string) :
	print >> sys.stderr, string

#main
def main(argv):

	#Parameters
	ratt_file = ''
	#name of output-file
	outputfile = 'output'
	#prints output to STOUT
	bool_OutputPipe = False

	#marks features as pseudo based on biopython-lib
	bool_BioPythonAlgo = False

	try:
		opts, args = getopt.getopt(argv,":l:fpbh:o",["rast=","ratt=","output=","pipe","biopython"])

		for opt, arg in opts:
			if opt == '-h':
				eprint ("For options see main.py")
				sys.exit()
			
			elif opt in ("b","--biopython"):
				bool_BioPythonAlgo = True

			elif opt in ("f","--ratt"):
				ratt_file = arg

			elif opt in ("-p","--pipe"):
				bool_OutputPipe = True

			elif opt in ("-o","--output"):
				outputfile = arg

		if (ratt_file == '') : 
			raise getopt.GetoptError ("Provide Files")

	except ( getopt.GetoptError, IndexError ) :
		print 'Wrong CLI-Arguments, see script.py -h'
		sys.exit(2)

	#CODE
	#######Functions
	#For checking if two features interfere in location (binary search with intervals with open upper limit

	###################
	#Completing Annotations with RAST-Annotations without interfering exsiting annotations
	#Algorithm for searching if a position is included in a range

	####################
	#correct RATT-Transfer

	#delete all features which are not CDS's
	feature_list_CDS = list()

	for ratt_record in SeqIO.parse(open(ratt_file,"r"), "embl") :

		for feature in ratt_record.features :
			
			#If true, feature qualifiers have to be changed
			bool_pseudo=False
			correction_message = "Location " + str(feature.location) + " was marked as /pseudo because of: "
			
			# biopython-Algo for marking genes as pseudo
			if ( bool_BioPythonAlgo ) :
				try :
					translation_seq = str(feature.extract(ratt_record.seq).translate(cds=True))
					#missuse of exceptions :(
				except TranslationError as e:
					correction_message = correction_message + "[Biopython:] " + str(e)
					if feature.type == 'CDS' :
						bool_pseudo=True
						
			#own criteria for marking as pseudo	
			else : 
				#1.Check if feature is a 'join'ed location
				if isinstance(feature.location, CompoundLocation) :
					bool_pseudo=True
					correction_message = correction_message + "|Joined feature|"
				#2.Check for Stop-Codons in feature-location (which don't mark the end of location)
				translation_seq = str(feature.extract(ratt_record.seq).translate())
				searchObj = re.search( r'^.*\*.+$', translation_seq , re.M|re.I)
				if searchObj :
					correction_message = correction_message + "|Interfering Stop-Codon|"
					bool_pseudo=True

			#Correct feature qualifiers
			if bool_pseudo :
				eprint (correction_message)
				#correct colour-code (to red)
				feature.qualifiers.update({'colour':'2'})
				#add qualifier '/pseudo'
				feature.qualifiers.update({'pseudo':''})
			else :
				feature.qualifiers.update({'colour':'3'})

			#save all features of type CDS's
			if feature.type == 'CDS' :
				feature_list_CDS.append(feature)



	#Substitute feature list
	ratt_record.features = feature_list_CDS

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

