from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Data.CodonTable import TranslationError
from exceptions import IndexError
from subprocess import call
import sys, getopt, re, os, shutil, inspect
import glob
from exceptions import IOError

	

#python-scripts
import RATT_Correction
import mergeFeatures

#print to stderr (for piping)
def eprint(string) :
	print >> sys.stderr, string

#main
def main(argv):
	'''
	Parameters

	--rast		<Path> (a embl-file containing a RAST-Annotation for the file specified in query)

	--reference	<Path> (a folder containing one or more embl-files for transfer)

	--query		<Path> (a fasta-file of the organism without annotations)

	--output	<*Path>/<filename> (without extension)

	#other options
	--biopython	Use Algorithm of biopython for the correction of Color-Code 
	--pipe		print final output-file to STOUT (for piping)

	---
	*is not necessary

	--------
	Note:
	- Output-Dir of all files will be in the current dir, for the final file a specific Dir can be specified
	- Since RATT outputs a not correctly formated header (see embl-file specifications), the header of the RAST-file will be used for the final output


	Example
	>python script.py --rast CUL_Rast.embl --reference /embl_files --query CUL.fasta --output $HOME/data/output/CUL_Final

	'''

	#Parameters
	rast_file = ''
	ratt_file = ''
	query_fasta = ''
	reference_folder = ''

	#name of output-file
	outputfile = ''
	#prints output to STOUT
	bool_OutputPipe = False

	#marks features as pseudo based on biopython-lib
	bool_BioPythonAlgo = False

	try:
		opts, args = getopt.getopt(argv,"hpbl:o:c:r:",["rast=","output=","pipe","biopython","query=","reference="])

		for opt, arg in opts:
			if opt == '-h':
				print main.__doc__
				sys.exit()

			elif opt in ("l","--rast"):
				rast_file = arg

			elif opt in ("b","--biopython"):
				bool_BioPythonAlgo = True

			elif opt in ("-p","--pipe"):
				bool_OutputPipe = True

			elif opt in ("-o","--output"):
				outputfile = arg
			
			elif opt in ("-c","--query"):
				query_fasta = arg

			elif opt in ("-r","--reference"):
				reference_folder = arg

		if (query_fasta == '' or rast_file == '' or reference_folder == '' ) : 
			raise getopt.GetoptError ("Provide query, rast-files and reference Folder")

	except ( getopt.GetoptError, IndexError ) as e:
		print 'Wrong CLI-Arguments, see script.py -h', str(e)
		sys.exit(2)
	
	#Check Installations
	try:
		os.environ['RATT_HOME']
	except KeyError:
		print "The ENV $RATT_HOME is not set"
		sys.exit(1)

	##Program:::::::::::::::::::::

	#-----Call Ratt
	
	#Test running ratt

	#Building the Ratt-Call
	tag = 'RATT_TAG'
	binary = os.environ['RATT_HOME']+'/start.ratt.sh'

	#call Ratt
	call([binary, reference_folder, query_fasta, tag , 'Strain'])

	#-----Obtaining Ratt-Ouput
	final_file_regex = tag+".*.final.embl"
	files = glob.glob(final_file_regex)
	
	#Check Output
	if not len(files) == 1 :
		 raise IOError ("Either the RAST-Call went wrong or you were using a fasta-file with several replicons")
	
	ratt_file = files[0]

	print 'This is the output of ratt : '+ratt_file


	#Correction1: Correcting embl-file to ensure biopython-compatibility: taking the header of rast-file

	ratt_file_corr = str(os.path.splitext(files[0])[0])+'_corr1.embl' #path
	print ratt_file_corr
	file_ratt_corr = open(ratt_file_corr , 'w')	#filehandle

	#Transfer of the header (rast->ratt) and writing to file_ratt_corr	
	with open(ratt_file) as file_ratt:
		header_ratt = True
		line_ratt = ''
		while header_ratt:
			line_ratt = file_ratt.readline()
			if line_ratt.startswith('FT'):
				header_ratt = False

		with open (rast_file) as file_rast:
			header_rast = True
			line_rast = ''
			while header_rast : 
				line_rast = file_rast.readline()
				if line_rast.startswith('FT'):
					header_rast = False
					break

				file_ratt_corr.write(line_rast)
		
		#Copy rest of the file
		file_ratt_corr.write(line_ratt)
		shutil.copyfileobj(file_ratt, file_ratt_corr)
		#flush!
		file_ratt_corr.close()
		file_ratt.close()
		file_rast.close()
	
	#switching namespace
	ratt_file = ratt_file_corr

	#Testing RATT-ouput for biopython
	try :
		ratt_record = SeqIO.read(open(ratt_file,"r"), "embl")
		print ratt_record
	except ValueError as e:
		print "The output of RATT doesn't seem to be compatible with [biopython] : "+str(e)

	#Correction2: Calling RATT_Correction.py --> See file for documentation

	ratt_file_corr = str(os.path.splitext(files[0])[0])+'_corr2.embl' #path


	if (bool_BioPythonAlgo):
		RATT_Correction.main(['--ratt',ratt_file,'--output',str(os.path.splitext(ratt_file_corr)[0]),'--biopython'])
	else :
		RATT_Correction.main(['--ratt',ratt_file,'--output',str(os.path.splitext(ratt_file_corr)[0])])

	#switching namespace
	ratt_file = ratt_file_corr

	#Merging Features of RATT with the Features of RATT

	if outputfile == '' :
		outputfile = str(os.path.splitext(files[0])[0])+'_corr3' #path

	mergeFeatures.main(['--rast',rast_file,'--ratt',ratt_file,'--output',outputfile])

if __name__ == "__main__":

	main(sys.argv[1:])


