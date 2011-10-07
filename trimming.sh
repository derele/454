#!/bin/bash

#EH Mailed form Stephen on 30 Mar 2010
#EH Changes to use Crossmatch instead of Blast
#EH Changes to use seqclean and make the perlies unnecessary

# This is from Sujai's email: 3 Mar 2010 
# Mint trimming: (sfffiles already have trim points set from 5 onwards
# to get rid of initial tcag, so used -n in this case to get the right
# trim points, otherwise would have to add 4 to each starting trim
# point)

# To extract the raw fasta
# sffinfo -s -n file.sff >file_raw.fa

# To trim sff:
#  sfffile -t file_trimpoints.txt -i file_trimpoints.txt -o file_trim.sff file.sff

# (To combine all file_trim.sff into all_trim.sff (-xlr option needed to convert everything to same flowgram space):
# sfffile -xlr -o all_trim.sff *trim.sff)

AdaptersFile="/home/ele/sequences/MINT_SMART_All_Versions.fna"

FastaFileRaw=temporary_file_raw.fna
QualFileRaw=$FastaFileRaw.qual # Crossmatch can use Quality Information
Minscore=20 # Option for Crossmatch: Minimum alignment score
Minmatch=10 # Option for Crossmatch: Minimum alignment length 
Masklevel=100 # Option for Crossmatch: report any match whose domain is not completely contained within a higher scoring match
TrimPointsFile=temporary_trimpoints.txt

sffdir=/home/ele/Data/454/named/
# sffdir=/home/ele/Data/454/Li68_assembly/sff/

# sfffiles[1]=Li68.sff

sfffiles[1]=10F.sff
sfffiles[2]=179F.sff
sfffiles[3]=L2R3.sff
sfffiles[4]=M175.sff
sfffiles[5]=KS4F.sff
sfffiles[6]=UW07F.sff


function exitOnError
  {
  local EXITCODE="$?"
  if [ $EXITCODE -ne "0" ]; then
    echo -e "\nERROR: '$1' failed with error code: '$EXITCODE'\n"
    exit 1
  fi
  }

for InfileSFF in ${sfffiles[@]}
do
    
    InfileSFF=$sffdir$InfileSFF
    Trimfix=$(basename $InfileSFF)
    TrimmedFileSFF=$Trimfix.trimmed.sff 

    echo -e "\nTrimming sff file '$InfileSFF', using adapter file '$AdaptersFile' and crossmatch minscore of $Minscore', and will output to '$TrimmedFileSFF'\n"

# Need to turn sffinfo's lower-case letters to Xes, so that sequclean knows they should be trimmed
    echo "  sffinfo -s -n $InfileSFF  | tr '[:lower:]' 'X' > $FastaFileRaw"
    /home/ele/tools/Newbler2.6/bin/sffinfo -s -n $InfileSFF  | tr '[:lower:]' 'X' > $FastaFileRaw;   exitOnError "sffinfo"

    echo "  sffinfo -q -n $InfileSFF > $QualFileRaw"
    /home/ele/tools/Newbler2.6/bin/sffinfo -q -n $InfileSFF > $QualFileRaw;   exitOnError "sffinfo"
    
    /home/ele/tools/phrap/cross_match.manyreads $FastaFileRaw $QualFileRaw $AdaptersFile -screen -minscore $Minscore -minmatch $Minmatch -masklevel $Masklevel &>/dev/null;   exitOnError "crossmatch" 
    
# Run seqclean to screen for poly-A and dust using 4 processors
    seqclean $FastaFileRaw.screen -c 4 2>&1
    
# With seqclean things get even better, we already have a trimmpoints file
# Here is what the seqclean-readme at http://compbio.dfci.harvard.edu/tgi/software/seqclean_README 
# says:

# Cleaning report format
# ----------------------
# Each line in the cleaning report file (*.cln) has 7 tab-delimited fields 
# as follows:

# 1. the name of the input sequence
# 2. the percentage of undetermined names in the clear range 
# 3. 5' coordinate after cleaning
# 4. 3' coordinate after cleaning
# 5. initial length of the sequence
# 6. trash code
# 7. trimming comments (contaminant names, reasons for trimming/trashing)

# The trash code field (6) should be empty if (part of) a sequence is 
# considered valid - so it can be found in the final filtered file (*.clean)

# The trash code field will be set to the file name of the last contaminant 
# database, if that determined the clear range to fall below the minimum value 
# (-l parameter, default 100). There are three reserved values of 
# the trash code:

#   "shortq" - assigned when the sequence length decreases 
#              below the minimum accepted length (-l) after polyA            
#              or low quality ends trimming;
# "low_qual" - assigned when the percentage of undetermined bases
#              is greater than 3% in the clear range;
#     "dust" - assigned when less than 40nt of the sequence 
#              is left unmasked by the "dust" low-complexity filter;

#   "short" is another option this field could get, that is not documented


# Use grep and awk to get this in a new file:
# echo 'grep -v  --perl "shortq|low_qual|dust" temporary_file_raw.fna.screen.cln | awk '{print $1"\t"$3"\t"$4}' > $TrimPointsFile'
# not sure how to echo this many quotes ...
    grep -v  --perl "short|shortq|low_qual|dust" $FastaFileRaw.screen.cln | awk '{print $1"\t"$3"\t"$4}' > $TrimPointsFile

 # To trim sff:
    echo "/home/ele/tools/Newbler2.5.3/DataAnalysis_2.6/sfffile  -t $TrimPointsFile -i $TrimPointsFile -o $TrimmedFileSFF $InfileSFF"
    /home/ele/tools/Newbler2.6/bin/sfffile -t $TrimPointsFile -i $TrimPointsFile -o $TrimmedFileSFF $InfileSFF;  exitOnError "sfffile"
    
    echo -e "\nFinished writing trimmed file: '$TrimmedFileSFF'.\n"

# cleaning up
    /bin/rm *$FastaFileRaw*  $TrimPointsFile  outparts_cln.sort 
    /bin/rm -r cleaning_* #seqcleans subprocess directories

done;


