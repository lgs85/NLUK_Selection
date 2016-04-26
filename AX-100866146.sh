#!/bin/bash


#First get AX-100866146 data for all NLUK individuals
plink --bfile ../Great\ tit\ HapMap/Plink/NLUK/NLUK --extract AX-100866146.txt --recode --out AX-100866146NLUK --autosome-num 36

#Now do the same for the most recent Dutch samples
plink --file NL_bill_recent --extract AX-100866146.txt --recode --out AX-100866146NLrecent --autosome-num 36 --allow-extra-chr
