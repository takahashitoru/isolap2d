#!/bin/bash

RUN=1 # Execute the binary
#RUN=0 # Does not execute the binary (dryrun); this is useful to create errfile 

MAKEFILE=Makefile
#MAKEFILE=Makefile_s.nuem

METHOD=CONV
#METHOD=FMM

if [ ${METHOD} = "CONV" ] ; then
    method=conv
elif [ ${METHOD} = "FMM" ] ; then
    method=fmm
else
    echo "Invalid METHOD"
    exit
fi


#KIND=ORIGINAL
#KIND=CBIE
#KIND

if [ "${KIND}" = "ORIGINAL" ] ; then
    EXE=${method}_org
    FLAG=-D_DUMMY_
elif [ "${KIND}" = "CBIE" ] ; then
    EXE=${method}_cbie
    FLAG=-DSET_FREE_TERM
else
    EXE=${method}
    KIND=_DUMMY_
    FLAG=-DSET_FREE_TERM
fi

#for P in 2 3 4
for P in 2 3 4
do
    
#   for NG in 1 3 4 6 9 16
#    for NG in 4 8 12
#    for NG in 4 8 12
#    for NG in 4
    for NG in 4 6 8 10 12
    do

#	for MAXEPC in 30
	for MAXEPC in 0 # dummy
	do

#	    for NTERM in 10 20 30
#	    for NTERM in 20 30 40
	    for NTERM in 0 # dummy
	    do

		if [ ${METHOD} = "CONV" ] ; then
		    errfile=${EXE}-p${P}ng${NG}.err
		elif [ ${METHOD} = "FMM" ] ; then
		    errfile=${EXE}-p${P}ng${NG}epc${MAXEPC}nterm${NTERM}.err
		fi

		rm -f ${errfile}

		#for CP in 100 200 400 800 1600 3200 6400
		#for CP in 1000 2000 4000 8000 16000 32000 64000
		#for CP in 1000 2000 4000 8000 16000
		#for CP in 1000 4000 16000 64000 256000
		#for CP in 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576
		for CP in 256 512 1024 2048 4096 8192 16384 32768
		#for CP in 32768 65536 131072 262144 524288 1048576
		do
		    
		    if [ ${RUN} = 1 ] ; then
			WORK_PLACE=$(pwd)
			cd ..
			make clean
			make -f ${MAKEFILE} -j EXE=${EXE}.out METHOD=${METHOD} KIND=${KIND} FLAG=${FLAG} P=${P} NG=${NG} MAXEPC=${MAXEPC} NTERM=${NTERM} CP=${CP}
			cd ${WORK_PLACE}
		    fi
		    
		    if [ ${METHOD} = "CONV" ] ; then
			tail=-p${P}ng${NG}cp${CP}
		    elif [ ${METHOD} = "FMM" ] ; then
			tail=-p${P}ng${NG}epc${MAXEPC}nterm${NTERM}cp${CP}
		    fi

		    resfile=${EXE}${tail}.res
		    logfile=${EXE}${tail}.log

		    if [ ${RUN} = 1 ] ; then
			rm -f ${resfile} ${logfile}
			(../${EXE}.out > ${resfile}) >& ${logfile}
		    fi

		    if [ -f ${resfile} ] && [ -f ${logfile} ] ; then

			L2error=$(grep L2error ${logfile} | cut -d"=" -f2)
			Maxerror=$(grep Maxerror ${logfile} | cut -d"=" -f2)
			residual=$(grep -e "resid = " ${logfile} | cut -d"=" -f2)
			nmv=$(grep -e "nmv   = " ${logfile} | cut -d"=" -f2) # larger by one than niter
			timer_all=$(grep -e "timer_all = " ${logfile} | cut -d"=" -f2)
			ave_mv=$(echo "${timer_all}/${nmv}" | bc -l)

			echo ${CP} ${L2error} ${Maxerror} ${residual} ${nmv} ${timer_all} ${ave_mv} >> ${errfile}

		    fi

		done
	    done
	done
    done
done
