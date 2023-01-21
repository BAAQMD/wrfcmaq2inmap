#!/bin/csh -f

# Load more recent netcdf for ncrcat
source ~/.module
module unload library/netcdf-4.1.3
module load library/netcdf-4.4.1.1

# Current work path that of wrfcmaq2inmap_2018
set work = $cwd:h

# Path to daily WRF outputs
set wrfpath = /wrk2/bkoo/WRFOUT/wrf-v4.1/wrfout_2018_NAM_NLCD/base

# Path to MCIP outputs -- just need DENS variable from daily METCRO3D
set mcippath = /wrk2/bkoo/CMAQv5.3/CMAQ-5.3.1/data/mcip4.5_d04/wrf4.1_NAM_NLCD_2018

# Path to output from CMAQ -- need daily concentration files
set cmaqpath = /scratch/bkoo/cmaq/outputs/CCTM_v531_pgi_saprc07tcx_ae6_aq_newngc2_base18

set vertmap  = $work/ancillary/vert_layers_50_to_28.csv
set griddesc = $work/ancillary/griddesc.txt
set chunkmap = $work/scripts/CHUNK_DATE_MAPPING_2018.csv

# Year to process
set year = 2018

set start_mm = 01
set start_dd = 01
set num_days = 365

# Starting julian date
set nday = 0
@ prevday = $nday - 1

mkdir -p $work/temp
mkdir -p $work/daily

while ($nday < $num_days)
    echo "CLOCK_BEG: `date`"

    set ymd = `date -d "${nday} days ${year}-${start_mm}-${start_dd}" +"%Y%m%d"`
    echo $ymd
    set wrfdate = "`echo $ymd | cut -c-4`-`echo $ymd | cut -c5-6`-`echo $ymd | cut -c7-8`"
    echo $wrfdate

    set mcip = "${mcippath}/METCRO3D_${ymd}.nc"
    set cmaq = "${cmaqpath}/CCTM_CONC_v531_pgi_saprc07tcx_ae6_aq_newngc2_base18_${ymd}.nc"

    set prevymd = `date -d "${prevday} days ${year}-${start_mm}-${start_dd}" +"%Y%m%d"`
    set prevdate = "`echo $prevymd | cut -c-4`-`echo $prevymd | cut -c5-6`-`echo $prevymd | cut -c7-8`"

    set chunk = `grep "${ymd}," ${chunkmap} | cut -d"," -f2`
    echo "chunk: $chunk"

    set wrf = wrfout_d04_${wrfdate}_00:00:00
    set wrf1 = wrfout_d04_${prevdate}_12:00:00
    set wrf2 = wrfout_d04_${wrfdate}_12:00:00

    if ( -e "$wrfpath/${chunk}/${wrf1}.gz") then
        echo "start unzip"
        gunzip -k $wrfpath/${chunk}/${wrf1}.gz
    else
        echo "already unzipped"
    endif

    if ( -e "$wrfpath/${chunk}/${wrf2}.gz") then
        echo "start unzip"
        gunzip -k $wrfpath/${chunk}/${wrf2}.gz
    else
        echo "already unzipped"
    endif

    set wrfout = $work/temp/${wrf}
    set wrffile1 = $wrfpath/${chunk}/${wrf1}
    set wrffile2 = $wrfpath/${chunk}/${wrf2}
    rm -f $wrfout

    ncrcat -d Time,12,35 $wrffile1 $wrffile2 $wrfout

    set outfile = $work/daily/wrfcmaq_$wrfdate.ncf
    rm -f $outfile

    $work/ancillary/gen_wrfcmaq_2018_1km_batch.py $wrfout $mcip $cmaq $ymd $outfile $vertmap $griddesc

    echo "CLOCK_END: `date`"
  
    @ nday++
    @ prevday++
end #while
