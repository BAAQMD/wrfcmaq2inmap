#!/bin/csh -f

# Current work path that of wrfcmaq2inmap_2018
set work = ~/Downloads/InMap_run/wrfcmaq2inmap_2018

# Path to daily WRF outputs
#set wrfpath = /path/to/WRF/v3.4_12US_2011_35L_NLCD40/GHRSST_gzipped_wrfouts
set wrfpath = ~/Downloads/InMap_run

# Path to MCIP outputs -- just need DENS variable from daily METCRO3D
#set mcippath = /path/to/mcip/metcro3d/2011_GHRSST
set mcippath = ~/Downloads/InMap_run

# Path to output from CMAQ -- need daily concentration files
#set cmaqpath = /path/to/camq/conc/12US2/2011/output
set cmaqpath = ~/Downloads/InMap_run

set vertmap  = $work/ancillary/vert_layers_50_to_28.csv
set griddesc = $work/ancillary/griddesc.txt

# Year to process
set year = 2018

# Starting julian date
set nday = 0
set prevday = -1
set chunk = "2017-12-29"

mkdir $work/temp
mkdir $work/daily

while ($nday <= 364)
    set ymd = `date -d "${nday} days ${year}-01-01" +"%Y%m%d"`
    echo $ymd
    set wrfdate = "`echo $ymd | cut -c-4`-`echo $ymd | cut -c5-6`-`echo $ymd | cut -c7-8`"
    echo $wrfdate

    set mcip = "${mcippath}/METCRO3D_${ymd}.nc"
    set cmaq = "${cmaqpath}/CCTM_CONC_v533_pgi_saprc07tcx_ae6_aq_r2_base18_${ymd}.nc"

    set prevymd = `date -d "${prevday} days ${year}-01-01" +"%Y%m%d"`
    set prevdate = "`echo $prevymd | cut -c-4`-`echo $prevymd | cut -c5-6`-`echo $prevymd | cut -c7-8`"

    if (($prevday % 4 == 0)) then
        echo "chunk"
        set chunk = $prevdate
        echo $chunk
    endif

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

    ncrcat -d Time,12,35 $wrffile1 $wrffile2 $wrfout

    set outfile = $work/daily/wrfcmaq_$wrfdate.ncf

    $work/ancillary/gen_wrfcmaq_2018_1km_batch.py $wrfout $mcip $cmaq $ymd $outfile $vertmap $griddesc
  
    @ nday++
    @ prevday++
end #while
