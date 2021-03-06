#!/bin/tcsh -f

if ( $#argv == 0 ) then
  echo "Usage: $0 obsids.lst"
  echo "  for each obsid, downloads the relevant data files from the Chandra archive"
  echo "  and runs processing algorithms and puts everything under /data/hrc/s/"
  exit 1
endif

set currwd = `pwd`
set obslst = $currwd/$1

cd /data/hrc/s/

foreach obsid ( `cat $obslst | grep -v ^\#` )
  echo $obsid
  set numobs = `echo $obsid | awk '{print $1*1}'`

  if ( -d $obsid ) then
    if ( ! -d $obsid/analysis ) mkdir -p $obsid/analysis
    if ( ! -f $obsid/analysis/hrcf${obsid}_evt2.fits && ! -f $obsid/analysis/hrcf${obsid}_evt2.fits.gz ) then

      if ( ! -f $obsid/analysis/obsid.lst ) then
        echo "$obsid" > $obsid/analysis/obsid.lst
      endif

      if ( ! -f $obsid/analysis/.fix_rangelev+widthres ) then
        perl /data/aschrc6/wilton/isobe/Project9/Scripts/Hrc_S/fix_rangelev+widthres.pl $obsid/analysis/obsid.lst
        touch $obsid/analysis/.fix_rangelev+widthres
      endif

      if ( ! -f $obsid/analysis/.mk_new_badpix ) then
        perl /data/aschrc6/wilton/isobe/Project9/Scripts/Hrc_S/mk_new_badpix.pl $obsid/analysis/obsid.lst
        touch $obsid/analysis/.mk_new_badpix
      endif

      if ( ! -f $obsid/analysis/.run_hpe ) then
        perl /data/aschrc6/wilton/isobe/Project9/Scripts/Hrc_S/run_hpe.pl $obsid/analysis/obsid.lst
        touch $obsid/analysis/.run_hpe
      endif

      if ( ! -f $obsid/analysis/.stat+gti_filter ) then
        perl /data/aschrc6/wilton/isobe/Project9/Scripts/Hrc_S/stat+gti_filter.pl $obsid/analysis/obsid.lst
        touch $obsid/analysis/.stat+gti_filter
      endif

      if ( ! -f $obsid/analysis/.fix_dtcor ) then
        perl /data/aschrc6/wilton/isobe/Project9/Scripts/Hrc_S/fix_dtcor.pl $obsid/analysis/obsid.lst
        touch $obsid/analysis/.fix_dtcor
      endif

      chgrp -R mta /data/hrc/s/$obsid
      chmod -R a+r,g+w /data/hrc/s/$obsid
      find /data/hrc/s/$obsid -type d -exec chmod og+rx {} \;
    endif
  endif
end
