#! /usr/bin/perl
   use List::Util qw[min max];
#--------------------------------------------------------------------------
#
# adas8#3    :  ADAS8#3 - Automated R-matrix calculations
#
# AUTHOR     :  Allan Whiteford, Mike Witthoeft, Stuart Loch and Connor Ballance
# DATE       :  26-09-2005
# 
# MODIFIED:
#
#       1.1     Allan Whiteford
#		26-09-2005 
#               First Release
#       1.2     Allan Whiteford
#               Fixed bug where the month was wrong by one.
#       1.3     Allan Whiteford
#               Allowed --outer=N for spefification of chunk.
#		Only re-run structure if necessary.
#       1.4     Allan Whiteford
#               Fixed bug in maximum allowed occupation
#		Added --example output
#       1.5     Allan Whiteford
#               Added --root option
#		Made --archive option more sensible
#
#       2011 :
#
#       1.6     Mike Witthoeft, Stuart Loch and Connor Ballance
#               Added --dip option for calculation D- and B- matrices
#                for the general case
#               Added OMEGA filtering for numerical failures
#               Replaced omgmerge with arrange       
#
#      2012
#
#       1.7     Connor Ballance and Stuart Loch
#               Now an 'accelerated' parallel approach in which 
#               every partial wave and diagonalisation is carried 
#               out concurrently.
#
#               Added utility codes : hfilter,dfilter (no input required)
#     
#
# VERSION:
#       1.1    26-09-2005 
#       1.2    22-01-2007 
#	1.3    09-07-2007
#	1.4    11-08-2008
#	1.5    20-08-2008
#       1.6    19-09-2011
#--------------------------------------------------------------------------

sub auto_maxc
{
	mkdir ("maxc", 0755);
	chdir ("maxc");
	symlink("../str/radout","radial");

#	return;

	$try=80;

        open(FD,"> dstg1") or die("Can't open max/dstg1");
        print FD "S.S. Automatically generated ADAS8#3 stg1 input\n";
        print FD " &STG1A RELOP='MVD' &END\n";
        print FD " &STG1B MAXLA=$maxla MAXLT=$maxlt MAXC=$try MAXE=$maxep &END\n";
	close(FD);

#	system("$stg1ex_exec < dstg1");
	$_=`tail -1 rout1r`;
	if (!/CPU/)
	{
		die ("Basis orbital calculation stg1 calculation didn't finish");
	}

	open(FD,"rout1r");

	while(<FD>)
	{
		if (/EIGENVALUE/)
		{
			$_=<FD>;
			$i=0;
			$index=0;
			while($_=<FD> and ! /BUTTLE/)
			{
				$i++;
				if (/^ +[-,\.,0-9]* +([\.,0-9]*) +[0-9]*/)
				{
					if ($1 > $maxe*2 and $index==0)
					{
						$index=$i;
					}
				}

			}
			print $index;
		}
	}
	close(FD);
	chdir("..");
}

sub make_dirs
{
	if (!-d "str"){ mkdir("str", 0755) || die "Cannot make directory"; }
	if (!-d "auger"){ mkdir("auger", 0755) || die "Cannot make directory"; }
	if (!-d "inner"){ mkdir("inner", 0755) || die "Cannot make directory"; }
	if (!-d "outer"){ mkdir("outer", 0755) || die "Cannot make directory"; }
	if (!-d "nonx"){ mkdir("nonx", 0755) || die "Cannot make directory"; }
	if (!-d "dip" && $rad_damp == 1){ mkdir("dip", 0755) || die "Cannot make directory"; }
	if (!-d "tcc"){ mkdir("tcc", 0755) || die "Cannot make directory"; }
	if (!-d "born"){ mkdir("born", 0755) || die "Cannot make directory"; }
	if (!-d "adas"){ mkdir("adas", 0755) || die "Cannot make directory"; }

	if (!-e "str/radial") { open(FD,"> str/radout"); close(FD); }
	if (!-e "inner/H.DAT") { open(FD,"> inner/H.DAT"); close(FD); }
	if (!-e "inner/NX1.DAT") { open(FD,"> inner/NX1.DAT"); close(FD); }
	if (!-e "inner/NX2.DAT") { open(FD,"> inner/NX2.DAT"); close(FD); }
	if (!-e "tcc/TCCDW.DAT") { open(FD,"> tcc/TCCDW.DAT"); close(FD); }

	if (!-e "inner/radial") { symlink("../str/radout","inner/radial")|| die "Cannot make symlink"; };
	if ($rad_damp == 1 && !-e "dip/radial") { symlink("../str/radout","dip/radial")|| die "Cannot make symlink"; };
	if ($rad_damp == 1 && !-e "outer/H.DAT") { symlink("../inner/H.DAT","outer/H.DAT")|| die "Cannot make symlink"; };
	if (!-e "tcc/radial") { symlink("../str/radout","tcc/radial")|| die "Cannot make symlink"; };
	if (!-e "nonx/TCCDW.DAT") { symlink("../tcc/TCCDW.DAT","nonx/TCCDW.DAT")|| die "Cannot make symlink"; };
	if (!-e "nonx/NX1.DAT") { symlink("../inner/NX1.DAT","nonx/NX1.DAT")|| die "Cannot make symlink"; };
	if (!-e "nonx/NX2.DAT") { symlink("../inner/NX2.DAT","nonx/NX2.DAT")|| die "Cannot make symlink"; };
}

sub write_inputs_inner
{
        $radeq="";
        open(FD,"> inner/dstg1") or die("Can't open inner/dstg1");
        print FD "S.S. Automatically generated ADAS8#3 stg1 input\n";
        print FD " &STG1A RELOP='MVD'$radeq &END\n";
        print FD " &STG1B MAXLA=$maxla MAXLT=$maxlt MAXC=$maxc MAXE=$maxep &END\n";
	close(FD);

	$maxorb=@orbvals;
	$nast=@TERMSFILE;
        
        if ($stg2ex_ncpu > 1)
        {
                $nprocstg1="NPROCSTG1=$stg1ex_ncpu";
        }
        open(FD,"> inner/dstg2") or die("Can't open inner/dstg2");
        print FD "S.S. Automatically generated ADAS8#3 stg2 input\n";
        print FD " &STG2A RELOP='MVD'$radeq ISORT=1 $nprocstg1 &END\n";
        print FD " &STG2B LNOEX=18 MAXORB=$maxorb NELC=$nelec NAST=$nast INAST=0 MINLT=0 MAXLT=$maxlt MINST=$minst MAXST=$maxst &END\n";


	foreach (@orbvals)
	{
        	print FD " @$_[0] @$_[1]";
        }
       	print FD "\n";
        
        $mxconf=@occ_number;
	print FD  ("  $mxconf\n ");

	for ($i=0;$i<@orbvals;$i++)
	{
		printf FD ("%2d",$min_orb[$i]);
	}
	print FD "\n ";
	for ($i=0;$i<@orbvals;$i++)
	{
		printf FD ("%2d",$max_orb[$i]);
	}
	print FD "\n";

 	foreach (@occ_number)
        {
        	@temp=@$_;
                print FD " ";
        	for ($i=0;$i<@orbvals;$i++)
                {
                	printf FD ("%2d",$temp[$i]);
                }
		print FD "  0\n";
        }
        print FD @TERMSFILE;
        
	write_nplus1();
               
	close(FD);

	$sqrt_cpu=sqrt($stg3ex_ncpu);
	$rem = $sqrt_cpu % 1;

	if (int($sqrt_cpu)==$sqrt_cpu)
	{
		$nprow=$sqrt_cpu;
		$npcol=$sqrt_cpu;
	}
	else
	{
		$nprow=int($sqrt_cpu);
        	while(int($stg3ex_ncpu / $nprow)*$nprow != $stg3ex_ncpu)
	        {
        		$nprow--;
	        }
		$npcol=$stg3ex_ncpu / $nprow
	}

        open(FD,"> inner/dstg3") or die("Can't open inner/dstg3");
        print FD "S.S. Automatically generated ADAS8#3 stg3 input\n";
        print FD " &STG3A$radeq &END\n";
        print FD " &STG3B INAST=0 NAST=0 &END\n";
	if ($nprow!=1 or $npcol!=1)
	{
	        print FD " &MATRIXDAT NPROW=$nprow NPCOL=$npcol &END\n";
	}
	close(FD);

}


sub write_inputs_dip
{
  if($rad_damp == 0) { return; }

  # determine parameters for input files
  $ming=99;
  $maxg=0;
  $minl=99;
  $maxl=0;
  foreach(@TERMSFILE)
  {
    @arr=split(/ */,$_);
    $g=@arr[1];
    $l=@arr[2];
    if($g < $ming) { $ming=$g-1; }
    if($ming == 0) { $ming=$ming+2; }
    if($g > $maxg) { $maxg=$g+1; }
    if($l < $minl) { $minl=$l; }
    if($l > $maxl) { $maxl=$l; }
  }
  $maxlt_dip=$orbl+$maxl;

  # determine how many processors pstg2_dip needs to have
  $nproc2_req=3*($maxl-$minl);
  if($minl > 0) { $nproc2_req++; }
#  system("echo $ming $maxg");
  $nproc2_req*=($maxg-$ming)/2+1;

  # write dstg1
  $radeq=" RAD='YES'";
  open(FD,"> dip/dstg1") or die("Can't open dip/dstg1");
  print FD "S.S. Automatically generated ADAS8#3 stg1 input\n";
  print FD " &STG1A RELOP='MVD'$radeq &END\n";
  print FD " &STG1B MAXLA=$maxla MAXLT=$maxlt_dip MAXC=$maxc MAXE=$maxep &END\n";
 	close(FD);

	 $maxorb=@orbvals;
	 $nast=@TERMSFILE;
        
  if ($stg2dip_ncpu > 1)
  {
    if($stg2dip_ncpu != $nproc2_req && $rad_damp == 1) { die("PSTG2R_DIP must have $nproc2_req processors"); }
    $nprocstg1="NPROCSTG1=$stg1dip_ncpu";
  }
  open(FD,"> dip/dstg2") or die("Can't open dip/dstg2");
  print FD "S.S. Automatically generated ADAS8#3 stg2 input\n";
  print FD " &STG2A RELOP='MVD'$radeq ISORT=1 $nprocstg1 &END\n";
  print FD " &STG2B MAXORB=$maxorb NELC=$nelec NAST=$nast INAST=0 MINLT=0 MAXLT=$maxlt_dip MINST=$minst MAXST=$maxst &END\n";
 	foreach (@orbvals)
 	{
   	print FD " @$_[0] @$_[1]";
  }
  print FD "\n";

  $mxconf=@occ_number;
	 print FD  ("  $mxconf\n ");

	 for ($i=0;$i<@orbvals;$i++)
	 {
	  	printf FD ("%2d",$min_orb[$i]);
	 }
	 print FD "\n ";
	 for ($i=0;$i<@orbvals;$i++)
	 {
	  	printf FD ("%2d",$max_orb[$i]);
	 }
	 print FD "\n";

 	foreach (@occ_number)
  {
   	@temp=@$_;
    print FD " ";
  	 for ($i=0;$i<@orbvals;$i++)
    {
      printf FD ("%2d",$temp[$i]);
    }
		  print FD "  0\n";
  }
  print FD @TERMSFILE;
        
	 write_nplus1();
               
	 close(FD);

 	$sqrt_cpu=sqrt($stg3dip_ncpu);
 	$rem = $sqrt_cpu % 1;

 	if (int($sqrt_cpu)==$sqrt_cpu)
 	{
	  	$nprow=$sqrt_cpu;
	  	$npcol=$sqrt_cpu;
	 }
	 else
	 {
		  $nprow=int($sqrt_cpu);
    while(int($stg3dip_ncpu / $nprow)*$nprow != $stg3dip_ncpu)
	   {
     	$nprow--;
	   }
		  $npcol=$stg3dip_ncpu / $nprow
 	}

  open(FD,"> dip/dstg3") or die("Can't open dip/dstg3");
  print FD "S.S. Automatically generated ADAS8#3 stg3 input\n";
  print FD " &STG3A $radeq &END\n";
  print FD " &STG3B INAST=0 NAST=0 &END\n";
 	if ($nprow!=1 or $npcol!=1)
 	{
	   print FD " &MATRIXDAT NPROW=$nprow NPCOL=$npcol &END\n";
	 }
	 close(FD);

}


sub write_inputs_stgb
{
  if($rad_damp == 0) { return; }

  $totalboundsym ==0;

  # determine parameters for input file
  $maxn=0;
  $orbl=0;
  foreach(@orbvals)
  {
    if(@$_[1] > $orbl) { $orbl=@$_[1]; }
    if(@$_[0] > $maxn) { $maxn=@$_[0]; }
  }

  $minn=99;
 	foreach (@occ_number)
  {
   	@temp=@$_;
  	 for ($i=0;$i<@orbvals;$i++)
    {
      $n=@orbvals[$i]->[0];
      $l=@orbvals[$i]->[1];
      $ne=$temp[$i];
      $nchk=4*@$_[1]+1;
      if($ne < $nchk) { last; }
    }
    if($n < $minn) { $minn=$n; }
  }
  $minn=$minn-0.5;
  $maxn=$maxn+0.5;

  $ming=99;
  $maxg=0;
  $minl=99;
  $maxl=0;
  foreach(@TERMSFILE)
  {
    @arr=split(/ */,$_);
    $g=@arr[1];
    $l=@arr[2];
    if($g < $ming) { $ming=$g; }
    if($g > $maxg) { $maxg=$g; }
    if($l < $minl) { $minl=$l; }
    if($l > $maxl) { $maxl=$l; }
  }
  $maxlt_dip=$orbl+$maxl;

  $ming=$ming-1;
  if($ming == 0) { $ming=$ming+2; }
  $maxg=$maxg+1;
  $minl=$minl-1;
  if($minl < 0) { $minl=0; }
  $maxl=$maxlt_dip+1;

  if(!-e "inner/sizeH.dat")
  {
    die("sizeH.dat file does not exist; run inner region first\n");
  }
  open(FD,"inner/sizeH.dat") or die("unable to open sizeH.dat\n");
  @dat=<FD>;
  close(FD);

  open(FD,"> outer/dstgb") or die("Can't open outer/dstgb");
  print FD " &STGB IOPT2=1 IRAD=1 IPERT=0 IPRINT=0 &END\n";
  foreach(@dat)
  {
    @arr=split(/\s+/);
    $p=pop(@arr);
    $l=pop(@arr);
    $g=pop(@arr);
    $mnp1=pop(@arr);
    next if $mnp1 == 0;
    next if $g < $ming || $g > $maxg;
    next if $l < $minl || $l > $maxl;
    print FD " $g  $l  $p\n";
    print FD "$minn $maxn 0.001\n";
    $totalboundsym=$totalboundsym + 1;
  }
  print FD "-1 -1 -1\n";
  close(FD);
  system("echo Total number of stgb symmetries equals $totalboundsym ");
}

sub write_omegafilter
{
	chdir("outer");
	open(FD,"> dadd") or die("Can't open outer/dadd");
	print FD " &SADD YMULT=1000.0 IBIGE=0 &END\n";
	close(FD);
	symlink("OMEGAZ","omadd1");
	system("$omadd_exec");
	unlink("OMEGAZ");
	system("mv omaddt omadd1");
	system("$omadd_exec");
	system("mv omaddt OMEGAZ-YMULT1000");
        system("rm omadd1 omadd2 oadd oadd2");
}

sub write_omegafilternon
{
	chdir("nonx");
	open(FD,"> dadd") or die("Can't open outer/dadd");
	print FD " &SADD YMULT=5.0 IBIGE=0 &END\n";
	close(FD);
	symlink("OMEGAZ","omadd1");
	system("$omadd_exec");
	unlink("OMEGAZ");
	system("mv omaddt omadd1");
	system("$omadd_exec");
	system("mv omaddt OMEGAZ-NON-YMULT5");
        system("rm omadd1 omadd2 oadd oadd2");
}


sub write_inputs_tcc
{

        open(FD,"> tcc/dstg1") or die("Can't open tcc/dstg1");
        print FD "S.S. Automatically generated ADAS8#3 stg1 input\n";
        print FD " &STG1A RELOP='TCC' IDWOUT=2 &END\n";
        print FD " &STG1B MAXLA=$maxla MAXLT=$maxlt MAXC=$maxc MAXE=$maxep &END\n";
	close(FD);


	$maxorb=@orbvals;
	$nast=@TERMSFILE;
        
        open(FD,"> tcc/dstg2") or die("Can't open tcc/dstg2");
        print FD "S.S. Automatically generated ADAS8#3 stg2 input\n";
        print FD " &STG2A RELOP='TCC' IDWOUT=2 &END\n";
        print FD " &STG2B MAXORB=$maxorb NELC=$nelec NAST=$nast INAST=0 MINLT=0 MAXLT=$maxlt MINST=$minst MAXST=$maxst &END\n";


	foreach (@orbvals)
	{
        	print FD " @$_[0] @$_[1]";
        }
       	print FD "\n";
        
        $mxconf=@occ_number;
	print FD  ("  $mxconf\n ");

	for ($i=0;$i<@orbvals;$i++)
	{
		printf FD ("%2d",$min_orb[$i]);
	}
	print FD "\n ";
	for ($i=0;$i<@orbvals;$i++)
	{
		printf FD ("%2d",$max_orb[$i]);
	}
	print FD "\n";

 	foreach (@occ_number)
        {
        	@temp=@$_;
                print FD " ";
        	for ($i=0;$i<@orbvals;$i++)
                {
                	printf FD ("%2d",$temp[$i]);
                }
		print FD "  0\n";
        }
        print FD @TERMSFILE;
        
	write_nplus1();
       
	close(FD);



	$jnast=@LEVELSFILE;


	open(FD,"> tcc/dstgjk") or die("Can't open tcc/dstgjk");
        print FD "S.S. Automatically generated ADAS8#3 stgjk input\n";
        print FD " &STGJA IDWOUT=2 RELOP='TCC' &END\n";
        print FD " &STGJB JNAST=$jnast &END\n";
        print FD @LEVELSFILE;
	close(FD);



}

sub write_nplus1
{

#	print FD "1\n";
#	print FD " 2 0 4 0 0 0\n";
#	print FD " 2 2 6 2 2 2\n";
#	print FD " 2 2 6 0 0 0  2\n";

#	print FD "1\n";
#	print FD " 2 0 0 0 0 0\n";
#	print FD " 2 2 2 2 2 2\n";
#	print FD " 2 1 0 0 0 0  2\n";

	@nplus_number=();	

 	foreach (@occ_number)
        {                
               @temp_hold=@$_;
               for ($i=0;$i<@orbvals;$i++)
               {
			@temp=@temp_hold;

                       	if ($temp[$i] < $maxpaul[$i] and $openshell[$i]==1)
                       	{
                               $temp[$i]=$temp[$i]+1;
                           
                               #print "@temp\n";
                               push @nplus_number, [ @temp ];
                       	}
               }
        }
        foreach (@nplus_number)
        {
               @temp=@$_;
               for ($i=0;$i<@orbvals;$i++)
               {
                       
                        if ($temp[$i] > $maxplus_orb[$i])
                        {
                               $maxplus_orb[$i]=$temp[$i];
                        }
                       
                        if ($minplus_orb[$i] > $temp[$i])
                        {
                               $minplus_orb[$i]=$temp[$i];
                        }
               }
        }
        
	$j=0;
 	foreach (@nplus_number)
        {
        	@temp=@$_;
		$nplus_string[$j]="";		
        	for ($i=0;$i<@orbvals;$i++)
                {
                	$nplus_string[$j] .= sprintf ("%2d",$temp[$i]);
                }
		$j++;
        }
	
	%seen = ();
	@nplus_uniq=();
	foreach $item (@nplus_string) {
	    push(@nplus_uniq, $item) unless $seen{$item}++;
	}
	
        $mxconf=@nplus_uniq;
	
        print FD  ("$mxconf\n");

	for ($i=0;$i<@orbvals;$i++)
	{
		printf FD ("%2d",$minplus_orb[$i]);
	}
	print FD "\n";
	for ($i=0;$i<@orbvals;$i++)
	{
		printf FD ("%2d",$maxplus_orb[$i]);
	}
	print FD "\n";


 	foreach (@nplus_uniq)
        {
		printf FD "$_  0\n";
#        	@temp=@$_;
#        	for ($i=0;$i<@orbvals;$i++)
#                {
#                	printf FD ("%2d",$temp[$i]);
#                }
#		print FD "  0\n";
        }

}

sub make_outer_sub
{
	mkdir("outer/$dir", 0755);
	open (FD,"> outer/$dir/dstgf");

        if($rad_damp == 1)
        { 
	print FD " &STGF IMESH=1 IQDT=2 PERT='YES' IPRINT=-2 IAUGER=$auger_damp \n ";
        }
        else
        {
        print FD " &STGF IMESH=1 IQDT=2 PERT='YES' IPRINT=-2\n ";
        }

	print FD " IPRKM=4 IRD0=101 NOMWRT=50 &END\n";
                
        if ($use_inter == 1)
        {
	        print FD " &MESH1 MXE=$points E0=$start EINCR=$eincr NSEQ=$nseq &END\n";
	}
        else
        {
 	        print FD " &MESH1 MXE=$points E0=$start EINCR=$eincr &END\n";
       
        }
        close (FD);

	open (FD,"> outer/$dir/dstgicf");
	print FD "C \n";
	
        if ($use_inter == 1)
        {
                $imode=-1;
        }
        else
        {
                $imode=0;
        }
        print FD " &STGIC ITCC=1 IMODE=$imode INOEXCH=0 NOMWRT=9999999 IRD0=50 PRINT='FORM' &END\n";
        
        if ($use_inter == 1)
        {
	        print FD " &MESH1 IEQ=$interpolation &END\n";
	}
        else
        {
	        print FD " &MESH1 MXE=$points E0=$start EINCR=$eincr &END\n";
        }

	close (FD);

	symlink("../../inner/H.DAT","outer/$dir/H.DAT");
	symlink("../../tcc/TCCDW.DAT","outer/$dir/TCCDW.DAT");
}

sub write_inputs_outer
{
	$points_outer=int($chunk_size);
        $mesh_fine=($last_thres_z2 - $first_thres_z2)/$points_outer;
	$running_total=0;

        system("echo OUTER EX fine mesh: $points_outer,$first_thres_z2,$last_thres_z2");
	$dir=1;
	$start=$first_thres_z2+$mesh_fine;
	$points=$chunk_size;
	$eincr=$mesh_fine;
        $emax=$start+$eincr*$points;
#        system("echo OUTER EX fine mesh: $start,$points,$emax");
	make_outer_sub();
	$running_total=$running_total+$chunk_size;
	
	$dir=$dir+1;

	$eincr=$mesh_fine*100;
        $mesh_coarse=$eincr;

	$start=$start+$points*$mesh_fine+$mesh_fine;
#        system("echo OUTER EX coarse mesh: $maxe_z2,$start,$mesh_coarse");

	$points=int((($maxe_z2-$start)/$mesh_coarse));
        $intmp=int($points/$stgfex_ncpu)+1;
        $points=$intmp*$stgfex_ncpu;
        
        $points=int($chunk_size_high_ex);
        $eincr=($maxe_z2-$start)/$points;
        $emax=$start+$eincr*$points;
        system("echo OUTER EX coarse mesh: $points,$maxe_z2,$start,$mesh_coarse");
#        system("echo exchange,$start,$points,$emax");
	make_outer_sub();
}

sub write_inputs_nonx
{


        $radeq="";
        open(FD,"> nonx/dstg1") or die("Can't open nonx/dstg1");
        print FD "S.S. Automatically generated ADAS8#3 stg1 input\n";
        print FD " &STG1A RELOP='MVD'$radeq &END\n";
        print FD " &STG1B MAXLA=$maxla MAXLT=$maxlt_nx MAXC=$maxc MAXE=$maxep &END\n";
	close(FD);

	$maxorb=@orbvals;
	$nast=@TERMSFILE;
        
        if ($stg2nx_ncpu > 1)
        {
                $nprocstg1="NPROCSTG1=$stg1nx_ncpu";
        }
        open(FD,"> nonx/dstg2") or die("Can't open nonx/dstg2");
        print FD "S.S. Automatically generated ADAS8#3 stg2 input\n";
        print FD " &STG2A RELOP='MVD'$radeq ISORT=1 $nprocstg1 &END\n";
        print FD " &STG2B LNOEX=13 MAXORB=$maxorb NELC=$nelec NAST=$nast INAST=0 MINLT=$minlt_nx MAXLT=$maxlt_nx MINST=$minst MAXST=$maxst &END\n";
#        print FD " &STG2B LNOEX=10 MAXORB=$maxorb NELC=$nelec NAST=$nast INAST=0 MINLT=$minlt_nx MAXLT=$maxlt_nx MINST=$minst MAXST=$maxst &END\n";


	foreach (@orbvals)
	{
        	print FD " @$_[0] @$_[1]";
        }
       	print FD "\n";
        
        $mxconf=@occ_number;
	print FD  ("  $mxconf\n ");

	for ($i=0;$i<@orbvals;$i++)
	{
		printf FD ("%2d",$min_orb[$i]);
	}
	print FD "\n ";
	for ($i=0;$i<@orbvals;$i++)
	{
		printf FD ("%2d",$max_orb[$i]);
	}
	print FD "\n";

 	foreach (@occ_number)
        {
        	@temp=@$_;
                print FD " ";
        	for ($i=0;$i<@orbvals;$i++)
                {
                	printf FD ("%2d",$temp[$i]);
                }
		print FD "  0\n";
        }
        print FD @TERMSFILE;
        
	write_nplus1();
               
	close(FD);

	$sqrt_cpu=sqrt($stg3ex_ncpu);
	$rem = $sqrt_cpu % 1;

	if (int($sqrt_cpu)==$sqrt_cpu)
	{
		$nprow=$sqrt_cpu;
		$npcol=$sqrt_cpu;
	}
	else
	{
		$nprow=int($sqrt_cpu);
        	while(int($stg3ex_ncpu / $nprow)*$nprow != $stg3ex_ncpu)
	        {
        		$nprow--;
	        }
		$npcol=$stg3ex_ncpu / $nprow
	}

        open(FD,"> nonx/dstg3") or die("Can't open nonx/dstg3");
        print FD "S.S. Automatically generated ADAS8#3 stg3 input\n";
        print FD " &prediag npw_per_subworld=152 &END\n ";
        print FD " &STG3A$radeq &END\n";
        print FD " &STG3B INAST=0 NAST=0 &END\n";
#	if ($nprow!=1 or $npcol!=1)
#	{
	        print FD " &MATRIXDAT NPROW=$nprow NPCOL=$npcol &END\n";
#	}
	close(FD);
#        system("echo NPROW,$nprow,$npcol,$stg3ex_ncpu");
        
        
	$points_nx=int((($maxe_z2-$first_thres_z2)/($mesh_fine*10)));
        $intmp=int($points_nx/$stgfnx_ncpu)+1;
        $points_nx=$intmp*$stgfnx_ncpu;
        $points_nx=int($chunk_size)/100;

        $points_nx=int($chunk_size_nex);
        $mesh_coarse = ($maxe_z2-$first_thres_z2)/$points_nx;

#	$points_nx = $points_nx + $stgfnx_ncpu - $points_nx % $stgfnx_ncpu;

        open(FD,"> nonx/dstgf") or die("Can't open nonx/dstgf");
	print FD " &STGF IMESH=1 IQDT=2 PERT='YES' IPRINT=-2 IPRKM=4 &END\n";
	print FD " &MESH1 MXE=$points_nx E0=$first_thres_z2 EINCR=$mesh_coarse &END\n";
	close(FD);
        open(FD,"> nonx/dstgicf") or die("Can't open nonx/dstgicf");
	print FD "C \n";
	print FD " &STGIC ITCC=1 IMODE=0 INOEXCH=1 JMNTWO=$jmntwo LRGLAM=999\n";
        print FD " NOMWRT=9999999 PRINT='FORM' &END\n";
	close(FD);
        $emax=$first_thres_z2+$mesh_coarse*$points_nx;
#        system("echo NONEX energy mesh:,$first_thres_z2,$points_nx,$emax");

}

sub occ_numbers
{
	%lvals= ( 'S' => '0' , 'P' => '1' , 
                  'D' => '2' , 'F' => '3' ,
                  'G' => '4' , 'H' => '5' );

	foreach $config (@configs)
        {
                $i=0;
		foreach(@orbitals)
		{
                	$config =~ /.*($_)([0-9]*)/i;
                        if ($1)
                        {
                                if ($2)
                                {
                                        $number=$2;
                        	}
                                else
                                {
                                        $number=1;
                                }
                                
                                $occn[$i]=$number;
				if ($max_orb[$i] < $number) { $max_orb[$i]=$number; }
				if ($min_orb[$i] > $number) { $min_orb[$i]=$number; }
                        } 
                        else
                        {
                        	$occn[$i]=0;
				if ($min_orb[$i] > 0) { $min_orb[$i]=0; }
                        }
                	$i++;
                }
                push @occ_number, [ @occn ];
        }
	
        $i=0;
        foreach (@orbitals)
        {
        	/([0-9])([A-Z])/i;
        	push @orbvals, [$1 , $lvals{uc($2)}];
                $maxpaul[$i++]=4*$lvals{uc($2)}+2;
        }

        $allan=$occ_number[0];
        
        $nelec=0;
	foreach (@$allan)
        {
        	$nelec=$nelec+$_;
        }

	$charge=$nzion-$nelec;

        if ($charge == 0)
        {
        $charge=1;
        }        

        $first=1;
	foreach (@occ_number)
        {
        	@temp=@$_;
        	for ($i=0;$i<@orbvals;$i++)
                {
			if ($first!=1)
			{
                        	if ($temp[$i]!=$oldtemp[$i] or $forceopen[$i]==1)
				{
                                	$openshell[$i]=1;
                                }
                        }
                        if ($temp[$i]!=$maxpaul[$i] or $forceopen[$i]==1)
                        {
                                $openpaul[$i]=1;
                        }
                }
                @oldtemp=@temp;
                $first=0;
        }
        
        $foundclosed=0;
        $korb1=0;
        $korb2=0;
        for ($i=0;$i<@openpaul;$i++)
        {
                if ($openpaul[$i]==0 and $foundclosed==0)
                {
                        $korb2=$i+1;
                        $korb1=1;
                }
                else
                {
                        $foundclosed=1;
                }
                if ($foundclosed==1)
                {
                        $openpaul[$i]=1;
                }
        }

}

sub run_prelim_autos
{
	$MXVORB=@orbvals-$korb2;
	$MXCONF=@occ_number;
        
       	chdir "str" || die "Can't change directory to str";
        open(FD,"> das");
        print FD "A.S. Automatically generated ADAS8#3 AS input\n";
        print FD " &SALGEB RAD='NO' CUP='ICM' MXVORB=$MXVORB MXCONF=$MXCONF KUTSO=0 KORB1=$korb1 KORB2=$korb2 &END\n";
        $i=0;
	foreach (@orbvals)
	{
                if ($openpaul[$i++])
                {
        	        print FD "  @$_[0] @$_[1]";
                }
        }
       	print FD "\n";
	foreach (@occ_number)
        {
        	@temp=@$_;
                print FD "  ";
        	for ($i=0;$i<@orbvals;$i++)
                {
                        if ($openpaul[$i])
                        {
                        	printf FD ("%2d   ",$temp[$i]);
                        }
                }
		print FD "\n";
        }
	$nlam=@scale;
	print FD " &SMINIM  NZION=$nzion INCLUD=0 NLAM=$nlam PRINT='FORM' RADOUT='YES' &END\n";
	foreach (@scale)
	{
	        print FD "$_ " ;
	}
	print FD "\n";
        close(FD);
        system("$autos_exec < das");
	$_=`tail -1 olg`;
	if (!/CPU/)
	{
		die ("Preliminary Autostructure run didn't finish");
	}
#
#       if auger damping is on create LSM rather than ICM, but do not run
#
        if($auger_damp == 1){
        open(FD,"> das2");
        print FD "A.S. Automatically generated ADAS8#3 AS input\n";
        print FD " &SALGEB RAD='NO' CUP='LSM' MXVORB=$MXVORB MXCONF=$MXCONF KUTSO=0 KORB1=$korb1 KORB2=$korb2 &END\n";
        $i=0;
	foreach (@orbvals)
	{
                if ($openpaul[$i++])
                {
        	        print FD "  @$_[0] @$_[1]";
                }
        }
       	print FD "\n";
	foreach (@occ_number)
        {
        	@temp=@$_;
                print FD "  ";
        	for ($i=0;$i<@orbvals;$i++)
                {
                        if ($openpaul[$i])
                        {
                        	printf FD ("%2d   ",$temp[$i]);
                        }
                }
		print FD "\n";
        }
	$nlam=@scale;
	print FD " &SMINIM  NZION=$nzion INCLUD=0 NLAM=$nlam PRINT='FORM' RADOUT='YES' &END\n";
	foreach (@scale)
	{
	        print FD "$_ " ;
	}
	print FD "\n";
        close(FD);
        }
#        
	unlink("ols");
	unlink("oic");
	unlink("olsu");
	unlink("oicu");
	unlink("opic");
	unlink("opicu");
	unlink("opls");
	unlink("oplsu");
	unlink("SHFTIC");
	unlink("SHFTLS");
	unlink("hffcin");
	unlink("RESTART");
	unlink("OVRLAP");
	unlink("radwin");
	unlink("adf04ic");
	unlink("adf04ls");
	unlink("OMGINFIC");
	unlink("OMGINFLS");

	chdir("..");

}


sub run_auger
{
        $AUGERFILE="./auger/AUGER"; 

        if(($auger_damp == 0) || (-e $AUGERFILE)) { return; }

       	chdir "str" || die "must run inp before calculating auger widths";
        chdir "..";
        chdir "auger" || die "calculating the auger widths";
#
#
        system("cp ../str/das2 das");
        system("cp ../src/misc/rap*.py .");
        system("cp ../src/misc/rt*.py .");
        system("cp ../src/misc/asaug.py .");
#
        open(FD,"> rtk_proc") or die("Can't open rtk_proc");
        print FD "PATH:        ./ \n";
        print FD "AUTOS:       ../../rap/clean/asdeck.exe < %INFILE \n";
        print FD "WIDTH:       ../bin/width.x \n "; 
        close(FD);
#
        system("../../rap/clean/asdeck.exe < das");
        system("./asaug.py -s $auger_damp olg");
# 
#       cleanup
#
#   
        unlink("rap_core.pyc");
        unlink("rap_core.py");
        unlink("rtk_core.pyc");
        unlink("rtk_core.py");
        unlink("rtc_filter.py ");
        unlink("rtc_genconf.py");
        unlink("rtc_filter.py");
        unlink("rtc_mix.py");
        unlink("rtk_proc");
	unlink("ols");
	unlink("oic");
        unlink("CONFIG.DAT");
        unlink("TERMS");
        unlink("LEVELS");
	unlink("olsu");
	unlink("oicu");
	unlink("opic");
	unlink("opicu");
	unlink("opls");
	unlink("oplsu");
	unlink("SHFTIC");
	unlink("SHFTLS");
	unlink("hffcin");
	unlink("RESTART");
	unlink("OVRLAP");
	unlink("radwin");
	unlink("adf04ic");
	unlink("adf04ls");
	unlink("OMGINFIC");
	unlink("OMGINFLS");

	chdir("..");

}


sub write_inputs_born
{
	$MXVORB=@orbvals-$korb2;
	$MXCONF=@occ_number;
        
       	chdir "born" || die "Can't change directory to born";
        open(FD,"> das");
        print FD "A.S. Automatically generated ADAS8#3 AS input\n";
        print FD " &SALGEB RAD='E3' CUP='ICM' MXVORB=$MXVORB MXCONF=$MXCONF KUTSO=0 BORN='INF' KORB1=$korb1 KORB2=$korb2 &END\n";
	$i=0;
        foreach (@orbvals)
	{
                if ($openpaul[$i++])
                {
        	        print FD "  @$_[0] @$_[1]";
                }
        }
       	print FD "\n";
	foreach (@occ_number)
        {
        	@temp=@$_;
                print FD "  ";
        	for ($i=0;$i<@orbvals;$i++)
                {
                        if ($openpaul[$i])
                        {
                 	        printf FD ("%2d   ",$temp[$i]);
                        }
                }
		print FD "\n";
        }
	$nlam=@scale;
	print FD " &SMINIM  NZION=$nzion INCLUD=0 NLAM=$nlam PRINT='FORM' RADOUT='YES' &END\n";
        foreach (@scale)
	{
	        print FD "$_ " ;
	}
	print FD "\n";

        close(FD);
        chdir("..");
}



sub read_physics
{
	open(FD,"$physicsfile");
        
	$mesh_fine=0.00001;
	$mesh_coarse=0.001;
        $j2max_ex=24;
        $j2max_nx=80;
	$maxc=30;
	$maxe_ionpot=3;
	$interpolation=1;
        $rad_damp=0;
        $auger_damp=0;

        $read_cfgs=0;
        $read_orbs=0;
        
        
        while(<FD>)
        {
                if (/^\#.*/)
                {
                        next;
                }
        
		if (/maxc[^0-9]*([0-9]*)/i)
                {
                	$maxc=$1;
                }
		if (/continuum[^0-9]*([0-9]*)/i)
                {
                	$maxc=$1;
                }
		if (/maxe\/ionpot[^0-9]*([0-9]*)/i)
                {
                	$maxe_ionpot=$1;
                }
		if (/rdamp *= *([0-9]*)/)
                {
                	$rad_damp=$1;
                }
                if (/adamp *= *([0-9]*)/)
                {
                        $auger_damp=$1;
                }
                if (/accel *= *([0-9]*)/)
                {
                        $accelerated=$1;
                }
		if (/interpolation[^0-9]*([0-9]*)/i)
                {
                	$interpolation=$1;
                }
		if (/2Jmax_ex[^0-9]*([0-9]*)/i)
                {
                	$j2max_ex=$1;
                }
		if (/max([^0-9]*)exchange[^0-9]*([0-9]*)/i)
                {
                        if (!($1 =~ /non/i))
                	{
                                $j2max_ex=$2;
                        
                        }
                }
		if (/max.*non[^0-9]*([0-9]*)/i)
                {
                	$j2max_nx=$1;
                }
		if (/2Jmax_nx[^0-9]*([0-9]*)/i)
                {
                	$j2max_nx=$1;
                }
		if (/fine[^0-9]*([0-9]*\.[0-9]*)/i)
                {
                	$mesh_fine=$1;
                }
		if (/coarse[^0-9]*([0-9]*\.[0-9]*)/i)
                {
                	$mesh_coarse=$1;
                }
                
        	if (/CONF/i)
                {
                	$read_cfgs=1; $read_orbs=0; $i=0;
                        next;
                }
        	if (/SCAL/i)
                {
                	$read_orbs=1; $read_cfgs=0; $i=0;
                        next;
                }
                if (/COMM/)
                {
                	break;
                }
                
        	if (/^(\s)*$/)
                {
                        $read_orbs=0; $read_cfgs=0; $i=0;
                        next;
                }
        	if ($read_cfgs==1)
        	{
                	@configs[$i++]=$_;
                }
        	if ($read_orbs==1)
        	{
                	/([0-9][A-Z]) +\= +([\.,0-9]*).*/i;
                        #chop();
                	$orbitals[$i]=$1;
                	$scale[$i]=$2;
                        if (/OP/i) { $forceopen[$i]=1; }
                        
                        $max_orb[$i]=0;
                        $min_orb[$i]=999;
                        $maxplus_orb[$i]=0;
                        $minplus_orb[$i++]=999;
                        
                }
        }
        
	@knownorbitals=("1S","2S","2P","3S","3P","3D","4S","4P","4D","4F","5S","5P","5D","5F","5G");
        $searchstr=join('',@configs);
        $i=0;
        foreach (@knownorbitals)
        {
                if ($searchstr =~ /$_/i)
                {
                        $maxfound=$i;
                }
                $i++;
        }
        
        
        $size= @orbitals;
        
        if ($maxfound >= $size)
        {
                for ($i=$size;$i <= $maxfound;$i++)
                {
                        $orbitals[$i]=$knownorbitals[$i];
                        $scale[$i]="1.00";
                        $max_orb[$i]=0;
                        $min_orb[$i]=999;
                        $maxplus_orb[$i]=0;
                        $minplus_orb[$i]=999;
                }
        }
        
	close(FD);

}

sub read_proc
{
	open(FD,"$procfile") or die ("Can't open proc file: $procfile");

        $chunk_size=2000;
        $internal_split=10;
        $matrix_type='K';

	while(<FD>)
 {
		if (/adf00_path *= *(.*)/)
  {
  	$adas_path=$1;
          chomp($adas_path);
  }
  if (/chunk_size *= *(.*)/)
  {
  	$chunk_size=$1;
          chomp($chunk_size);
  }
  if (/chunk_size_high_ex *= *(.*)/)
  {
  	$chunk_size_high_ex=$1;
          chomp($chunk_size_high_ex);
  }
  if (/chunk_size_nex *= *(.*)/)
  {
  	$chunk_size_nex=$1;
          chomp($chunk_size_nex);
  }
  if (/matrix_type *= *(.*)/)
  {
  	$matrix_type=$1;
          chomp($matrix_type);
          $matrix_type=uc($matrix_type);
  }
  if (/internal_split *= *(.*)/)
  {
  	$internal_split=$1;
          chomp($internal_split);
  }
  if (/auto *([0-9]*) *(.*)/)
  {
  	$autos_ncpu=$1;
  	$autos_exec=$2;
          chomp($autos_exec);
  }
  if (/stg1_ex *([0-9]*) *(.*)/)
  {
  	$stg1ex_ncpu=$1;
  	$stg1ex_exec=$2;
          chomp($stg1ex_exec);
  }
  if (/stg2_ex *([0-9]*) *(.*)/)
  {
  	$stg2ex_ncpu=$1;
  	$stg2ex_exec=$2;
          chomp($stg2ex_exec);
  }
  if (/stg3_ex *([0-9]*) *(.*)/)
  {
  	$stg3ex_ncpu=$1;
  	$stg3ex_exec=$2;
          chomp($stg3ex_exec);
  }
  if (/stgb *([0-9]*) *(.*)/)
  {
  	$stgb_ncpu=$1;
  	$stgb_exec=$2;
          chomp($stgb_exec);
  }

  if (/stg1_dip *([0-9]*) *(.*)/)
  {
  	$stg1dip_ncpu=$1;
  	$stg1dip_exec=$2;
   chomp($stg1dip_exec);
  }
  if (/stg2_dip *([0-9]*) *(.*)/)
  {
  	$stg2dip_ncpu=$1;
  	$stg2dip_exec=$2;
   chomp($stg2dip_exec);
  }
  if (/stg3_dip *([0-9]*) *(.*)/)
  {
  	$stg3dip_ncpu=$1;
  	$stg3dip_exec=$2;
   chomp($stg3dip_exec);
  }
  if (/stgd_dip *([0-9]*) *(.*)/)
  {
  	$stgddip_ncpu=$1;
  	$stgddip_exec=$2;
   chomp($stgddip_exec);
  }

  if (/stg1_tcc *([0-9]*) *(.*)/)
  {
  	$stg1tcc_ncpu=$1;
  	$stg1tcc_exec=$2;
          chomp($stg1tcc_exec);
  }
  if (/stg2_tcc *([0-9]*) *(.*)/)
  {
  	$stg2tcc_ncpu=$1;
  	$stg2tcc_exec=$2;
          chomp($stg2tcc_exec);
  }
  if (/stgjk_tcc *([0-9]*) *(.*)/)
  {
  	$stgjktcc_ncpu=$1;
  	$stgjktcc_exec=$2;
          chomp($stgjktcc_exec);
  }
  if (/stg1_nx *([0-9]*) *(.*)/)
  {
  	$stg1nx_ncpu=$1;
  	$stg1nx_exec=$2;
          chomp($stg1nx_exec);
  }
  if (/stg2_nx *([0-9]*) *(.*)/)
  {
  	$stg2nx_ncpu=$1;
  	$stg2nx_exec=$2;
          chomp($stg2nx_exec);
  }
  if (/stg3_nx *([0-9]*) *(.*)/)
  {
  	$stg3nx_ncpu=$1;
  	$stg3nx_exec=$2;
          chomp($stg3nx_exec);
  }
  if (/stgf_nx *([0-9]*) *(.*)/)
  {
  	$stgfnx_ncpu=$1;
  	$stgfnx_exec=$2;
          chomp($stgfnx_exec);
  }
  if (/stgicf_nx *([0-9]*) *(.*)/)
  {
  	$stgicfnx_ncpu=$1;
  	$stgicfnx_exec=$2;
          chomp($stgicfnx_exec);
  }
  if (/stgf_ex *([0-9]*) *(.*)/)
  {
  	$stgfex_ncpu=$1;
  	$stgfex_exec=$2;
          chomp($stgfex_exec);
  }
  if (/stgicf_ex *([0-9]*) *(.*)/)
  {
  	$stgicfex_ncpu=$1;
  	$stgicfex_exec=$2;
          chomp($stgicfex_exec);
  }
  if (/stgfdamp_ex *([0-9]*) *(.*)/)
  {
  	$stgfdampex_ncpu=$1;
  	$stgfdampex_exec=$2;
          chomp($stgfdampex_exec);
  }
  if (/stgicfdamp_ex *([0-9]*) *(.*)/)
  {
  	$stgicfdampex_ncpu=$1;
  	$stgicfdampex_exec=$2;
          chomp($stgicfdampex_exec);
  }

  if (/omgmrgp *([0-9]*) *(.*)/)
  {
  	$omgmrgp_ncpu=$1;
  	$omgmrgp_exec=$2;
          chomp($omgmrgp_exec);
  }
  if (/omgmrg\s+([0-9]+)\s+(.*)/)
  {
  	$omgmrg_ncpu=$1;
  	$omgmrg_exec=$2;
          chomp($omgmrg_exec);
  }
  if (/omadd *([0-9]*) *(.*)/)
  {
  	$omadd_ncpu=$1;
  	$omadd_exec=$2;
          chomp($omadd_exec);
  }
  if (/om2omu *([0-9]*) *(.*)/)
  {
  	$om2omu_ncpu=$1;
  	$om2omu_exec=$2;
          chomp($om2omu_exec);
  }
  if (/omr2omc *([0-9]*) *(.*)/)
  {
  	$omr2omc_ncpu=$1;
  	$omr2omc_exec=$2;
          chomp($omr2omc_exec);
  }

  if (/dfilter *([0-9]*) *(.*)/)

  {
        $dfilter_ncpu=$1;
        $dfilter_exec=$2;
          chomp($dfilter_exec);
  }

  if (/arrange *([0-9]*) *(.*)/)

  {
        $arrange_ncpu=$1;
        $arrange_exec=$2;
          chomp($arrange_exec);
  }

  if (/adasexj *([0-9]*) *(.*)/)
  {
  	$adasexj_ncpu=$1;
  	$adasexj_exec=$2;
          chomp($adasexj_exec);
  }
  
  #etc
 }


	close(FD);


}

sub get_ip
{		
	open(FD,"$adas_path/$sym[$nzion-1].dat") or die ("Can't find adf00 file: $adas_path/$sym[$nzion-1].dat");
	for ($i=0;$i<$nzion-$nelec+1;$i++) { $junk=<FD>; }
	$line=<FD>;
	$line =~ / *[0-9]* *(.*) 1s.*/i;
	$ionpot = $1;
	$ionpot =~ s/d/e/i;
	
	close(FD);

        $maxe=$maxe_ionpot*$ionpot/13.606;

        $maxep=sprintf("%i",$maxe+1);
	$maxe_z2=$maxe/$charge/$charge;
}

sub get_prelim_autos
{
        
        chdir("str");
        
        open(FD,"TERMS");
        $junk=<FD>;
        @TERMSFILE=<FD>;
        close(FD);
        open(FD,"LEVELS");
        $junk=<FD>;
        @LEVELSFILE=<FD>;
        close(FD);
        $tmp1=pop(@LEVELSFILE);
        $tmp2=pop(@TERMSFILE);

	$maxla=0; $maxwa=0; $minwa=999;
	foreach(@TERMSFILE)
        {
        	/ *([0-9]*) *([0-9]*).*/;
                if ($1 > $maxwa) { $maxwa = $1; }
                if ($1 < $minwa) { $minwa = $1; }
                if ($2 > $maxla) { $maxla = $2; }
        }
#	print "@TERMSFILE \n;
#	print "@TERMSFILE $minwa $maxwa $maxla\n";
	
        # ADW - check with Mike
        $maxlt= int(($j2max_ex / 2.0) + ($maxwa)/2.0);
        $minlt_nx=int(($j2max_ex / 2.0) - ($maxwa)/2.0)+1;
        $maxlt_nx=int(($j2max_nx / 2.0) + ($maxwa)/2.0)+1;
	$minst=1+2*((($minwa-1)/2.0) - 0.5);
        $maxst=1+2*((($maxwa-1)/2.0) + 0.5);
        if ($minst==0) { $minst=2; }

        $jcutex=$j2max_ex - ($j2max_ex - $maxwa) % 2;
        $jcutnx=$j2max_nx - ($j2max_ex - $maxwa) % 2;
        $jmntwo=$j2max_ex + 2 - ($j2max_ex - $maxwa) % 2;
# maybe, check J ranges carefully

	open(FD,"olg");
	while(<FD>)
	{
		if (/EK-E1/)
		{
			$_=<FD>;
			$_=<FD>;
			/([0-9]*\.[0-9]*)$/;
			$first_thres=$1;
			while($_=<FD> and !/LIST/)
			{
				if (!/LIST/)
				{
					/([0-9]*\.[0-9]*)$/;
					if ($1>$last_thres)
                                          {$last_thres=$1;
#                                        print "READING THRESHOLDS $first_thres $last_thres \n";
					}
				}
	
			}
		
		}
	}
	close(FD);
	
	$first_thres_z2=$first_thres/$charge/$charge;
	$last_thres_z2=$last_thres/$charge/$charge;
	
	print "$first_thres $last_thres $first_thres_z2 $last_thres_z2\n";
	
	chdir("..");

}

sub del_prelim_autos
{


}

sub write_inputs_tot
{


}

sub run_code_born
{
	chdir("born");
	system("$autos_exec < das");
	$_=`tail -1 olg`;
	if (!/CPU/)
	{
		die ("Born/Radiative rate calculation didn't finish");
	}
	unlink("ols");
	unlink("oic");
	unlink("olsu");
	unlink("oicu");
	unlink("opic");
	unlink("opicu");
	unlink("opls");
	unlink("oplsu");
	unlink("SHFTIC");
	unlink("SHFTLS");
#	unlink("adasexj.in.form");
#	unlink("adasex.in.form");
	unlink("hffcin");
	unlink("RESTART");
	unlink("TERMS");
#	unlink("LEVELS");
	unlink("OVRLAP");
	unlink("radout");
	unlink("radwin");
	chdir("..");
}


sub run_code_tcc
{
	chdir("tcc");
	system("$stg1tcc_exec < dstg1");
	$_=`tail -1 rout1r`;
	if (!/CPU/)
	{
		die ("TCC stg1 calculation didn't finish");
	}

	system("$stg2tcc_exec < dstg2");
	$_=`tail -1 rout2r`;
	if (!/CPU/)
	{
		die ("TCC stg2 calculation didn't finish");
	}


	system("$stgjktcc_exec < dstgjk");
	$_=`tail -1 routjk`;
	if (!/CPU/)
	{
		die ("TCC stgjk calculation didn't finish");
	}


	chdir("..");
}

sub run_code_inner
{
	chdir("inner");

        system("echo Version 1.7 of iscript ");

	system("$stg1ex_exec < dstg1");
	
	$_=`tail -1 rout1r`;
	
	if (!/CPU/)
	{
		die ("Inner region stg1 calculation didn't finish");
	}
	system("echo $stg1ex_exec FINISHED");

	system("$stg2ex_exec < dstg2");
	
	system("echo $stg2ex_exec FINISHED");	

        if($accelerated == 1)
        {
        system("../../bin/hfilter.x");  
        }

        open(FD,"sizeH.dat") or die("unable to open sizeH.dat\n");
        @dat=<FD>;
        close(FD);

        $totalexsym == 0;

        foreach(@dat)
         {
         $totalexsym=$totalexsym + 1;
         }

        system("echo $totalexsym syms for inner");
        system("echo acceleration = $accelerated");

        if($accelerated == 1)
        {
        $np_per_subworld=$stg3ex_ncpu/$totalexsym;
        $nprow=sqrt($np_per_subworld);
        $npcol=sqrt($np_per_subworld);
        }
        else
        {    

        $sqrt_cpu=sqrt($stg3ex_ncpu);
        $rem = $sqrt_cpu % 1;

        if (int($sqrt_cpu)==$sqrt_cpu)
        {
                $nprow=$sqrt_cpu;
                $npcol=$sqrt_cpu;
        }
        else
        {
                $nprow=int($sqrt_cpu);
                while(int($stg3ex_ncpu / $nprow)*$nprow != $stg3ex_ncpu)
                {
                        $nprow--;
                }
                $npcol=$stg3ex_ncpu / $nprow
        }
        }

        if($npcol * $nprow * $totalexsym == $stg3ex_ncpu)
        {
        $totalexsym = 1;
        } 

        open(FD,"> dstg3") or die("Can't open  dstg3");
        print FD "S.S. Automatically generated ADAS8#3 stg3 input\n";
        print FD " &prediag npw_per_subworld=$totalexsym &END\n ";
        print FD " &STG3A$radeq &END\n";
        print FD " &STG3B INAST=0 NAST=0 &END\n";
#        if ($accelerated == 1)
#        {
        print FD " &MATRIXDAT NPROW=$nprow NPCOL=$npcol NDIV=1 &END\n";
#        }
        close(FD);

	system("$stg3ex_exec ");
	
        if ($accelerated == 1 || $stg3ex_ncpu > 1 )
        {
                system("cat H.DAT[0-9][0-9][0-9] > H.DAT ");
        }
	system("echo $stg3ex_exec FINISHED");	
        
	chdir("..");
}

sub run_code_dip
{

        if($rad_damp == 0) { return; }

	chdir("dip");
	system("$stg1dip_exec < dstg1");

	$_=`tail -1 rout1r`;
	if (!/CPU/)
	{
		die ("Dipole region stg1 calculation didn't finish");
	}
	
	system("echo $stg1dip_exec  FINISHED ");	

	system("$stg2dip_exec < dstg2");
	
	system("echo $stg2dip_exec  FINISHED ");

        if ($stg2dip_ncpu > 1){
	
        system("$dfilter_exec");

        system("echo dfilter  FINISHED ");

                              }
	
	system("$stg3dip_exec < dstg3");
	
	system("echo $stg3dip_exec  FINISHED ");	
	
	$_=`tail -1 rout3r`;
	if (!/CPU/)
	{
		die ("Dipole region stg3 (maybe stg2) calculation didn't finish");

	}
 if ($stg3dip_ncpu > 1)
 {
        open(FD,"sizeH.dat") or die("unable to open sizeH.dat\n");
        @dat=<FD>;
        close(FD);

        $npairsym == 0;

        foreach(@dat)
         {
         $npairsym=$npairsym + 1;
         }

        $npairsym = $npairsym/2;

        system("echo $npairsym syms for pstgd");
        system("echo acceleration = $accelerated");


        open(FD,"> dstgd") or die("Can't open  dstgd");
        print FD " &prediag nproc_per_subworld=4 ndip=1 ndipole=$npairsym &END\n ";
        close(FD);

        $stgddip_ncpu=$npairsym*4;

        system("echo $stgddip_exec");
        system("$stgddip_exec");

   system("echo $stgddip_exec  FINISHED ");
   
 }
        
	chdir("..");

}

sub run_code_stgb
{
  if($rad_damp == 0) { return; }
 	if (!-e "outer/dstgb") { write_inputs_stgb(); }

  chdir("outer");
  symlink("../inner/H.DAT","H.DAT");
  system("$stgb_exec < dstgb");

  if($stgb_ncpu == 1){
  $_=`tail -1 routb`;
 	if (!/CPU/)
  {
    die ("stgb calculation didn't finish");
 	}
  }

   if($stgb_ncpu > 1){
    if($totalboundsym  >  $stgb_ncpu){
    die ("pstgb calculation input.dat mismatch"); 
   }  
  }

   system("echo $stgb_exec  FINISHED ");
	 chdir("..");
}


sub run_code_outer
{
	chdir("outer");
#	@list = <*>;
        @list = <*>;
	foreach (@list)
	{
		next if ($outer_region_dir and $outer_region_dir!=$_);
		if (-d)
		{
			chdir($_);
#                        $_=`tail -1 routicfdamp`;
#                        if (/CPU/)
#                        {
#                system("echo 'seems complete... exiting dir'");
#                        goto label1
#                        }
#
#
#             For auger damping ... file auger.inp   in AUGER directory
#
                        if($auger_damp > 0){
                        symlink("../../auger/AWIDTH","AUGER");
                                            }
#
#
			$dir=$_;
                        $label=$dir-1;
#
                        $requiremerge=0;
                        
                        if ($rad_damp==1)
                        {
                                symlink("dstgf","dstgfdamp");
                                symlink("dstgicf","dstgicfdamp");
                                
                                chdir("../../dip");                                
                                @dbfiles=<[D][0-9][0-9]>;
                                chdir("../outer/$dir");
                                foreach (@dbfiles)
                                {
                                        symlink("../../dip/$_","$_");
                                }

                                chdir("../../outer");                                
                                @dbfiles=<[B]*>;
                                chdir("../outer/$dir");
                                foreach (@dbfiles)
                                {
                                        symlink("../$_","$_");
                                }
                                                                
        	        	system("$stgfdampex_exec < dstgf");
	        		$_=`tail -1 routfdamp`;
        			if (!/CPU/)
	        		{
        				die ("Outer region stgfdamp (directory $dir) calculation didn't finish");
	        		}

                		system("$stgicfdampex_exec < dstgicf");

	        		$_=`tail -1 routicfdamp`;
        			if (!/CPU/)
		        	{
	        			die ("Outer region stgicfdamp (directory $dir) calculation didn't finish");
        			}

                                foreach (@dbfiles)
                                {
                                        unlink();
                                }
                                

		        	if ($stgicfdampex_ncpu > 1)
                                {
                                	$requiremerge=1;
                	        }
                        }
                        else
                        {
        	        	system("$stgfex_exec < dstgf");
	        		$_=`tail -1 routf`;
        			if (!/CPU/)
	        		{
        				die ("Outer region stgf (directory $dir) calculation didn't finish");
	        		}

                		system("$stgicfex_exec < dstgicf");

	        		$_=`tail -1 routicf`;
        			if (!/CPU/)
		        	{
	        			die ("Outer region stgicf (directory $dir) calculation didn't finish");
        			}

		        	if ($stgicfex_ncpu > 1)
                                {
                                	$requiremerge=1;
                	        }
                        }

			unlink("adasexj.in.formg");
			unlink("adasexj.in.form");
			unlink("strength.dat");
			unlink("term.dat");
                        unlink("AUGER");

                        if ($requiremerge==1)
                        {
                                system("$arrange_exec");
#  			        system("rm -f OMEGA[0-9]*");
 			        system("rm -f kmtls* jbinls* zkmls*");
				if ($label < 10) 
				{
                                system("mv OMEGAZ ../OMEGA00$label");
				}
				elsif($label < 100)
				{
				system("mv OMEGAZ ../OMEGA0$label");
				}
				else
				{
				system("mv OMEGAZ ../OMEGA$label");
				}
                        }
			else
			{
 			        system("rm -f kmtls* jbinls* zkmls*");
				if ($label < 10) 
				{
                                system("mv OMEGA ../OMEGA00$label");
				}
				elsif($label < 100)
				{
				system("mv OMEGA ../OMEGA0$label");
				}
				else
				{
				system("mv OMEGA ../OMEGA$label");
				}
			}
                        label1:
                        chdir("..");
		}
	}
#
#       will now merge the OMEGAXXX from each directory 1-9
#
        system("$arrange_exec");
        write_omegafilter();
 	chdir("..");
       

}

sub run_code_nonxinn
{
        
	chdir("nonx") or die ("Can't change directory to nonx");
	system("cp ../str/radout radial");
	system("$stg1nx_exec ");
	$_=`tail -1 rout1r`;
	if (!/CPU/)
	{
		die ("Non-exchange stage1 calculation didn't finish");
	}
	system("echo $stg1nx_exec  FINISHED ");

	system("$stg2nx_exec ");
	
	system("echo $stg2nx_exec  FINISHED ");	

        if($accelerated == 1)
        {
        system("../bin/hfilter.x");  
        }

        open(FD,"sizeH.dat") or die("unable to open sizeH.dat\n");
        @dat=<FD>;
        close(FD);

        $totalnexsym == 0;

        foreach(@dat)
         {
         $totalnexsym=$totalnexsym + 1;
         }

        system("echo $totalnexsym syms for nonex inner");
        system("echo acceleration = $accelerated");

        if($accelerated == 1)
        {
        $np_per_subworld=$stg3nx_ncpu/$totalnexsym;
        $nprow=sqrt($np_per_subworld);
        $npcol=sqrt($np_per_subworld);
        }
        else
        {    

        $sqrt_cpu=sqrt($stg3nx_ncpu);
        $rem = $sqrt_cpu % 1;

        if (int($sqrt_cpu)==$sqrt_cpu)
        {
                $nprow=$sqrt_cpu;
                $npcol=$sqrt_cpu;
        }
        else
        {
                $nprow=int($sqrt_cpu);
                while(int($stg3nx_ncpu / $nprow)*$nprow != $stg3nx_ncpu)
                {
                        $nprow--;
                }
                $npcol=$stg3nx_ncpu / $nprow
        }
        }

        if($npcol * $nprow * $totalnexsym == $stg3nx_ncpu)
        {
        $totalnexsym = 1;
        } 

        open(FD,"> dstg3") or die("Can't open  dstg3");
        print FD "S.S. Automatically generated ADAS8#3 stg3 input\n";
        print FD " &prediag npw_per_subworld=$totalnexsym &END\n ";
        print FD " &STG3A$radeq &END\n";
        print FD " &STG3B INAST=0 NAST=0 &END\n";
        print FD " &MATRIXDAT NPROW=$nprow NPCOL=$npcol &END\n";
        close(FD);

	system("$stg3nx_exec ");
	
        if ($accelerated == 1 || $stg3nx_ncpu > 1)
        {
           system("cat H.DAT[0-9][0-9][0-9] > H.DAT ");
        }

	system("echo $stg3nx_exec FINISHED");
	
	chdir("..");
}

sub run_code_nonxout
{
	chdir("nonx") or die ("Can't change directory to nonx");

        if($stgfnx_ncpu == 1)
	
	{


	system("$stgfnx_exec < dstgf");
	$_=`tail -1 routf`;
	if (!/CPU/)
	{
		die ("Non-exchange stgf calculation didn't finish");
	}

	system("$stgicfnx_exec < dstgicf");
	$_=`tail -1 routicf`;
	if (!/CPU/)
	{
		die ("Non-exchange stgicf calculation didn't finish");
	}
	
	}
	else
	{
		system("$stgfnx_exec");
	$_=`tail -1 routf`;	
	if (!/CPU/)
	{
		die ("Non-exchange pstgf calculation didn't finish");
	}		
		system("$stgicfnx_exec");
		
	$_=`tail -1 routicf`;	
	if (!/CPU/)
	{
		die ("Non-exchange pstgicf calculation didn't finish");	
	}
        }

        if ($stgicfnx_ncpu == 1)	
	{
	system("mv OMEGA OMEGA000");
	}
	system("$arrange_exec");
        write_omegafilternon();
#	system("rm -f OMEGA[0-9]*");
	system("rm -f kmtls* jbinls*");
	chdir("..");
}

sub process_arguments
{

        use Getopt::Long;
        use Cwd;


	$actionhelp=0;
        $actiondir=0;	# Make directories
        $actionpre=0;
        $actionaug=0;
        $actioninp=0;
        $actiontcc=0;
        $actioninn=0;
        $actionout=0;
        $actionnoi=0;
        $actionnoo=0;
        $actionadf=0;
        $actionmrg=0;
        $actionbrn=0;
        $actioncln=0;
        $actionvcl=0;
        $actiondel=0;
        $actionrep=0;
        $actionxmp=0;
	$rerun_str=0;
        $nooption=1;

	$initialdir=getcwd();


 if ( @ARGV > 0 ) {
        
        $ok=GetOptions(\%options,       "proc=s",
                                        "root=s",
					"help",
                                        "report",
                                        "example",
                                        "dir",
                                        "veryclean",
                                        "clean",
                                        "archive",
                                        "delete",
                                        "all",
                                        "inp",
                                        "run",
                                        "inner",
                                        "dip",
                                        "tcc",
                                        "stgbinp",
                                        "stgb",
                                        "outer:i",
                                        "noninn",
                                        "nonout",
                                        "born",
                                        "merge",
                                        "adf04"
                                        );
    }

	($name, $pass, $uid, $gid, $quota, $name, $gcos, $dir, $shell)  = getpwuid($<);

        if ($options{proc})
        {
                $procfile=$options{proc};
        }
        else
        {
        	$procfile="$dir/.adas803proc";

        }
        if ($options{root})
        {
		$options{root} =~ s/PID/$</;
		if (!-d $options{root}){ mkdir($options{root}, 0755) || die "Cannot make directory $options{root}"; }
		chdir($options{root}) || die("Can't change directory to adas");
	}
	if ($options{"help"})
	{
		$actionhelp=1;
		return;
        }
	if ($options{"report"})
	{
		$actionrep=1;
		return;
        }
	if ($options{"example"})
	{
		$actionxmp=1;
		return;
        }
	if ($options{"veryclean"})
	{
		$actionvcl=1;
		return;
        }
	if ($options{"archive"})
	{
		$actionarc=1;
        }
	if ($options{"clean"})
	{
		$actioncln=1;
		return;
        }
	if ($options{"delete"})
	{
		$actiondel=1;
		return;
        }

        $count=@ARGV;

        $nzion=$ARGV[1];
        
        if ($count < 1)
        {

                $actionhelp=1;
                return;        
        }
        
        $physicsfile=$ARGV[0];
        
	$physicsfile = $initialdir . "/" . $physicsfile if ($options{root} && !($physicsfile =~ /^\//));
        	
        $defaultz=1;
        
        if ($count > 1)
        {
                $nzion=$ARGV[1];
                
		
		if ($nzion =~ /[A-Z,a-z]/)
		{
			$nzion_c=$nzion;
			for ($i=0;$i<@sym;$i++)
			{
				if ($nzion_c eq $sym[$i])
				{
					$nzion=$i+1;
				}
			}
		}
		
		$defaultz=0;
        }
        else
        {
                if (-e "str/das")
                {
                       	open (FD,"str/das") or die ("Can't open str/das");
                        while(<FD>)
                        {
                                if (/NZION=([0-9]*)/)
                                {
                                        $nzion=$1;
                                        $defaultz=0;
                                }
                        }
                
                }
        }


        if ($defaultz)
        {
                print "Ion charge not specified, using 26\n";
                $nzion=26;
        }
        if (!-e $physicsfile)
        {
                $actionhelp=1;
                return;
        }

	if ($options{"dir"})
	{
		$actiondir=1;
                $nooption=0;
        }

	if ($options{"inp"})
	{
	        $actiondir=1;
	        $actionpre=1;
	        $actioninp=1;
                $nooption=0;
		$rerun_str=1;
	}
	if ($options{"run"})
	{
        	$actiontcc=1;
	        $actioninn=1;
                $actiondip=1;
                $actionstgbinp=1;
                $actionstgb=1;
	        $actionout=1;
	        $actionnoi=1;
	        $actionnoo=1;
	        $actionadf=1;
	        $actionmrg=1;
	        $actionbrn=1;
                $nooption=0;
	}
	if ($options{"inner"})
	{
	        $actionpre=1;
	        $actioninn=1;
                $nooption=0;
	}
	if ($options{"dip"})
	{
	        $actionpre=1;
	        $actiondip=1;
                $nooption=0;
	}
	if ($options{"tcc"})
	{
	        $actionpre=1;
	        $actiontcc=1;
                $nooption=0;
	}
	if ($options{"stgbinp"})
	{
	        $actionpre=1;
	        $actionstgbinp=1;
                $nooption=0;
	}
	if ($options{"stgb"})
	{
	        $actionpre=1;
	        $actionstgb=1;
                $nooption=0;
	}
	if (defined($options{"outer"}))
	{
		$outer_region_dir=$options{"outer"};
	        $actionpre=1;
	        $actionout=1;
                $nooption=0;
	}
	if ($options{"noninn"})
	{
	        $actionpre=1;
	        $actionnoi=1;
                $nooption=0;
	}
	if ($options{"nonout"})
	{
	        $actionpre=1;
	        $actionnoo=1;
                $nooption=0;
	}
	if ($options{"born"})
	{
	        $actionpre=1;
	        $actionbrn=1;
                $nooption=0;
	}
	if ($options{"merge"})
	{
	        $actionpre=1;
	        $actionmrg=1;
                $nooption=0;
	}
	if ($options{"adf04"})
	{
	        $actionpre=1;
	        $actionadf=1;
                $nooption=0;
	}
	if ($options{"all"} or $nooption)
	{
	        $actiondir=1;
	        $actionpre=1;
	        $actioninp=1;
	        $actioninn=1;
                $actiondip=1;
                $actionstgbinp=1;
                $actionstgb=1;
        	$actiontcc=1;
	        $actionout=1;
	        $actionnoi=1;
	        $actionnoo=1;
	        $actionadf=1;
	        $actionmrg=1;
	        $actionbrn=1;
	}
      
}


sub clean_up_very
{
	clean_up();
	@files=`find . -name 'rout*'`;
	foreach (@files)
	{
		chomp();
		unlink();
	}	
	@files=`find ./outer -name 'OMEGA*'`;
	foreach (@files)
	{
		chomp();
		unlink();
	}
        @files=`find ./outer -name 'B*'`;
        foreach (@files)
        {
                chomp();
                unlink();
        }

        if(rad_damp == 1)
        {

        @files=`find ./dip -name 'D*'`;
        foreach (@files)
        {
                chomp();
                unlink();
        }
        @files=`find ./dip -name 'symm*'`;
        foreach (@files)
        {
                chomp();
                unlink();
        }
        @files=`find ./dip -name 'RK*'`;
        foreach (@files)
        {
                chomp();
                unlink();
        }
        @files=`find ./dip -name 'STG2*'`;
        foreach (@files)
        {
                chomp();
                unlink();
        }

        }

	@files=`find ./nonx -name 'OMEGA*'`;
	foreach (@files)
	{
		chomp();
		unlink();
	}	
	unlink("inner/H.DAT");
	unlink("nonx/H.DAT");
	unlink("inner/NX1.DAT");
	unlink("inner/NX2.DAT");
	unlink("tcc/TCCDW.DAT");
	unlink("born/adf04ic");
	unlink("born/adf04ls");
#	unlink("born/olg");
#	unlink("str/olg");
	unlink("str/adasexj.in.form");
	unlink("str/adasex.in.form");
	unlink("adas/adasexj.in");
	unlink("adas/adasexj.out");




}

sub clean_up
{
	@files=`find . -name 'fort.*'`;
	foreach (@files)
	{
		chomp();
		unlink();
	}	
	unlink("inner/STG1.DAT");
	unlink("inner/STG1g.DAT");
	unlink("tcc/STG1.DAT");
	unlink("tcc/TCC.DAT");
	unlink("tcc/RECUPH.DAT");
	unlink("tcc/RK.DAT");
	unlink("tcc/AMAT.DAT");
	unlink("inner/sizeH.dat");
	unlink("nonx/sizeNX.dat");
	unlink("nonx/term.dat");
	unlink("nonx/ANG1.DAT");
	unlink("nonx/ANG2.DAT");
	unlink("nonx/ANG3.DAT");
	unlink("adas/adf04");
	@files=<inner/STG2H*>;
	foreach (@files)
	{
		unlink();
	}
	@files=<tcc/STG2H*>;
	foreach (@files)
	{
		unlink();
	}
	@files=<inner/RK*>;
	foreach (@files)
	{
		unlink();
	}
	@files=<inner/time*>;
	foreach (@files)
	{
		unlink();
	}
	@files=<nonx/time*>;
	foreach (@files)
	{
		unlink();
	}
	@files=<inner/NX?g.DAT>;
	foreach (@files)
	{
		unlink();
	}
	@files=<nonx/RAD*>;
	foreach (@files)
	{
		unlink();
	}

}


sub display_help
{
	print "          ADAS8#3 - Automated R-matrix calculations\n";
	print "          -----------------------------------------\n";	
	print "Usage: $0 [options] input.dat [Z]\n";
	print "\n";	
	print " input.dat     - ADF41 file to control calculation.\n";	
	print " Z             - Nuclear charge (or element symbol) for calculation.\n";	
	print "\n";	
	print "                 Valid Options\n";	
	print "                 +++++++++++++\n";	
	print " --help        - display this message\n";
	print " --example     - print out example input file\n";
	print " --report      - analyse directories for calculation progress\n";	
	print " --dir         - create directory structure\n";	
	print " --clean       - remove large (unnecssary) passing files\n";	
	print " --veryclean   - leave only inputs, collision strenths and adf04\n";	
	print " --delete      - delete all subdirectories\n";
	print " --archive     - copy files to parent directory\n";
	print " --inp         - generate input files (implies dir)\n";	
	print " --inner       - do only inner region (exchange)\n";	
	print " --tcc         - create TCCDW.DAT file\n";	
	print " --dip         - do only inner region (dipole, if damped)\n";	
	print " --stgbinp     - make input file for STGB (if inner finished)\n";	
	print " --stgb        - run STGB\n";	
	print " --outer       - do only outer region (exchange)\n";	
	print " --noninn      - do inner region non-exchange calculation\n";	
	print " --nonout      - do outer region non-exchange calculation\n";	
	print " --born        - calculate born limits and non-dipole A-values\n";	
	print " --merge       - merge collision strengths\n";	
	print " --adf04       - generate adf04 file\n";	
	print " --run         - run whole calculation (from inner to adf04)\n";	
	print " --all         - [default] '--inp' then '--run'\n";	
	print " --proc=file   - alternative excecutables file, defaults to ~/.adas803proc\n";	
	print " --root=path   - top level directory to run in, defaults to current\n";	

}

sub status_report
{


}

sub write_inputs_adas
{
	open (FD,"born/adasexj.in.form") or die ("Can't open str/adasexj.in.form");
	open (FO,"> adas/adasexj.in") or die ("Can't open adas/adasexj.in");
        
        $iel=uc($sym[$nzion-1]);
        $ionpot_cm=$ionpot * 8065.54445;
        
        @ADASEXJIN=<FD>;
        splice(@ADASEXJIN,0,1); 
#        $_=<FD>;

	@t_grid=(2.0e+02,5.0e+02,1.0e+03,2.0e+03,5.0e+03,1.0e+04,2.0e+04,5.0e+04,
		 1.0e+05,2.0e+05,5.0e+05,1.0e+06,2.0e+06);
 
        open(FD,"born/LEVELS");
        $junk=<FD>;
        @LEVELSFILE=<FD>;
        close(FD);
        $levels=@LEVELSFILE-1;

	$numtmp=@t_grid;
#        s/&END/ NLEVS=$levels FIPOT=$ionpot_cm IEL='$iel' NUMTMP=$numtmp IRDTMP=1 MXTMP=1 &END/;
        print FO "&ADASEX NLEVS=$levels FIPOT=$ionpot_cm IEL='$iel' NUMTMP=$numtmp IRDTMP=1 MXTMP=1 &END\n";
	
	$charge=$nzion-$nelec;
        if($charge == 0)
        {
        $charge =1;
        }

	foreach (@t_grid)
	{
		printf FO "%7.2e ",$_ * ($charge+1)*($charge+1);
	}
	print FO "\n";

        while (<FD>)
        {
        	print FO;
        }
        print FO @ADASEXJIN;
	close(FD);
	close(FO);
}

sub make_omega_file
{
	
        
        chdir("adas") || die("Can't change directory to adas");
        
        symlink ("../outer/OMEGAZ-YMULT1000","../adas/omadd1");
        
        symlink("../nonx/OMEGAZ-NON-YMULT5","omadd2");
        system("$omadd_exec");
	$_=`tail -1 oadd`;
	if (!/CPU/)
	{
		die ("Addition of non-exchange cross-sections failed");
	}
	rename("omaddt","OMEGA000");
        
        symlink("../born/OMGINFIC","OMEGA001");
        system("$arrange_exec");

#	unlink("omadd1");
#	unlink("omadd2");
#	unlink("OMEGA001");
	rename("OMEGAZ","omega");
        
	chdir("..");
}

sub make_adf_file
{
	chdir("adas") or die("Can't change directory to adas");
                
        system("$adasexj_exec");
	
        open(FD,"../born/adf04ic") or die("Can't open file ../born/adf04ic");
        
        while(<FD>)
        {
                
                if (/^ +([0-9]*) +([0-9]*) ([0-9]\.[0-9][0-9][\+,-][0-9][0-9]).*/)
        	{
                        
                        $acoeff[$1][$2]=$3;
                
                }
        
        }
           
        ($name, $pass, $uid, $gid, $quota, $name, $gcos, $dir, $shell)  = getpwuid($<);
                        
        if ($gcos =~ /Allan.*White/)
        {
        	$initials="adw";
        }
        if ($gcos =~ /Mike.*Witt/)
        {
        	$initials="mcw";
        }
        if ($gcos =~ /Nigel.*Badn/)
        {
        	$initials="nrb";
        }
        if ($gcos =~ /Hugh.*Sum/)
        {
		$initials="hps";
        }
        if ($gcos =~ /Stuart.*Loch/)
        {
		$initials="sdl";
        }
        if ($gcos =~ /Connor.*Ballance/)
        {
		$initials="cpb";
        }
        if ($gcos =~ /Di.*Wu/)
        {
		$initials="diw";
        }
        if ($gcos =~ /Rebecca.*Rogers/)
        {
		$initials="rar";
        }
        if ($gcos =~ /Jonathan.*Pearce/)
        {
		$initials="jap";
        }
        if ($gcos =~ /Nathalia.*Alzate/)
        {
		$initials="nat";
        }
        
        
        if (!$initials)
        {
        	$initials = $gcos;
                $initials =~ s/^([A-Z,a-z]).* ([A-Z,a-z]).*/\1\2/; 
        	$initials = lc($initials);
        }
        
        
        $sequence=$sym[$nelec-1];
        
        ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

	if ($year > 100) { $year -= 100 };

	$charge=$nzion-$nelec;
        if($charge == 0)
        {
         $charge ==1;
        }
        $filename=$sequence . "like_" . $initials . sprintf("%02d",$year) . "#" . lc($sym[$nzion-1]) . $charge . ".dat";

	$configlen=18;
	open(FD,"adf04") or die ("Can't open adf04 file");
	$_=<FD>;
	while( <FD> )
	{
		if (s/(.*)\([0-9][A-Z]\)(.*) /\1\2/)
		{
			/ *([0-9]*)  *([0-9].*) \((.*)\)(.*)\((.*)\) *(.*$)/;
			$config=$2;
			while ($config =~ s/([0-9])([a-z])([0-9])([a-z])/\1\2 1 \3\4/) { }
			$config =~ s/([0-9])([a-z])$/\1\2 1/;
			while ($config =~ s/([a-z]) ([0-9])/\1\2/) { }
			while ($config =~ s/([a-z])([0-9])([0-9])([a-z])/\1\2 \3\4/) { }
			$config =~ s/ *$//;
			if (length($config) > $configlen) { $configlen=length($config); }  
		}
	}
	
	close(FD);
	
	
        	
	open(FA,"> $filename") or die ("Can't open $filename");
        
        $firstpass=1;
	open(FD,"adf04") or die ("Can't open adf04 file");
        while(<FD>)
        {
        	if (!/^C/ or $firstpass)
      		{
                	$printed=0;
                        if (s/(.*)\([0-9][A-Z]\)(.*) /\1\2/)
                        {
					
                                / *([0-9]*)  *([0-9].*) \((.*)\)(.*)\((.*)\) *(.*$)/;
                                
				
                                $idx=$1;
                                $config=$2;
                                $mult=$3;
                                $lquantum=$4;
                                $jquantum=$5;
                                $energy=$6;
			
				
                                while ($config =~ s/([0-9])([a-z])([0-9])([a-z])/\1\2 1 \3\4/) { }
                                $config =~ s/([0-9])([a-z])$/\1\2 1/;
                                while ($config =~ s/([a-z]) ([0-9])/\1\2/) { }
                                while ($config =~ s/([a-z])([0-9])([0-9])([a-z])/\1\2 \3\4/) { }
                                $config =~ s/ *$//;
                                
				#print "++$config++\n";
                                
				$spaces=" "x($configlen-length($config));
				
                                printf FA "%5d %s%s(%1d)%1d(%4.1f)    %12.1f\n",$idx,uc($config),$spaces,$mult,$lquantum,$jquantum,$energy;
                        	
				# Note, I was just using %-18s rather than the spaces nonsense which worked on some versions of perl but not others.
				
				$printed=1;
                        }
                        
                        # check for data line and supplement a-value
			
                        if(/^ +([0-9]*) +([0-9]*) ([0-9]\.[0-9][0-9][\+,-][0-9][0-9])(.*)/)
                        {
                        	
                                if (!$acoeff[$1][$2] or $acoeff[$1][$2] eq "0.00+00")
                                {
                                	$acoeff[$1][$2]="1.00-30";
                                }
                                printf FA "%4d%4d %s%s\n",$1,$2,$acoeff[$1][$2],$4; 
				$printed=1;
                        }

                        
                        if ($printed==0)
                        {
                		print FA;
			}
        		$firstpass=0;
                }
        }
        
	$mesh_fine_ev=$mesh_fine*$charge*$charge*13.606;
	$mesh_coarse_ev=$mesh_coarse*$charge*$charge*13.606;        
	$maxe_ev=$maxe*13.606;
	$maxe_z2_p=sprintf("%8f",$maxe_z2);
	$mesh_fine_p=sprintf("%8f",$mesh_fine);
	$mesh_coarse_p=sprintf("%8f",$mesh_coarse);
	chdir("../inner");
	get_smallest_max();
	chdir("../adas");
	$smallest_max_ev=sprintf("%8f",$smallest_max*13.606);
	$type='ICFT';
	$formatted_date=sprintf("%02d/%02d/%02d",$mday,$mon+1,$year);

        $comments  = "C------------------------------------------------------------------------\n";
        $comments .= "C\n"; 
	if ($type eq 'ICFT')
        {
	$comments .= "C        Result of ICFT R-matrix calculation\n"; 
        $comments .= "C        +++++++++++++++++++++++++++++++++++\n";
	}
	$comments .= "C\n"; 
	$comments .= "C Calculation automated by offline ADAS code ADAS8\#3\n"; 
	$comments .= "C\n"; 
	$comments .= "C     Summary Information\n"; 	
	$comments .= "C     +++++++++++++++++++\n"; 	
	$comments .= "C\n"; 
	$comments .= "C *  Ionisation energy taken from ADAS adf00 files\n"; 
	$comments .= "C\n"; 
	$comments .= "C *  Energy levels calculated by AUTOSTRUCTURE\n"; 
	$comments .= "C\n"; 
	$comments .= "C *  Radiative rates calculated by AUTOSTRUCTURE\n"; 
	$comments .= "C\n"; 
	$comments .= "C *  Effective collision strengths calculated by ICFT\n"; 
	$comments .= "C    R-matrix calculation using target structure from\n"; 
	$comments .= "C    AUTOSTRUCTURE.\n"; 
	$comments .= "C\n"; 
	$comments .= "C     Calculation details\n"; 	
	$comments .= "C     +++++++++++++++++++\n";
	$comments .= "C\n"; 	
	$comments .= "C The following scaling parameters on bound orbitals were used:\n"; 	
	$comments .= "C\n";
	$i=0;
	foreach (@orbitals)
	{
		$comments .= "C         $_ = $scale[$i++]\n"
	}
	$comments .= "C\n";
	$comments .= "C An exchange calculation was performed up to J=" . $jcutex /2 . "\n"; 	
	$comments .= "C A non-exchange calculation was then performed up to J=" . $jcutnx /2 . "\n"; 	
	$comments .= "C Followed by a top-up calculation to J -> infinity.\n"; 	
	$comments .= "C\n"; 	
	$comments .= "C $maxc continuum basis orbitals were used, giving a smallest maximum\n";
	$comments .= "C basis-orbital energy of $smallest_max Rydbergs ($smallest_max_ev eV)\n"; 	
	$comments .= "C\n"; 	
	$comments .= "C A fine energy mesh of $mesh_fine_p Z**2 Rydbergs ($mesh_fine_ev eV) was used\n";
	$comments .= "C between the first and last thresholds for the excitation calculation.\n"; 	
	$comments .= "C\n"; 	
	$comments .= "C A coarse energy mesh of $mesh_coarse_p Z**2 Rydbergs ($mesh_coarse_ev eV) was used\n";
	$comments .= "C from the last threshold up to $maxe_z2_p Z**2 Rydbergs ($maxe_ev eV) in the\n";
	$comments .= "C excitation calculation.\n"; 	
	$comments .= "C\n"; 	
	$comments .= "C A coarse energy mesh of $mesh_coarse_p Z**2 Rydbergs ($mesh_coarse_ev eV) was used\n";
	$comments .= "C from first threshold up to $maxe_z2_p Z**2 Rydbergs ($maxe_ev eV) in the\n";
	$comments .= "C non-exchange calculation.\n"; 	
	$comments .= "C\n"; 	
	$comments .= "C Dipole and Born limits were used from the last calculated energy points\n";
	$comments .= "C up to infinite energy to complete the effective collision strengths\n"; 	
#	$comments .= "C\n"; 	
#	$comments .= "C     Code versions\n";
#	$comments .= "C     +++++++++++++\n";
#	$comments .= "C\n"; 	
#	$comments .= "C  ADAS8#3       : V 0.9B\n"; 	
#	$comments .= "C  Autostructure : $autos_version\n"; 	
#	$comments .= "C  Stage 1   EX  : V 0.9B\n"; 	
#	$comments .= "C  Stage 2   EX  : V 0.9B\n"; 	
#	$comments .= "C  Stage 3   EX  : V 0.9B\n";   
#	$comments .= "C  Stage 1   NX  : V 0.9B\n";   
#	$comments .= "C  Stage 2   NX  : V 0.9B\n";   
#	$comments .= "C  Stage 3   NX  : V 0.9B\n";   
#	$comments .= "C  Stage 1   TCC : V 0.9B\n";   
#	$comments .= "C  Stage 2   TCC : V 0.9B\n";   
#	$comments .= "C  Stage jk  TCC : V 0.9B\n";   
#	$comments .= "C  Stage f   EX  : V 0.9B\n";   
#	$comments .= "C  Stage icf EX  : V 0.9B\n"; 	
#	$comments .= "C  Stage f   NX  : V 0.9B\n"; 	
#	$comments .= "C  Stage icf NX  : V 0.9B\n"; 	
#	$comments .= "C  ADASEXJ       : V 0.9B\n"; 	
	$comments .= "C\n"; 	
	$comments .= "C Produced by: $gcos\n";
        $comments .= "C Date:        $formatted_date\n";
        $comments .= "C------------------------------------------------------------------------\n";
	
        print FA $comments;
        
        close(FA);
        close(FD);
        chdir("..");

}

sub get_smallest_max
{
	if (!open(F1,"rout1r")) { $smallest_max=0; return;  }
	$smallest_max=999999999;
	$onebefore="-";
	while(<F1>)
	{
		if (/BUTTLE/)
		{
			$twobefore =~ /^ +[-,\.,0-9]* +([\.,0-9]*) +[0-9]*/;
			if ($1 < $smallest_max) { $smallest_max=$1; }
		}
		$twobefore=$onebefore;
		$onebefore=$_;
	}
	close(F1);
}

sub clean_delete
{
	system("rm -rf auger born inner adas nonx outer str tcc dip");

}

sub example_input
{
        print "# Example Input for ADAS8#3\n";
        print "# Most options have sensible defaults and\n";
        print "# don't need to be specified. The only mandatory\n";
        print "# input is a configuration list\n";
        print "# All input is flexible format\n";
	print "\n";
        print "GENERAL \n";
        print "2Jmax_ex = 23\n";
        print "2Jmax_nx = 90\n";
        print "maxc = 30\n";
        print "mesh_fine = 0.00001\n";
        print "mesh_coarse = 0.001\n";
        print "maxe/ionpot = 3\n";
#        print "YMULT      = 1000\n";
#        print "YMULT (NX) = 1.1\n";
	print "\n";
	print "CONFIGURATION LIST\n";
	print "1s2\n";
	print "1s 2s\n";
	print "1s 2p\n";
	print "1s 3s\n";
	print "1s 3p\n";
	print "1s 3d\n";
	print "\n";
	print "SCALING PARAMETERS\n";
	print "1s = 1.0\n";
	print "2s = 1.0\n";
	print "2p = 1.0\n";
	print "3s = 1.0\n";
	print "3p = 1.0\n";
	print "3d = 1.0\n";
}

sub archive_output
{
	use File::Copy;
        copy("tcc/TCCDW.DAT","$initialdir/TCCDW.DAT") if -e "tcc/TCCDW.DAT";
        copy("adas/OMEGA","$initialdir/OMEGA") if -e "adas/OMEGA";
        copy("born/olg","$initialdir/olg") if -e "born/olg";
        copy("str/das","$initialdir/adf27.dat") if -e "str/das";
        chdir("adas");
        @filelist=<*like*>;
        copy($filelist[0],"$initialdir/$filelist[0]") if (@filelist && -e $filelist[0]);
}

##################################################################################

@sym = ('h' ,'he','li','be','b' ,'c' ,'n' ,'o' ,'f' ,'ne','na',
    	'mg','al','si','p' ,'s' ,'cl','ar','k' ,'ca','sc','ti',
    	'v' ,'cr','mn','fe','co','ni','cu','zn','ga','ge','as',
    	'se','br','kr','rb','sr','y' ,'zr','nb','mo','tc','ru',
    	'rh','pd','ag','cd','in','sn','sb','te','i' ,'xe','cs',
    	'ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
    	'ho','er','tm','yb','lu','hf','ta','w', 're','os','ir',
    	'pt','au','hg','tl','pb','bi','po','at','rn','fr','ra',
    	'ac','th','pa','u' );



process_arguments();

if ($actionhelp)
{
	display_help();
	exit 0;
}
if ($actionrep)
{
	status_report();
	exit 0;
}
if ($actionxmp)
{
	example_input();
	exit 0;
}
if ($actionvcl)
{
	clean_up_very();
	exit 0;
}
if ($actioncln)
{
	clean_up();
	exit 0;
}
if ($actiondel)
{
	clean_delete();
	exit 0;
}



read_physics();
read_proc();
occ_numbers();
get_ip();

if ($actiondir)
{
  make_dirs();
  $rerun_str=1;
}

if ($actionpre)
{
	run_prelim_autos() if (!-e "str/olg" or $rerun_str==1);
	get_prelim_autos();
        if($auger_damp > 0){
        run_auger();
        }
}

if ($actioninp)
{
#	if ($maxc eq 'auto')
#	{
#		auto_maxc();
#	}
	write_inputs_tcc();
	write_inputs_inner();
	write_inputs_dip();
	write_inputs_nonx();
	write_inputs_outer();
	write_inputs_tot();
	write_inputs_born();
#	write_inputs_adas();
}


if ($actiontcc)
{
	run_code_tcc();
}

if ($actioninn)
{
	run_code_inner();
}

if ($actiondip)
{
	run_code_dip();
}

if ($actionstgbinp)
{
  write_inputs_stgb();
}

if ($actionstgb)
{
  run_code_stgb();
}

if ($actionout)
{
	run_code_outer();
}

if ($actionnoi)
{
	run_code_nonxinn();
}

if ($actionnoo)
{
	run_code_nonxout();
}

if ($actionbrn)
{
	run_code_born();
	write_inputs_adas();
}

if ($actionmrg)
{
	make_omega_file();
}

if ($actionadf)
{
	make_adf_file();
}

if ($actionarc)
{
        archive_output();
}
