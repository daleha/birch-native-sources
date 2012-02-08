#!/usr/bin/env perl
use Env;
use FileHandle;
use Cwd;
use File::Path;
use Sys::Hostname;

########################## GLOBAL VARIABLES ################
our $PIDCHILD;
our $ERROR_DONE;
our @TMPFILE_LIST;
our $EXIT_FAILURE=1;
our $EXIT_SUCCESS=0;

our $REFDIR=getcwd;
our $EXIT_SUCCESS=0;
our $EXIT_FAILURE=1;

our $PROGRAM="tc_generic_method.pl";
our $CL=$PROGRAM;

our $CLEAN_EXIT_STARTED;
our $debug_lock=$ENV{"DEBUG_LOCK"};
our $debug_generic_method=$ENV{"DEBUG_GENERIC_METHOD"};
our $LOCKDIR=$ENV{"LOCKDIR_4_TCOFFEE"};
if (!$LOCKDIR){$LOCKDIR=getcwd();}
our $ERRORDIR=$ENV{"ERRORDIR_4_TCOFFEE"};
our $ERRORFILE=$ENV{"ERRORFILE_4_TCOFFEE"};
&set_lock ($$);
if (isshellpid(getppid())){lock4tc(getppid(), "LLOCK", "LSET", "$$\n");}
      
########################## GLOBAL VARIABLES ################


#Initiation


our $BLAST_MAX_NRUNS=2;
our $COMMAND;
our $PIDCHILD;

$REF_EMAIL="";
$tmp_dir="";
$init_dir="";


$test=0;
if ($test==1)
  {
    $SERVER="NCBI";
    $query=$ARGV[0];
    $hitf=$ARGV[1];
    %s=read_fasta_seq($query);
    @sl=keys(%s);
    &blast_xml2profile ("xx", $s{$sl[0]}{seq},$maxid,$minid,$mincov, $hitf);
    myexit ($EXIT_FAILURE);
  }

foreach $v(@ARGV){$cl.="$v ";}
$COMMAND=$cl;
($mode)=&my_get_opt ( $cl, "-mode=",1,0);

($A)=(&my_get_opt ( $cl, "-name1=",0,0));
($B)=(&my_get_opt ( $cl, "-name2=",0,0));
($TMPDIR)=(&my_get_opt ( $cl, "-tmpdir=",0,0));
($CACHE)=(&my_get_opt ( $cl, "-cache=",0,0));
($SERVER)=((&my_get_opt ( $cl, "-server=",0,0)));
($EMAIL)=((&my_get_opt ( $cl, "-email=",0,0)));

if (!$A){$A="A";}
if (!$B){$B="B";}


if (!$TMPDIR)
  {
    $HOME=$ENV{HOME};
    if ($ENV{TMP_4_TCOFFEE}){$TMPDIR=$ENV{TMP_4_TCOFFEE};}
    else{$TMPDIR="$HOME/.t_coffee/tmp/";}
  }
if ( ! -d $TMPDIR)
  {
    mkdir $TMPDIR;
  }
if ( ! -d $TMPDIR)
  {
    print "ERROR: Could not create temporary dir: $TMPDIR\n";
    myexit ($EXIT_FAILURE);
  }

$EMAIL=~s/XEMAILX/\@/g;
if (!$EMAIL)
  {
    if ($ENV{EMAIL_4_TCOFFEE}){$EMAIL=$ENV{EMAIL_4_TCOFFEE};}
    elsif ($ENV{EMAIL}){$EMAIL=$ENV{EMAIL};}
    else {$EMAIL=$REF_EMAIL;}
  }

($maxid,$minid,$mincov)=(&my_get_opt ( $cl, "-maxid=",0,0, "-minid=",0,0,"-mincov=",0,0));
if (!$cl=~/\-maxid\=/){$maxid=95;}
if (!$cl=~/\-minid\=/){$minid=35;}
if (!$cl=~/\-mincov\=/){$mincov=80;}




if ($mode eq "seq_msa")
  {
    &seq2msa($mode,&my_get_opt ( $cl, "-infile=",1,1, "-method=",1,2, "-param=",0,0,"-outfile=",1,0, "-database=",0,0));
  }
elsif ( $mode eq "tblastx_msa")
  {
    &seq2tblastx_lib ($mode,&my_get_opt ( $cl, "-infile=",1,1, "-outfile=",1,0));
  }
elsif ( $mode eq "tblastpx_msa")
  {
    &seq2tblastpx_lib ($mode,&my_get_opt ( $cl, "-infile=",1,1, "-outfile=",1,0));
  }
elsif ( $mode eq "thread_pair")
  {
    &seq2thread_pair($mode,&my_get_opt ( $cl, "-infile=",1,1, "-pdbfile1=",1,1, "-method=",1,2,"-param=",0,0, "-outfile=",1,0, ));
  }
elsif ( $mode eq "pdbid_pair")
  {
    &seq2pdbid_pair($mode,&my_get_opt ( $cl, "-pdbfile1=",1,0, "-pdbfile2=",1,0, "-method=",1,2,"-param=",0,0, "-outfile=",1,0, ));
  }
elsif ( $mode eq "pdb_pair")
  {
    &seq2pdb_pair($mode,&my_get_opt ( $cl, "-pdbfile1=",1,1, "-pdbfile2=",1,1, "-method=",1,2,"-param=",0,0, "-outfile=",1,0, ));
  }
elsif ( $mode eq "profile_pair")
  {
     &seq2profile_pair($mode,&my_get_opt ( $cl, "-profile1=",1,1, "-profile2=",1,1, "-method=",1,2,"-param=",0,0, "-outfile=",1,0 ));
  }
elsif ($mode eq "pdb_template_test")
  {
    &blast2pdb_template_test ($mode,&my_get_opt ( $cl, "-infile=",1,1));

  }
elsif ($mode eq "psi_template_test")
  {
    &psiblast2profile_template_test ($mode,&my_get_opt ( $cl, "-seq=",1,1,"-blast=",1,1));

  }

elsif ( $mode eq "pdb_template")
  {
    &blast2pdb_template ($mode,&my_get_opt ( $cl, "-infile=",1,1, "-database=",1,0, "-method=",1,0, "-outfile=",1,0,"-pdb_type=",1,0));
  }

elsif ( $mode eq "profile_template")
  {
    
    &psiblast2profile_template ($mode,&my_get_opt ( $cl, "-infile=",1,1, "-database=",1,0, "-method=",1,0, "-outfile=",1,0));
  }
elsif ( $mode eq "psiprofile_template")
  {
    &psiblast2profile_template ($mode,&my_get_opt ( $cl, "-infile=",1,1, "-database=",1,0, "-method=",1,0, "-outfile=",1,0));
  }
elsif ( $mode eq "RNA_template")
  {
    &seq2RNA_template ($mode,&my_get_opt ( $cl, "-infile=",1,1, "-outfile=",1,0));
  }
elsif ( $mode eq "tm_template")
  {
    &seq2tm_template ($mode, "", &my_get_opt ( $cl, "-infile=",1,1,"-arch=",1,1,"-psv=",1,1, "-outfile=",1,0,));
  }
elsif ( $mode eq "psitm_template")
  {
    &seq2tm_template ($mode,&my_get_opt ( $cl, "-database=",1,0, "-infile=",1,1, "-arch=",1,1,"-psv=",1,1, "-outfile=",1,0,));
  }
elsif ( $mode eq "ssp_template")
  {
    &seq2ssp_template ($mode,&my_get_opt ( $cl, "-infile=",1,1,"-seq=",1,1,"-obs=",1,1, "-outfile=",1,0));
  }
elsif ( $mode eq "psissp_template")
  {
    &seq2ssp_template ($mode,&my_get_opt ( $cl, "-infile=",1,1,"-seq=",1,1,"-obs=",1,1, "-outfile=",1,0));
  }

elsif ( $mode eq "rna_pair")
{
    &seq2rna_pair($mode,&my_get_opt ( $cl, "-pdbfile1=",1,1, "-pdbfile2=",1,1, "-method=",1,2,"-param=",0,0, "-outfile=",1,0, ));
}
elsif ( $mode eq "calc_rna_template")
{
    &calc_rna_template($mode,&my_get_opt ( $cl, "-infile=",1,1,"-pdbfile=",1,1, "-outfile=",1,0));
}
else
  {
    myexit(flush_error( "$mode is an unknown mode of tc_generic_method.pl"));
  }
myexit ($EXIT_SUCCESS);


sub seq2ssp_template
  {
  my ($mode, $infile,$gor_seq,$gor_obs,$outfile)=@_;
  my %s, %h;
  my $result;
  my (@profiles);
  &set_temporary_dir ("set",$infile,"seq.pep");
  %s=read_fasta_seq ("seq.pep");

  
  open (R, ">result.aln");
  
  #print stdout "\n";
  foreach $seq (keys(%s))
    {
      
      open (F, ">seqfile");
      $s{$seq}{seq}=uc$s{$seq}{seq};
      print (F ">$s{$seq}{name}\n$s{$seq}{seq}\n");
      close (F);
      $lib_name="$s{$seq}{name}.ssp";
      $lib_name=&clean_file_name ($lib_name);
      
      if ($mode eq "ssp_template"){&seq2gor_prediction ($s{$seq}{name},$s{$seq}{seq}, "seqfile", $lib_name,$gor_seq, $gor_obs);}
      elsif ($mode eq "psissp_template")
	{
	  &seq2msa_gor_prediction ($s{$seq}{name},$s{$seq}{seq},"seqfile", $lib_name,$gor_seq, $gor_obs);
	}
    
      if ( !-e $lib_name)
	{
	  myexit(flush_error("GORIV failed to compute the secondary structure of $s{$seq}{name}"));
	  myexit ($EXIT_FAILURE);
	}
      else
	{
	  print stdout "\tProcess: >$s{$seq}{name} _E_ $lib_name \n";
	  print R ">$s{$seq}{name} _E_ $lib_name\n";
	}
      unshift (@profiles, $lib_name);
    }
  close (R);
  &set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile, @profiles);
}

sub seq2tm_template
  {
  my ($mode, $db, $infile,$arch,$psv,$outfile)=@_;
  my %s, %h;
  my $result;
  my (@profiles);
  &set_temporary_dir ("set",$infile,"seq.pep");
  %s=read_fasta_seq ("seq.pep");

  
  open (R, ">result.aln");
  
  #print stdout "\n";
  foreach $seq (keys(%s))
    {
      open (F, ">seqfile");
      print (F ">$s{$seq}{name}\n$s{$seq}{seq}\n");
      close (F);
      $lib_name="$s{$seq}{name}.tmp";
      $lib_name=&clean_file_name ($lib_name);

      if ($mode eq "tm_template")
	{
	  &safe_system ("t_coffee -other_pg fasta_seq2hmmtop_fasta.pl -in=seqfile -out=$lib_name -arch=$arch -psv=$psv");
	}
      elsif ( $mode eq "psitm_template")
	{
	  &seq2msa_tm_prediction ($s{$seq}{name},$s{$seq}{seq}, $db, "seqfile", $lib_name,$arch, $psv);
	}
      if ( !-e $lib_name)
	{
	  myexit(flush_error("RNAplfold failed to compute the secondary structure of $s{$seq}{name}"));
	  myexit ($EXIT_FAILURE);
	}
      else
	{
	  print stdout "\tProcess: >$s{$seq}{name} _T_ $lib_name\n";
	  print R ">$s{$seq}{name} _T_ $lib_name\n";
	}
      unshift (@profiles, $lib_name);
    }
  close (R);
  &set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile, @profiles);
}

sub seq2RNA_template
  {
  my ($mode, $infile,$outfile)=@_;
  my %s, %h, ;
  my $result;
  my (@profiles);
  &set_temporary_dir ("set",$infile,"seq.pep");
  %s=read_fasta_seq ("seq.pep");

  
  open (R, ">result.aln");
  
  #print stdout "\n";
  foreach $seq (keys(%s))
    {
      open (F, ">seqfile");
      print (F ">$s{$seq}{name}\n$s{$seq}{seq}\n");
      close (F);
      $lib_name="$s{$seq}{name}.rfold";
      $lib_name=&clean_file_name ($lib_name);
      &safe_system ("t_coffee -other_pg RNAplfold2tclib.pl -in=seqfile -out=$lib_name");
      
      if ( !-e $lib_name)
	{
	 myexit(flush_error("RNAplfold failed to compute the secondary structure of $s{$seq}{name}"));
	  myexit ($EXIT_FAILURE);
	}
      else
	{
	  print stdout "\tProcess: >$s{$seq}{name} _F_ $lib_name\n";
	  print R ">$s{$seq}{name} _F_ $lib_name\n";
	}
      unshift (@profiles, $lib_name);
    }
  close (R);
  &set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile, @profiles);
}
sub psiblast2profile_template_test
  {
  my ($mode, $seq,$blast)=@_;
  my %s, %h, ;
  my ($result,$psiblast_output,$profile_name,@profiles);
  my $trim=0;
  my $maxid=100;
  my $minid=0;
  my $mincov=0;
  my $maxcov=100;
  
  %s=read_fasta_seq ($seq);
  open (R, ">result.aln");
  
  #print stdout "\n";
  foreach $seq (keys(%s))
    {
      
      open (F, ">seqfile");
      print (F ">$A\n$s{$seq}{seq}\n");
      close (F);
      $psiblast_output=$blast;
      if ( -e $psiblast_output)
	{
	  %profile=blast_xml2profile($s{$seq}{name}, $s{$seq}{seq},$maxid, $minid,$mincov,$psiblast_output);


	  
	  $profile_name="$s{$seq}{name}.prf";
	  $profile_name=&clean_file_name ($profile_name);
	  unshift (@profiles, $profile_name);
	  output_profile ($profile_name, \%profile, $trim);
	  print stdout "\tProcess: >$s{$seq}{name} _R_ $profile_name [$profile{n} Seq.] [$SERVER/blast/$db][$CACHE_STATUS]\n";
	  print R ">$s{$seq}{name} _R_ $profile_name\n";
	}
    }
  close (R);
  
  die;
}
sub psiblast2profile_template 
  {
  my ($mode, $infile, $db, $method, $outfile)=@_;
  my %s, %h, ;
  my ($result,$psiblast_output,$profile_name,@profiles);
  my $trim=0;
  &set_temporary_dir ("set",$infile,"seq.pep");
  %s=read_fasta_seq ("seq.pep");
  open (R, ">result.aln");
  
  #print stdout "\n";
  foreach $seq (keys(%s))
    {
      open (F, ">seqfile");
      print (F ">$A\n$s{$seq}{seq}\n");
      close (F);
      $psiblast_output=&run_blast ($s{$seq}{name},$method, $db, "seqfile","outfile");
      
if ( -e $psiblast_output)
	{
	  %profile=blast_xml2profile($s{$seq}{name}, $s{$seq}{seq},$maxid, $minid,$mincov,$psiblast_output);
	  unlink ($psiblast_output);
	  
	  $profile_name="$s{$seq}{name}.prf";
	  $profile_name=&clean_file_name ($profile_name);
	  unshift (@profiles, $profile_name);
	  output_profile ($profile_name, \%profile, $trim);
	  print stdout "\tProcess: >$s{$seq}{name} _R_ $profile_name [$profile{n} Seq.] [$SERVER/blast/$db][$CACHE_STATUS]\n";
	  print R ">$s{$seq}{name} _R_ $profile_name\n";
	}
    }
  close (R);
  &set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile, @profiles);
}
sub blast2pdb_template_test
    {
      my ($mode,$infile)=@_;
      my ($maxid,$minid,$mincov);
      $maxid=100;
      $minid=0;
      $mincov=0;
      
      print "$infile\n";
      
      %p=blast_xml2profile($s{$seq}{name}, $s{$seq}{seq},$maxid, $minid,$mincov,$infile);
      $c=1;
      print stdout "\tProcess: >$s{$seq}{name} [$SERVER/blast/$db][$CACHE_STATUS]\n";
      while (!$found && $c<$p{n})
	{
	  $pdbid=&id2pdbid($p{$c}{identifyer});
	  if ( length ($pdbid)>5){$pdbid=id2pdbid($p{$c}{definition});}
	  
	  if ( length ($pdbid)>5)
	    {
	      myexit(add_error (EXIT_FAILURE,$$,$$,getppid(), "BLAST_FAILURE::Could Not Parse PDBID ($p{$c}{identifyer},$p{$c}{definition})"));
	    }
	  
	  
	  if (!&pdb_is_released($pdbid))
	    {
	      print stdout "\t\t**$pdbid [PDB NOT RELEASED or WITHDRAWN]\n";
	      $c++;
	    }
	  elsif (!&pdb_has_right_type ($pdbid,$type))
	    {
	      my $ptype=&pdb2type ($pdbid);
	      my $etype=&type2etype($type);
	      
	      print stdout "\t\t**$pdbid [$ptype cannot be used (expected: $etype)]\n";
	      $c++;
	    }
	  else
	    {
	      $found=1;
	    }
	}

      if ($found)
	{
	  print stdout "\t\t >$s{$seq}{name} _P_ $pdbid\n";
	}
      else
	{
	  print stdout "\t\t >$s{$seq}{name} No Template Selected\n";
	}
      die;
    }
sub blast2pdb_template 
  {
  my ($mode, $infile, $db, $method, $outfile,$type)=@_;
  my %s, %h, ;
  my ($result,$blast_output);
  &set_temporary_dir ("set",$infile,"seq.pep");
  %s=read_fasta_seq ("seq.pep");
  open (R, ">result.aln");
  
 
  #print stdout "\n";
  foreach $seq (keys(%s))
    {
      my $c;
      my $found;
      
      open (F, ">seqfile");
      print (F ">$A\n$s{$seq}{seq}\n");
      close (F);
     
      $blast_output=&run_blast ($s{$seq}{name},$method, $db, "seqfile","outfile");
     
      %p=blast_xml2profile($s{$seq}{name}, $s{$seq}{seq},$maxid, $minid,$mincov,$blast_output);
      unlink ($blast_output);
      
      $c=1;
      print stdout "\tProcess: >$s{$seq}{name} [$SERVER/blast/$db][$CACHE_STATUS]\n";
      while (!$found && $c<$p{n})
	{
	  $pdbid=&id2pdbid($p{$c}{identifyer});
	  if ( length ($pdbid)>5){$pdbid=id2pdbid($p{$c}{definition});}

	  if ( length ($pdbid)>5)
	    {
	      myexit(add_error (EXIT_FAILURE,$$,$$,getppid(), "BLAST_FAILURE::Could Not Parse PDBID ($p{$c}{identifyer},$p{$c}{definition})"));
	    }
	  

	  if (!&pdb_is_released($pdbid))
	    {
	      print stdout "\t\t**$pdbid [PDB NOT RELEASED or WITHDRAWN]\n";
	      $c++;
	    }
	  elsif (!&pdb_has_right_type ($pdbid,$type))
	    {
	      my $ptype=&pdb2type ($pdbid);
	      my $etype=&type2etype($type);
	      
	      print stdout "\t\t**$pdbid [$ptype cannot be used (expected: $etype)]\n";
	      $c++;
	    }
	  else
	    {
	      $found=1;
	    }
	}

      if ($found)
	{
	  print R ">$s{$seq}{name} _P_ $pdbid\n";
	  print stdout "\t\t >$s{$seq}{name} _P_ $pdbid\n";
	}
      else
	{
	  print R ">$s{$seq}{name}\n";
	  print stdout "\t\t >$s{$seq}{name} No Template Selected\n";
	}
    }
  close (R);
  &set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile);
}
sub type2etype
  {
    my $type=shift;
    my $etype;
    
    if ( $type=~/n/){$etype.="NMR ";}
    if ( $type=~/d/){$etype.="diffraction ";}
    if ( $type=~/m/){$etype.="model ";}
    return $etype;
  }
sub pdb2type
  {
     my $pdb=shift;
     my $f=vtmpnam();
     
     my $value= &safe_system ("t_coffee -other_pg extract_from_pdb -model_type $pdb > $f");
     my $r=&file2string ($f);
     chomp($r);
     return $r;
   }
sub pdb_has_right_type
  {
    my $pdb=shift;
    my $type=shift;
    
    my $f=vtmpnam();
    
    my $value= &safe_system ("t_coffee -other_pg extract_from_pdb -model_type $pdb > $f");
    my $r=&file2string ($f);
    chomp($r);

        
    if ( $r eq "NMR" && $type=~/n/){return 1;}
    elsif ( $r eq "diffraction" && $type=~/d/){return 1;}
    elsif ( $r eq "model" && $type=~/m/){return 1;}
    else {return 0;}
  }
sub pdb_is_released
  {
    my $pdb=shift;
    my $f=vtmpnam();
    
    $value= &safe_system ("t_coffee -other_pg extract_from_pdb -is_released_pdb_name $pdb > $f");
    my $r=&file2string ($f);
    chomp($r);
    return $r;
  }
sub blast_msa
  {
    my ($blast,$infile,$db,$outfile)=@_;
    my ($a, %s1, %s, %qs, %qs1);
    my $seqfile;
    my $SEQ=new FileHandle;
    my $seqfile="seqfile";
    my @txt;
    
    
    %s1=&read_fasta_seq ($db);
    %s=&fasta_hash2index_hash(%s1);
    %qs1=&read_fasta_seq ($infile);
    %qs=&fasta_hash2index_hash(%qs1);
    
    
    #&safe_system ("formatdb -i $db");
    if ($blast eq "blastp"){&safe_system  ("blastall -i $infile -d $db -m7 -p blastp -o io");}
    elsif ($blast eq "blastn"){&safe_system  ("blastn -query $infile -db $db -outfmt 5 -word_size 4 -out io");}

    &set_blast_type ("io");
    

    my %FB=&xml2tag_list ("io", "Iteration");
    open (F, ">$outfile");
    print F "! TC_LIB_FORMAT_01\n";
    print F "$s{n}\n";
    for ( my $a=0; $a<$s{n}; $a++)
      {
	print F "$s{$a}{name} $s{$a}{len} $s{$a}{seq}\n";
      }


    for ( my $a=0; $a<$FB{n}; $a++)
      {
	my %p=blast_xml2profile ($qs{$a}{name}, $qs{$a}{seq},100, 0, 0, $FB{$a}{body});
	my $query=$p{0}{name};
	my $i= $s1{$query}{order}+1;
	for (my $b=1; $b<$p{n}; $b++)
	  {
	    my $l=length ($p{$b}{Qseq});
	    my $hit=$p{$b}{definition};
	    my $Qstart=$p{$b}{Qstart};
	    my $Hstart=$p{$b}{Hstart};
	    my $identity=$p{$b}{identity};
	    my @lrQ=split (//,$p{$b}{Qseq});
	    my @lrH=split (//,$p{$b}{Hseq});
	    
	    my $j= $s1{$hit}{order}+1;
	    #if ( $j==$i){next;}
	    printf F "# %d %d\n", $i, $j;
	    #  print  F "\n$p{$b}{Qseq} ($Qstart)\n$p{$b}{Hseq} ($Hstart)";
	    for ($c=0; $c<$l; $c++)
	      {
		my $rQ=$lrQ[$c];
		my $rH=$lrH[$c];
		my $n=0;
		
		if ($rQ ne "-"){$n++, $Qstart++;}
		if ($rH ne "-"){$n++; $Hstart++;}
		
		if ( $n==2)
		  {
		    printf F "\t%d %d %d\n", $Qstart-1, $Hstart-1,$identity;
		  }
	      }
	  }
      }
    print F "! SEQ_1_TO_N\n";
    close (F);
    return $output;
  }

sub blast_msa_old
  {
    my ($infile,$outfile)=@_;
    my ($a, %seq);
    %s1=&read_fasta_seq ($infile);
    foreach $s (keys (%s1))
      {
	$i=$s1{$s}{order};
	$s{$i}{name}=$s;
	$s{$i}{seq}=$s1{$s}{seq};
	$s{$i}{len}=length( $s{$i}{seq});
	$s{n}++;
      }
    &safe_system ("formatdb -i $infile");
    &safe_system ("blastall -i $infile -d $infile -m7 -o io");
    &set_blast_type ("io");
    
    %FB=&xml2tag_list ("io", "Iteration");
    
    open (F, ">$outfile");
    print F "! TC_LIB_FORMAT_01\n";
    print F "$s{n}\n";
    for ( $a=0; $a<$s{n}; $a++)
      {
	print F "$s{$a}{name} $s{$a}{len} $s{$a}{seq}\n";
      }
    for ( $a=0; $a<$FB{n}; $a++)
      {
	%p=blast_xml2profile ($s{$a}{name}, $s{$a}{seq},100, 0, 0, $FB{$a}{body});
	for ($b=1; $b<$p{n}; $b++)
	  {
	    my $l=length ($p{$b}{Qseq});
	    my $hit=$p{$b}{definition};
	    my $Qstart=$p{$b}{Qstart};
	    my $Hstart=$p{$b}{Hstart};
	    my $identity=$p{$b}{identity};
	    my @lrQ=split (//,$p{$b}{Qseq});
	    my @lrH=split (//,$p{$b}{Hseq});
	    my $i= $s1{$s{$a}{name}}{order}+1;
	    my $j= $s1{$hit}{order}+1;
	    #if ( $j==$i){next;}
	    printf F "# %d %d\n", $i, $j;
	    #  print  F "\n$p{$b}{Qseq} ($Qstart)\n$p{$b}{Hseq} ($Hstart)";
	    for ($c=0; $c<$l; $c++)
	      {
		my $rQ=$lrQ[$c];
		my $rH=$lrH[$c];
		my $n=0;
		
		if ($rQ ne "-"){$n++, $Qstart++;}
		if ($rH ne "-"){$n++; $Hstart++;}
		
		if ( $n==2)
		  {
		    printf F "\t%d %d %d\n", $Qstart-1, $Hstart-1,$identity;
		  }
	      }
	  }
      }
    print F "! SEQ_1_TO_N\n";
    close (F);
    return $output;
  
  }

sub seq2msa
  {
    my ($mode, $infile, $method, $param, $outfile,$database)=@_;
    &set_temporary_dir ("set",$infile,"seq.pep", $database, "db.pep");
    $param.=" >/dev/null 2>&1 ";
    
    
    #make sure test.pep is in FASTA
    &safe_system ("t_coffee -other_pg seq_reformat -in seq.pep -output fasta_seq > x");
    `mv x seq.pep`;
    
    if ( $method eq "blastp")
      {
	&blast_msa ("blastp","seq.pep",$database,"result.aln");
      }
    elsif ( $method eq "blastn")
      {
	&blast_msa ("blastn","seq.pep",$database,"result.aln");
      }
    
    elsif ( $method eq "muscle")
      {
	`muscle -in seq.pep -out result.aln $param`;
      }
    elsif ( $method eq "probcons")
      {
	`probcons seq.pep >result.aln 2>/dev/null`;
      }
    elsif ( $method eq "mafft")
      {
	`mafft --quiet --localpair --maxiterate 1000 seq.pep> result.aln  2>/dev/null`
      }
    elsif ( $method=~/prank/)
      {
	`$method -d=seq.pep -o=result.aln -quiet 2>/dev/null`;
	`mv result.aln.1.fas result.aln`;
      }
    elsif ($method eq "clustalo")
      {
	`clustalo -i seq.pep > result.aln`;
      }
    else
      {
	`$method -infile=seq.pep -outfile=result.aln`;
      }
    
    &set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile);
    myexit ($EXIT_SUCCESS);
  }

sub seq2thread_pair
  {
    my ($mode, $infile, $pdbfile1, $method, $param, $outfile)=@_;
    &set_temporary_dir ("set",$infile,"seq.pep",$pdbfile1,"struc.pdb");
    if ($method eq "fugueali")
      {
	#Env Variable that need to be defined for Fugue
	if (!$ENV{FUGUE_LIB_LIST}){$ENV{FUGUE_LIB_LIST}="DUMMY";}
	if (!$ENV{HOMSTRAD_PATH})  {$ENV{HOMSTRAD_PATH}="DUMMY";}
	if (!$ENV{HOMS_PATH}){$ENV{HOMS_PATH}="DUMMY";}
	
	`joy struc.pdb >x 2>x`;
	&check_file("struc.tem", "Joy failed [FATAL:$PROGRAM/$method]");
	`melody -t struc.tem >x 2>x`;
	&check_file("struc.tem", "Melody failed [FATAL:$PROGRAM/$method]");
	`fugueali -seq seq.pep -prf struc.fug -print > tmp_result.aln`;
	
	&check_file("tmp_result.aln", "Fugue failed [FATAL:$PROGRAM/$method]");
	&safe_system ("t_coffee -other_pg seq_reformat -in tmp_result.aln -output fasta_aln >result.aln");
      }
    elsif ( $method eq "t_coffee")
      {
	&safe_system ("t_coffee -in Pstruc.pdb Sseq.pep Mslow_pair -outfile result.aln -quiet");
      }
    else
      {
	&safe_system ("$method -infile=seq.pep -pdbfile1=struc.pdb -outfile=result.aln $param>x 2>x");
      }
    &set_temporary_dir ("unset",$mode,$method,"result.aln",$outfile);
    myexit ($EXIT_SUCCESS);
  }
sub seq2pdbid_pair
  {
    my ($mode, $pdbfile1, $pdbfile2, $method, $param, $outfile)=@_;
    my ($name);

    
    &set_temporary_dir ("set");
    $name=$pdbfile1." ".$pdbfile2;

    if (    &cache_file("GET","","$name","$method","dali",$outfile,"EBI"))
      {return $outfile;}
    else
      {
	if ($method eq "daliweb")
	  {
	    $pdbfile1=~/(....)(.)/;
	    $id1=$1; $c1=$2;
	    
	    $pdbfile2=~/(....)(.)/;
	    $id2=$1; $c2=$2;
	    
	    $command="t_coffee -other_pg dalilite.pl --pdb1 $id1 --chainid1 $c1 --pdb2 $id2 --chainid2 $c2 --email=$EMAIL  >dali_stderr 2>dali_stderr";
	    $dali=`$command`;
	    
	    open (F, "dali_stderr");
	    while (<F>)
	      {
		if ( /JobId: dalilite-(\S+)/)
		{
		  $jobid=$1;
		}
	      }
	    close (F);
	    unlink ("dali_stderr");
	    
	    $output1="dalilite-$jobid.txt";
	    if ( -e $output1)
	      {
		unlink ($output1);
		&url2file ("http://www.ebi.ac.uk/Tools/es/cgi-bin/jobresults.cgi/dalilite/dalilite-$jobid/aln.html", "output2");
		
		if ( -e "output2")
		  {
		    my ($seq1, $seq2);
		    $seq1=$seq2="";
		    
		    open (F, "output2");
		    while (<F>)
		      {
			$l=$_;
			if ( $l=~/Query\s+(\S+)/)
			  {
			    $seq1.=$1;
			  }
			elsif ( $l=~/Sbjct\s+(\S+)/)
			  {
			    $seq2.=$1;
			  }
		      }
		    close (F);
		    unlink ("output2");
		    if ($seq1 ne "" && $seq2 ne "")
		      {
			$output3=">$A\n$seq1\n>$B\n$seq2\n";
			$output3=~s/\./-/g;
			open (F, ">result.aln");
			print F "$output3";
			close (F);
		      }
		  }
	      }
	  }
      }
    &cache_file("SET","","$name","$method","dali","result.aln","EBI");
    &set_temporary_dir ("unset",$mode, $method, "result.aln",$outfile);
    myexit ($EXIT_SUCCESS);
  }
sub seq2pdb_pair
  {
    my ($mode, $pdbfile1, $pdbfile2, $method, $param, $outfile)=@_;
    
    &set_temporary_dir ("set",$pdbfile1,"pdb1.pdb",$pdbfile2,"pdb2.pdb");
    if ($method eq "t_coffee")
      {
	&safe_system ("t_coffee -in Ppdb1.pdb Ppdb2.pdb -quiet -outfile=result.aln");
      }
    elsif ( $method eq "DaliLite")
      {
	if ( &safe_system ("DaliLite -pairwise pdb1.pdb pdb2.pdb >tmp1")==$EXIT_SUCCESS)
	  {
	     my ($seq1, $seq2);
	     $seq1=$seq2="";
		    
	     open (F, "tmp1");
	     while (<F>)
	       {
		 $l=$_;
		 if ( $l=~/Query\s+(\S+)/)
		   {
		     $seq1.=$1;
		   }
		 elsif ( $l=~/Sbjct\s+(\S+)/)
		   {
		     $seq2.=$1;
		   }
	       }
	     close (F);
	     unlink ("tmp1");
	     if ($seq1 ne "" && $seq2 ne "")
	       {
		 my $output3=">$A\n$seq1\n>$B\n$seq2\n";
		 $output3=~s/\./-/g;
		 open (F, ">result.aln");
		 print F "$output3";
		 close (F);
	       }
	   }
	else
	  {
	    print "ERROR: DalLite failed to align the considered structures[tc_generic_method.pl]\n";
	  }    
      }
    elsif ( $method eq "TMalign")
      {
	if ( &safe_system ("TMalign pdb1.pdb pdb2.pdb >tmp1")==$EXIT_SUCCESS)
	  {
	    `tail -4 tmp1 > tmp2`;
	    
	    open (F, "tmp2");
	    while (<F>)
	      {
		unshift(@l, $_);
	      }
	    close (F);
	    open (F, ">result.aln");
	    $l[3]=~s/[^a-zA-Z0-9-]/\-/g;
	    $l[1]=~s/[^a-zA-Z0-9-]/\-/g;
	    print F ">$A\n$l[3]\n>$B\n$l[1]\n";
	    close (F);
	  }
	else
	  {
	    print "ERROR: TMalign failed to align the considered structures[tc_generic_method.pl]\n";
	    `rm result.aln >/dev/null 2>/dev/null`;
	  }
      }
    elsif ( $method eq "mustang")
      {
	if ( &safe_system ("mustang -i pdb1.pdb pdb2.pdb -F fasta >/dev/null 2>/dev/null")==$EXIT_SUCCESS)
	  {
	    `mv results.afasta result.aln`;
	  }
	else
	  {
	    print "ERROR: mustang failed to align the considered structures[tc_generic_method.pl]\n";
	    `rm result.aln >/dev/null 2>/dev/null`;
	  }
      }
    else
      {
	if ( &safe_system ("$method -pdbfile1=pdb1.pep -pdbfile2=pdb2.pdb -outfile=result.aln $param>x 2>x")==$EXIT_SUCCESS)
	  {
	    `mv results.afasta result.aln`;
	  }
	else
	  {
	    print "ERROR: $method failed to align the considered structures[tc_generic_method.pl]\n";
	    `rm result.aln >/dev/null 2>/dev/null`;
	  }
      }
    &set_temporary_dir ("unset",$mode, $method, "result.aln",$outfile);
    myexit ($EXIT_SUCCESS);
  }

sub seq2profile_pair
  {
    my ($mode, $profile1, $profile2, $method, $param, $outfile)=@_;
    
    
    if ($method eq "clustalw")
      {
	&set_temporary_dir ("set",$profile1,"prf1.aln",$profile2,"prf2.aln");
	`clustalw -profile1=prf1.aln -profile2=prf2.aln -outfile=result.aln`;
	&set_temporary_dir ("unset",$mode, $method, "result.aln",$outfile);
      }
    elsif ( $method eq "hhalign")
      {
	hhalign ( $profile1,$profile2,$outfile,$param);
      }
    else
      {
	
	`$method -profile1=prf1.aln -profile2=prf2.aln -outfile=result.aln $param>x 2>x`;
      }
   
    myexit ($EXIT_SUCCESS);
  }

sub pg_is_installed
  {
    my @ml=@_;
    my ($r, $p, $m);
    my $supported=0;
    
    my $p=shift (@ml);
    if ($p=~/::/)
      {
	if (safe_system ("perl -M$p -e 1")==$EXIT_SUCCESS){return 1;}
	else {return 0;}
      }
    else
      {
	$r=`which $p 2>/dev/null`;
	if ($r eq ""){$r=0;}
	else {$r=1;}
	
	if ($r==0 && is_blast_package ($p)){return pg_is_installed ("legacy_blast.pl");}
	else {return $r;}
      }
  }

sub is_blast_package
  {
    my $p=shift;
    if ( $p=~/blastp/){return 1;}
    elsif ($p=~/blastall/){return 1;}
    elsif ($p=~/blastn/){return 1;}
    elsif ($p=~/blastx/){return 1;}
    elsif ($p=~/formatdb/){return 1;}
    else {return 0;}
  }
    
sub check_internet_connection
  {
    my $internet;
    my $tmp;
    &check_configuration ( "wget"); 
    
    $tmp=&vtmpnam ();
    
    if     (&pg_is_installed    ("wget")){`wget www.google.com -O$tmp >/dev/null 2>/dev/null`;}
    elsif  (&pg_is_installed    ("curl")){`curl www.google.com -o$tmp >/dev/null 2>/dev/null`;}
    
    if ( !-e $tmp || -s $tmp < 10){$internet=0;}
    else {$internet=1;}
    if (-e $tmp){unlink $tmp;}

    return $internet;
  }
sub check_pg_is_installed
  {
    my @ml=@_;
    my $r=&pg_is_installed (@ml);
    if (!$r && $p=~/::/)
      {
	print STDERR "\nYou Must Install the perl package $p on your system.\nRUN:\n\tsudo perl -MCPAN -e 'install $pg'\n";
      }
    elsif (!$r)
      {
	myexit(flush_error("\nProgram $p Supported but Not Installed on your system"));
      }
    else
      {
	return 1;
      }
  }
sub set_temporary_dir
  {
    my @list=@_;
    my $dir_mode, $a, $mode, $method;
  
    $dir_mode=shift (@list);

    
    if ( $dir_mode eq "set")
      {
	$initial_dir=cwd();
	if ( !$tmp_dir)
	  {
	    $rand=rand (100000);
	    $tmp_dir="$TMPDIR/tmp4tcoffee_profile_pair_dir_$$\_P_$rand";
	  }
	if ( !-d $tmp_dir)
	  {
	    push (@TMPDIR_LIST, $tmp_dir);
	    `mkdir $tmp_dir`;
	  }
	
	for ( $a=0; $a<=$#list; $a+=2)
	      {
		if (-e $list[$a]){ `cp $list[$a] $tmp_dir/$list[$a+1]`;}
	      }
	chdir $tmp_dir;
      }
    elsif ( $dir_mode eq "unset")
      {
	$mode=shift (@list);
	$method=shift (@list);
	
	if (!-e $list[0])
	  {
	   myexit(flush_error("Program $method failed to produce $list[1]" ));
	    myexit ($EXIT_FAILURE);
	  }
	else
	  {
	    chdir $initial_dir;
	    # `t_coffee -other_pg seq_reformat -in $tmp_dir/$list[0] -output fasta_aln -out $tmp_dir/result2.aln`;
	    `cp $tmp_dir/$list[0] $tmp_dir/result2.aln`;
	    if ( $list[1] eq "stdout")
	      {
		open (F, "$tmp_dir/result2.aln");
		while (<F>){print $_;}close(F);
	      }
	    else
	      {
		`mv $tmp_dir/result2.aln $list[1]`;
	      }
	    shift (@list); shift (@list);
	    foreach $f (@list)
	      {
		if (-e ("$tmp_dir/$f")){`mv $tmp_dir/$f .`;}
	      }
	  }
      }
  }




sub my_get_opt
  {
    my @list=@_;
    my $cl, $a, $argv, @argl;
    
    @argl=();
    $cl=shift @list;
    for ( $a=0; $a<=$#list; $a+=3)
      {
	$option=$list[$a];
	$optional=$list[$a+1];
	$status=$list[$a+2];
	$argv="";
	if ($cl=~/$option(\S+)/){$argv=$1;}
	@argl=(@argl,$argv);
	
	
	#$optional:0=>optional
	#$optional:1=>must be set
	#$status: 0=>no requirement
	#$status: 1=>must be an existing file
	#$status: 2=>must be an installed package
	

	if ($optional==0){;}
	elsif ( $optional==1 && $argv eq "")
	  {
	    myexit(flush_error( "ERROR: Option $option must be set"));
	    myexit ($EXIT_FAILURE);
	  }
	if ($status==0){;}
	elsif ($status ==1 && $argv ne "" && !-e $argv)
	  {
	    myexit(flush_error( "File $argv must exist"));
	    myexit ($EXIT_FAILURE);
	  }
	elsif ( $status==2 && $argv ne "" && &check_pg_is_installed ($argv)==0)
	  {
	    myexit(flush_error( " $argv is not installed"));
	    myexit ($EXIT_FAILURE);
	  }
      }

    return @argl;
    }

sub check_file 
  {
    my ($file, $msg)=@_;

    if ( !-e $file)
      {
	myexit(flush_error("$msg"));
      }
    }
sub hhalign
  {
    my ($aln1, $aln2, $outfile, $param)=@_;
    my $h1, $h2;
    
    $h{0}{index}=0;
    $h{1}{index}=1;
    
    $h{0}{aln}=$aln1;
    $h{1}{aln}=$aln2;

   

    %{$h{0}}=aln2psi_profile (%{$h{0}});
    %{$h{1}}=aln2psi_profile (%{$h{1}});

    $param=~s/#S/ /g;
    $param=~s/#M/\-/g;
    $param=~s/#E/\=/g;
    

    
    $command="hhalign -i $h{0}{a3m} -t $h{1}{a3m} -tc $outfile.tmp -rank 1 -mapt 0 $param";
    `$command`;
    
  #  `hhalign -i $h{0}{a3m} -t $h{1}{a3m} -tc $outfile.tmp -rank 1 -mapt 0 -gapf 0.8 -gapg 0.8`;
    

    # To run global use the following
    
    open (I, "$outfile.tmp");
    open (O, ">$outfile");
    $h{0}{cons}=s/\./x/g;
    $h{1}{cons}=s/\./x/g;

    print O "! TC_LIB_FORMAT_01\n2\n$h{0}{name} $h{0}{len} $h{0}{seq}\n$h{1}{name} $h{1}{len} $h{1}{seq}\n#1 2\n";
    
    while (<I>)
      {
	if (/(\d+)\s+(\d+)\s+(\d+)/)
	  {
	    print O "\t$h{0}{$1}\t$h{1}{$2}\t$3\n";
	  }
      }
    print O "! SEQ_1_TO_N\n";

    close (O);
    close (I);
  }

sub aln2psi_profile
  {
    my (%h)=@_;
    my ($aln,$i,$hv, $a, @c, $n);
   
    
    $i=$h{index};
    $aln=$h{aln};
    
    `cp $aln $$.hhh_aln`;
    $command="t_coffee -other_pg seq_reformat -in $aln -output hasch";
    $hv=`$command`;chomp ($hv);
    
    $h{a2m}="$tmp/$hv.tmp4hhpred.a2m";
    $h{a3m}="$tmp/$hv.tmp4hhpred.a3m";
    if ( -e $h{a3m}){;}
    else
      {
	$x=`which hhconsensus`;
	`hhconsensus  -M 50 -i $h{aln} -oa2m $h{a2m}`;
	if (!-e $h{a2m})
	  {
	    print STDERR "Program tc_generic_method.pl FAILED to run:\n\thhconsensus  -M 50 -i $h{aln} -oa2m $h{a2m}";
	    myexit ($EXIT_FAILURE);
	  }
	
	`hhconsensus  -M 50 -i $h{aln} -oa3m $h{a3m}`;
	if (!-e $h{a3m})
	  {
	    print STDERR "Program tc_generic_method.pl FAILED to run:\n\thhconsensus  -M 50 -i $h{aln} -oa3m $h{a3m}";
	    myexit ($EXIT_FAILURE);
	  }
       `buildali.pl $h{a3m} -n 1`;
      }
    
    
    $h{a2m_seq}=`head -n 2 $h{a2m} | grep -v ">"`;chomp ($h{a2m_seq});
    $h{a3m_seq}=`head -n 2 $h{a3m} | grep -v ">"`;chomp ($h{a3m_seq});
    $h{cons}=$h{a2m_seq};
    $h{seq}=`head -n 2 $h{aln} | grep -v ">"`;chomp ($h{seq});
    
    

    @c=split (//, $h{cons});
    $h{len}=$#c+1;
    for ($n=0,$a=0, $b=0; $a<$h{len};$a++)
      {
	if ( $c[$a]=~/[A-Z]/)
	  {
	    $h{++$n}=++$b;

	  }
	elsif ( $c[$a]=~/[a-z\.]/)
	  {
	    ++$b;
	  }
      }
    
    $name=`head -n 2 $h{aln} | grep ">"`;
    $name=~/\>(\S+)/;
    $h{name}=$1;
    
    `cp $h{a2m} $i.a2m`;
    `cp $h{a3m} $i.a3m`;
    `cp $h{aln} $i.hh_aln`;
    
    return %h;
  }

sub read_fasta_seq 
  {
    my $f=@_[0];
    my %hseq;
    my (@seq, @com, @name);
    my ($a, $s,$nseq);

    open (F, $f);
    while (<F>)
      {
	$s.=$_;
      }
    close (F);

    
    @name=($s=~/>(\S*).*\n[^>]*/g);
    
    @seq =($s=~/>.*.*\n([^>]*)/g);
    @com =($s=~/>\S*(.*)\n([^>]*)/g);

    
    $nseq=$#name+1;
    
    for ($a=0; $a<$nseq; $a++)
      {
	my $s;
	my $n=$name[$a];
	$hseq{$n}{name}=$n;
	$seq[$a]=~s/[^A-Za-z]//g;
	$hseq{$n}{order}=$a;
	$hseq{$n}{seq}=$seq[$a];
	$hseq{$n}{com}=$com[$a];
	
      }
    return %hseq;
  }
sub fasta_hash2index_hash
  {
    my %s1=@_;
    my %s;
    foreach my $s (keys (%s1))
      {
	my $i=$s1{$s}{order};
	$s{$i}{name}=$s;
	$s{$i}{seq}=$s1{$s}{seq};
	$s{$i}{len}=length( $s{$i}{seq});
	$s{n}++;
      }
    return %s;
  }
sub file_contains 
  {
    my ($file, $tag, $max)=(@_);
    my ($n);
    $n=0;
    
    if ( !-e $file && ($file =~/$tag/)) {return 1;}
    elsif ( !-e $file){return 0;}
    else 
      {
	open (FC, "$file");
	while ( <FC>)
	  {
	    if ( ($_=~/$tag/))
	      {
		close (FC);
		return 1;
	      }
	    elsif ($max && $n>$max)
	      {
		close (FC);
		return 0;
	      }
	    $n++;
	  }
      }
    close (FC);
    return 0;
  }
	    
	  
sub file2string
  {
    my $f=@_[0];
    my $string, $l;
    open (F,"$f");
    while (<F>)
      {

	$l=$_;
	#chomp ($l);
	$string.=$l;
      }
    close (F);
    $string=~s/\r\n//g;
    $string=~s/\n//g;
    return $string;
  }


sub my_get_opt
  {
    my @list=@_;
    my $cl, $a, $argv, @argl;
    
    @argl=();
    $cl=shift @list;
    for ( $a=0; $a<=$#list; $a+=3)
      {
	$option=$list[$a];
	$optional=$list[$a+1];
	$status=$list[$a+2];
	$argv="";
	if ($cl=~/$option(\S+)/){$argv=$1;}
	@argl=(@argl,$argv);
	
	
	#$optional:0=>optional
	#$optional:1=>must be set
	#$status: 0=>no requirement
	#$status: 1=>must be an existing file
	#$status: 2=>must be an installed package
	

	if ($optional==0){;}
	elsif ( $optional==1 && $argv eq "")
	  {

	    myexit(flush_error("Option $option must be set"));
	   
	  }
	if ($status==0){;}
	elsif ($status ==1 && $argv ne "" && !-e $argv)
	  {
	     myexit(flush_error("File $argv must exist"));
	   
	  }
	elsif ( $status==2 && $argv ne "" && &check_pg_is_installed ($argv)==0)
	  {
	    myexit(flush_error("$argv is not installed"));
	   
	  }
      }

    return @argl;
    }

sub tag2value 
  {
    
    my $tag=(@_[0]);
    my $word=(@_[1]);
    my $return;
    
    $tag=~/$word="([^"]+)"/;
    $return=$1;
    return $return;
  }
      
sub hit_tag2pdbid
  {
    my $tag=(@_[0]);
    my $pdbid;
       
    $tag=~/id="(\S+)"/;
    $pdbid=$1;
    $pdbid=~s/_//;
    return $pdbid;
  }
sub id2pdbid
  {
    my $in=@_[0];
    my $id;
    
    $in=~/(\S+)/;
    $id=$in;
    $id=~s/PDB/pdb/g;
    
    if ($id =~/pdb(.*)/){$id=$1;}
    elsif ( $id=~/(\S+)\s+mol:protein/){$id=$1;}
    $id=~s/[:|��_]//g;
    return $id;
  }
sub set_blast_type 
  {
    my $file =@_[0];
    if (&file_contains ($file,"EBIApplicationResult",100)){$BLAST_TYPE="EBI";}
    elsif (&file_contains ($file,"NCBI_BlastOutput",100)) {$BLAST_TYPE="NCBI";}
    else
      {
	$BLAST_TYPE="";
      }
    return $BLAST_TYPE;
  }
sub is_valid_blast_xml
    {
      my $file=shift;
      my $line;
      
      
      if ( !-e $file) {return 0;}
      $line=&file2tail ($file,100);
      
      if ( $line=~/<\/EBIApplicationResult/ || $line=~/<\/NCBI_BlastOutput/ || $line=~/<\/BlastOutput/ ){return 1;}
      return 0;
    }
sub file2blast_flavor
      {
	my $file=shift;
	if (&file_contains ($file,"EBIApplicationResult",100)){return "EBI";}
	elsif (&file_contains ($file,"NCBI_BlastOutput",100)){return "NCBI";}
	else {return "UNKNOWN";}
      }
sub blast_xml2profile 
  {
    my ($name,$seq,$maxid, $minid, $mincov, $file)=(@_);
    my (%p, $a, $string, $n);
    
    

    if ($BLAST_TYPE eq "EBI" || &file_contains ($file,"EBIApplicationResult",100)){%p=ebi_blast_xml2profile(@_);}
    elsif ($BLAST_TYPE eq "NCBI" || &file_contains ($file,"NCBI_BlastOutput",100)){%p=ncbi_blast_xml2profile(@_);}
    else 
      {
	myexit(add_error ( $$,$$,getppid(), "BLAST_FAILURE::unkown XML",$CL));
      }
    for ($a=0; $a<$p{n}; $a++)
      {
	my $name=$p{$a}{name};
	$p{$name}{seq}=$p{$a}{seq};
	$p{$name}{index}=$a;
      }
    return %p;
  }
sub ncbi_tblastx_xml2lib_file 
  {
    my  ($outlib,$string)=(@_);
    my ($L,$l, $a,$b,$c,$d,$i,$nhits,@identifyerL);
    my (%ITERATION);
      
    open (F, ">>$outlib");
    
    $seq=~s/[^a-zA-Z]//g;
    $L=length ($seq);
    
    %ITERATION=xml2tag_list ($string, "Iteration");
    for ($i=0; $i<$ITERATION{n};$i++)
      {
	my ($qindex, $qlen, %hit, $string);
	$string=$ITERATION{$i}{body};

	$qindex=xmltag2value($string,"Iteration_iter-num");
	$qlen  =xmltag2value($string,"Iteration_query-len");
	%hit=&xml2tag_list  ($string, "Hit");

	for ($a=0; $a<$hit{n}; $a++)
	  {
	    my ($string);
	    $string=$hit{$a}{body};
	 
	    $hindex=xmltag2value($string,"Hit_accession")+1;
	    if ($hindex<=$qindex){next;}
	    else  {print F  "# $qindex $hindex\n";}
		   
	   
	    $hlen=xmltag2value  ($string,"Hit_len");
	    %HSP=&xml2tag_list  ($string, "Hsp");
	   
	    for ($b=0; $b<$HSP{n}; $b++)
	      {
		my ($string, $qs,$qe,$qf,$hs,$he,$hf,$s, $d, $e);
		$string=$HSP{$b}{body};
	
		$qs=xmltag2value  ($string,"Hsp_query-from");
		$qe=xmltag2value  ($string,"Hsp_query-to");
		$qf=xmltag2value  ($string,"Hsp_query-frame");

		$hs=xmltag2value  ($string,"Hsp_hit-from");
		$he=xmltag2value  ($string,"Hsp_hit-to");
		$hf=xmltag2value  ($string,"Hsp_hit-frame");
		
		$s=xmltag2value  ($string,"Hsp_identity");
		$l=xmltag2value  ($string,"Hsp_align-len");
		$s=int(($s*100)/$l);
		
		if ($qf>0)
		  {$rqs=$qs; $rqe=$qe;}
		else
		  {
		    $rqe=($qlen-$qs)+1;
		    $rqs=($qlen-$qe)+1;
		  }
		
		if ($hf>0)
		  {$rhs=$hs; $rhe=$he;}
		else
		  {
		    $rhe=($hlen-$hs)+1;
		    $rhs=($hlen-$he)+1;
		  }
		for ($d=0,$e=$rqs; $e<$rqe; $e++,$d++)
		  {
		    my ($r1,$r2);
		    $r1=$e;
		    $r2=$rhs+$d;
		    print F " $r1 $r2 $s 0\n";
		  }
	      }
	  }
      }
    print F "! SEQ_1_TO_N\n";
    
    close (F);
    return %lib;
  }

sub ncbi_tblastpx_xml2lib_file 
  {
    my  ($outlib,$string,%s)=(@_);
    my ($L,$l, $a,$b,$c,$d,$i,$nhits,@identifyerL);
    my (%ITERATION,%hdes, %qdes);
      
    open (F, ">>$outlib");
    
    $seq=~s/[^a-zA-Z]//g;
    $L=length ($seq);
    
    %ITERATION=xml2tag_list ($string, "Iteration");
    for ($i=0; $i<$ITERATION{n};$i++)
      {
	my ($qindex, $qlen, %hit, $string);
	$string=$ITERATION{$i}{body};

	$qdef=xmltag2value($string,"Iteration_query-def");
	%qdes=&tblastpx_name2description($qdef,%s);
	$qlen  =xmltag2value($string,"Iteration_query-len");
	%hit=&xml2tag_list  ($string, "Hit");

	for ($a=0; $a<$hit{n}; $a++)
	  {
	    my ($string);
	    $string=$hit{$a}{body};
	    $hdef=xmltag2value($string,"Hit_def");
	    %hdes=&tblastpx_name2description($hdef,%s);
	    if ($hdes{index}<=$qdes{index}){next;}
	    else  {print F  "# $qdes{index} $hdes{index}\n";}
		   
	   
	    $hlen=xmltag2value  ($string,"Hit_len");
	    %HSP=&xml2tag_list  ($string, "Hsp");
	   
	    for ($b=0; $b<$HSP{n}; $b++)
	      {
		my ($string, $l,$qs,$qe,$qf,$hs,$he,$hf,$s, $d, $e, @s1, @s2);
		$string=$HSP{$b}{body};
	
		$qs=xmltag2value  ($string,"Hsp_query-from");
		$qe=xmltag2value  ($string,"Hsp_query-to");
		$qf=$qdes{frame};
		$qseq=xmltag2value  ($string,"Hsp_qseq");
		
		$hs=xmltag2value  ($string,"Hsp_hit-from");
		$he=xmltag2value  ($string,"Hsp_hit-to");
		$hf=$hdes{frame};
		$hseq=xmltag2value  ($string,"Hsp_hseq");
		
		$s=xmltag2value  ($string,"Hsp_identity");
		$l=xmltag2value  ($string,"Hsp_align-len");
		$s=int(($s*100)/$l);
		@s1=tblastpx_hsp2coordinates($qseq,$qs,$qe,%qdes);
		@s2=tblastpx_hsp2coordinates($hseq,$hs,$he,%hdes);
		
		
		for ($f=0; $f<=$#s1; $f++)
		  {
		    if ($s1[$f]==-1 || $s2[$f]==-1){next;}
		    else 
		      {
			print F " $s1[$f] $s2[$f] $s 0\n";
		      }
		  }
	      }
	  }
      }
    print F "! SEQ_1_TO_N\n";
    
    close (F);
    return %lib;
  }
sub tblastpx_hsp2coordinates
  {
    my ($seq, $s, $e, %des)=@_;
    my @list;
    my @sa;
    my @gap=(-1,-1,-1);
    
    $s=$des{start}+3*($s-1);
  
    if ($des{strand} eq "d"){;}
    else {$s=($des{length}-$s)+1;}
    
    foreach $c (split (//,$seq))
      {
	if ( $c eq '-'){push (@list,@gap);}
	elsif ($des{strand} eq "d")
	  {
	    push(@list,$s++,$s++,$s++);
	  }
	else
	  {
	    push(@list, $s--,$s--,$s--);
	  }
      }
    return @list;
  }

sub tblastpx_name2description
  {
    my ($name, %s)=@_;
    my @at=split("__", $name);
    my %des;

    $des{name}=$at[0];
    $des{strand}=$at[1];
    
    $des{start}=$at[2];
    $des{end}=$at[3];
    $des{length}=$at[4];
    $des{index}=$s{$at[0]}{order}+1;
    return %des;
  }  
sub ncbi_blast_xml2profile 
  {
    my ($name,$seq,$maxid, $minid, $mincov, $string)=(@_);
    my ($L,$l, $a,$b,$c,$d,$nhits,@identifyerL);
    
    
    $seq=~s/[^a-zA-Z]//g;
    $L=length ($seq);
   
    #This is causing the NCBI parser to fail when Iteration_query-def is missing 
    #%query=&xml2tag_list ($string, "Iteration_query-def");
    #$name=$query{0}{body};
    
    %hit=&xml2tag_list ($string, "Hit");
    
    
    for ($nhits=0,$a=0; $a<$hit{n}; $a++)
      {
	my ($ldb,$id, $identity, $expectation, $start, $end, $coverage, $r);
	my (%ID,%DE,%HSP);
	
	$ldb="";

	%ID=&xml2tag_list ($hit{$a}{body}, "Hit_id");
	$identifyer=$ID{0}{body};
	
	%DE=&xml2tag_list ($hit{$a}{body}, "Hit_def");
	$definition=$DE{0}{body};
	
	%HSP=&xml2tag_list ($hit{$a}{body}, "Hsp");
	for ($b=0; $b<$HSP{n}; $b++)
	  {
	    my (%START,%END,%E,%I,%Q,%M);

	 
	    %START=&xml2tag_list ($HSP{$b}{body}, "Hsp_query-from");
	    %HSTART=&xml2tag_list ($HSP{$b}{body}, "Hsp_hit-from");
	    
	    %LEN=  &xml2tag_list ($HSP{$b}{body}, "Hsp_align-len");
	    %END=  &xml2tag_list ($HSP{$b}{body}, "Hsp_query-to");
	    %HEND=  &xml2tag_list ($HSP{$b}{body}, "Hsp_hit-to");
	    %E=&xml2tag_list     ($HSP{$b}{body}, "Hsp_evalue");
	    %I=&xml2tag_list     ($HSP{$b}{body}, "Hsp_identity");
	    %Q=&xml2tag_list     ($HSP{$b}{body}, "Hsp_qseq");
	    %M=&xml2tag_list     ($HSP{$b}{body}, "Hsp_hseq");
	    
	    for ($e=0; $e<$Q{n}; $e++)

	      {
		$qs=$Q{$e}{body};
		$ms=$M{$e}{body};
		
		$expectation=$E{$e}{body};
		$identity=($LEN{$e}{body}==0)?0:$I{$e}{body}/$LEN{$e}{body}*100;
		$start=$START{$e}{body};
		$end=$END{$e}{body};
		$Hstart=$HSTART{$e}{body};
		$Hend=$HEND{$e}{body};
	
		$coverage=($L)?(($end-$start)*100)/$L:0;
	
		if ($identity>$maxid || $identity<$minid || $coverage<$mincov){next;}
		@lr1=(split (//,$qs));
		@lr2=(split (//,$ms));
		$l=$#lr1+1;
		for ($c=0;$c<$L;$c++){$p[$nhits][$c]="-";}
		for ($d=0,$c=0; $c<$l; $c++)
		  {
		    $r=$lr1[$c];
		    if ( $r=~/[A-Za-z]/)
		      {
			
			$p[$nhits][$d + $start-1]=$lr2[$c];
			$d++;
		      }
		  }
		$Qseq[$nhits]=$qs;
		$Hseq[$nhits]=$ms;
		$QstartL[$nhits]=$start;
		$HstartL[$nhits]=$Hstart;
		$identityL[$nhits]=$identity;
		$endL[$nhits]=$end;
		$definitionL[$nhits]=$definition;
		$identifyerL[$nhits]=$identifyer;
		$comment[$nhits]="$ldb|$identifyer [Eval=$expectation][id=$identity%][start=$Hstart end=$Hend]";
		$nhits++;
	      }
	  }
      }
    
    
    $profile{n}=0;
    $profile{$profile{n}}{name}=$name;
    $profile{$profile{n}}{seq}=$seq;
    $profile {n}++;
    
    for ($a=0; $a<$nhits; $a++)
      {
	$n=$a+1;
	
	$profile{$n}{name}="$name\_$a";
	$profile{$n}{seq}="";
	$profile{$n}{Qseq}=$Qseq[$a];
	$profile{$n}{Hseq}=$Hseq[$a];
	$profile{$n}{Qstart}=$QstartL[$a];
	$profile{$n}{Hstart}=$HstartL[$a];
	$profile{$n}{identity}=$identityL[$a];
	$profile{$n}{definition}=$definitionL[$a];
	$profile{$n}{identifyer}=$identifyerL[$a];
	$profile{$n}{comment}=$comment[$a];

	for ($b=0; $b<$L; $b++)
	  {
	    if ($p[$a][$b])
	      {
		$profile{$n}{seq}.=$p[$a][$b];
	      }
	    else
	      {
		$profile{$n}{seq}.="-";
	      }
	  }
      }
    
    $profile{n}=$nhits+1;
    return %profile;
  }
sub ebi_blast_xml2profile 
  {
    my ($name,$seq,$maxid, $minid, $mincov, $string)=(@_);
    my ($L,$l, $a,$b,$c,$d,$nhits,@identifyerL,$identifyer);
    

    
    $seq=~s/[^a-zA-Z]//g;
    $L=length ($seq);
    %hit=&xml2tag_list ($string, "hit");
    
    for ($nhits=0,$a=0; $a<$hit{n}; $a++)
      {
	my ($ldb,$id, $identity, $expectation, $start, $end, $coverage, $r);
	my (%Q,%M,%E,%I);
	
	$ldb=&tag2value ($hit{$a}{open}, "database");
	$identifyer=&tag2value ($hit{$a}{open}, "id");

	$description=&tag2value ($hit{$a}{open}, "description");
	
	%Q=&xml2tag_list ($hit{$a}{body}, "querySeq");
	%M=&xml2tag_list ($hit{$a}{body}, "matchSeq");
	%E=&xml2tag_list ($hit{$a}{body}, "expectation");
	%I=&xml2tag_list ($hit{$a}{body}, "identity");
	

	for ($b=0; $b<$Q{n}; $b++)
	  {

	    $qs=$Q{$b}{body};
	    $ms=$M{$b}{body};
	    
	    $expectation=$E{$b}{body};
	    $identity=$I{$b}{body};
	    
	    	    
	    $start=&tag2value ($Q{$b}{open}, "start");
	    $end=&tag2value ($Q{$b}{open}, "end");
	    $startM=&tag2value ($M{$b}{open}, "start");
	    $endM=&tag2value ($M{$b}{open}, "end");
	    $coverage=(($end-$start)*100)/$L;
	    
	   # print "$id: ID: $identity COV: $coverage [$start $end]\n";
	    
	    if ($identity>$maxid || $identity<$minid || $coverage<$mincov){next;}
	    # print "KEEP\n";

	    
	    @lr1=(split (//,$qs));
	    @lr2=(split (//,$ms));
	    $l=$#lr1+1;
	    for ($c=0;$c<$L;$c++){$p[$nhits][$c]="-";}
	    for ($d=0,$c=0; $c<$l; $c++)
	      {
		$r=$lr1[$c];
		if ( $r=~/[A-Za-z]/)
		  {
		    
		    $p[$nhits][$d + $start-1]=$lr2[$c];
		    $d++;
		  }
	      }
	  
	    $Qseq[$nhits]=$qs;
	    $Hseq[$nhits]=$ms;
	    $QstartL[$nhits]=$start;
	    $HstartL[$nhits]=$Hstart;
	    $identityL[$nhits]=$identity;
	    $endL[$nhits]=$end;
	    $definitionL[$nhits]=$definition;
	    $identifyerL[$nhits]=$identifyer;
	    $comment[$nhits]="$ldb|$identifyer [Eval=$expectation][id=$identity%][start=$startM end=$endM]";
	    $nhits++;
	  }
      }
    
    $profile{n}=0;
    $profile{$profile{n}}{name}=$name;
    $profile{$profile{n}}{seq}=$seq;
    $profile {n}++;
    
    for ($a=0; $a<$nhits; $a++)
      {
	$n=$a+1;
	$profile{$n}{name}="$name\_$a";
	$profile{$n}{seq}="";
	$profile{$n}{Qseq}=$Qseq[$a];
	$profile{$n}{Hseq}=$Hseq[$a];
	$profile{$n}{Qstart}=$QstartL[$a];
	$profile{$n}{Hstart}=$HstartL[$a];
	$profile{$n}{identity}=$identityL[$a];
	$profile{$n}{definition}=$definitionL[$a];	
	$profile{$n}{identifyer}=$identifyerL[$a];
	$profile{$n}{comment}=$comment[$a];

	for ($b=0; $b<$L; $b++)
	  {
	    if ($p[$a][$b])
	      {
		$profile{$n}{seq}.=$p[$a][$b];
	      }
	    else
	      {
		$profile{$n}{seq}.="-";
	      }
	  }
      }
    $profile{n}=$nhits+1;
    
    return %profile;
  }
sub output_profile
  {
    my ($outfile,$profileR, $trim)=(@_);
    my ($a);
    my %profile=%$profileR;
    my $P= new FileHandle;
    my $tmp=vtmpnam();
    
    open ($P, ">$tmp");
    for ($a=0; $a<$profile{n}; $a++)
      {
	print $P ">$profile{$a}{name} $profile{$a}{comment}\n$profile{$a}{seq}\n";
      }
    close ($P);

    if ( $trim)
      {
	&safe_system ("t_coffee -other_pg seq_reformat -in $tmp -action +trim _aln_%%$trim\_K1 -output fasta_aln -out $outfile");
      }
    else
      {
	&safe_system ("mv $tmp $outfile");
      }
    return;
  }
sub blast_xml2hit_list
  {
    my $string=(@_[0]);
    return &xml2tag_list ($string, "hit");
  }
sub xmltag2value
  {
    my ($string_in, $tag)=@_;
    my %TAG;
    %TAG=xml2tag_list ($string_in, $tag);
    return $TAG{0}{body};
  }
      
sub xml2tag_list  
  {
    my ($string_in,$tag)=@_;
    my $tag_in, $tag_out;
    my %tag;
    
    if (-e $string_in)
      {
	$string=&file2string ($string_in);
      }
    else
      {
	$string=$string_in;
      }
    $tag_in1="<$tag ";
    $tag_in2="<$tag>";
    $tag_out="/$tag>";
    $string=~s/>/>##1/g;
    $string=~s/</##2</g;
    $string=~s/##1/<#/g;
    $string=~s/##2/#>/g;
    @l=($string=~/(\<[^>]+\>)/g);
    $tag{n}=0;
    $in=0;$n=-1;
  
 

    foreach $t (@l)
      {

	$t=~s/<#//;
	$t=~s/#>//;
	
	if ( $t=~/$tag_in1/ || $t=~/$tag_in2/)
	  {
	 
	    $in=1;
	    $tag{$tag{n}}{open}=$t;
	    $n++;
	    
	  }
	elsif ($t=~/$tag_out/)
	  {
	    

	    $tag{$tag{n}}{close}=$t;
	    $tag{n}++;
	    $in=0;
	  }
	elsif ($in)
	  {
	   
	    $tag{$tag{n}}{body}.=$t;
	  }
      }
  
    return %tag;
  }


sub seq2gor_prediction 
  {
    my ($name, $seq,$infile, $outfile, $gor_seq, $gor_obs)=(@_);
    my ($l);
    
    `gorIV -prd $infile -seq $gor_seq -obs $gor_obs > gor_tmp`;
    open (GR, ">$outfile");
    open (OG, "gor_tmp");

    while (<OG>)
      {
	
	$l=$_;
	if ($l=~/\>/){print GR "$l";}
	elsif ( $l=~/Predicted Sec. Struct./)
	  {
	    $l=~s/Predicted Sec. Struct\.//;
	    print GR "$l";
	  }
      }
    close (GR);
    close (OG);
    return;
  }
sub seq2msa_tm_prediction 
  {
    my ($name, $seq, $db, $infile, $outfile, $arch, $psv)=(@_);
    my (%p,%gseq,%R, $blast_output, %s, $l);
    my $R2=new FileHandle;
    my $db="uniprot";
    my $method="psitm";
    my $SERVER="EBI";
    
    $blast_output=&run_blast ($name,"blastp", $db, $infile, "outfile");
    
    if (&cache_file("GET",$infile,$name,$method,$db,$outfile,$SERVER))
      {
	print "\tPSITM: USE Cache\n";
	return $outfile;
      }
    else
      {
	$CACHE_STATUS="COMPUTE CACHE";
	%p=blast_xml2profile($name,$seq,$maxid, $minid,$mincov,$blast_output);
	
	
	open (F, ">tm_input");
	for (my $a=0; $a<$p{n}; $a++)
	  {
	    my $s;
	    
	    $s=$p{$a}{seq};
	    $s=uc($s);
	    print F ">$p{$a}{name}\n$s\n";
	    #print stdout ">$p{$a}{name}\n$s\n";
	  }
	close (F);
	print "\tPSITM: kept  $p{n} Homologues for Sequence $p{0}{name}\n";
	&safe_system ("t_coffee -other_pg fasta_seq2hmmtop_fasta.pl -in=tm_input -out=$outfile -output=cons -cov=70 -trim=95 -arch=$arch -psv=$psv");
	unlink ("tm_input");
	&cache_file("SET",$infile,$name,$method,$db,$outfile,$SERVER);
	return;
      }
  }


sub seq2msa_gor_prediction 
  {
    my ($name, $seq,$infile, $outfile, $gor_seq, $gor_obs)=(@_);
    my (%p,%gseq,%R, $blast_output, %s, $l);
    my $R2=new FileHandle;
    my $db="uniprot";
    my $method="psigor";
    my $SERVER="EBI";
    
    $blast_output=&run_blast ($name,"blastp", "uniprot", $infile, "outfile");
    
    if (&cache_file("GET",$infile,$name,$method,$db,$outfile,$SERVER))
      {
	print "\tPSIGOR: USE Cache\n";
	return $outfile;
      }
    else
      {
	$CACHE_STATUS="COMPUTE CACHE";
	%p=blast_xml2profile($name,$seq,$maxid, $minid,$mincov,$blast_output);
	
	
	open (F, ">gor_input");
	for (my $a=0; $a<$p{n}; $a++)
	  {
	    my $s;
	    
	    $s=$p{$a}{seq};
	    $s=uc($s);
	    print F ">$p{$a}{name}\n$s\n";
	    #print stdout ">$p{$a}{name}\n$s\n";
	  }
	close (F);
	print "\tGORTM: kept  $p{n} Homologues for Sequence $p{0}{name}\n";
	&safe_system ("t_coffee -other_pg fasta_seq2hmmtop_fasta.pl -in=gor_input -out=$outfile -output=cons -cov=70 -trim=95 -gor_seq=$gor_seq -gor_obs=$gor_obs -mode=gor");
	unlink ("tm_input");
	&cache_file("SET",$infile,$name,$method,$db,$outfile,$SERVER);
	return;
      }
  }



sub run_blast
  {
    my ($name, $method, $db, $infile, $outfile, $run)=(@_);
    if (!$run){$run=1;}
    my $error_log=vtmpnam();
    
    if (&cache_file("GET",$infile,$name,$method,$db,$outfile,$SERVER) && is_valid_blast_xml ($outfile))
      {return $outfile;}
    else
      {
	$CACHE_STATUS="COMPUTE CACHE";
	if ( $SERVER eq "EBI_SOAP")
	  {
	    &check_configuration ("EMAIL","SOAP::Light","INTERNET");
	    
	    $cl_method=$method;
	    if ($cl_method =~/wu/)
	      {
		$cl_method=~s/wu//;
		if ( $cl_method eq "psiblast")
		  {
		    add_warning($$,$$,"PSI BLAST cannot be used with the wuBLAST Client. Use server=EBI Or server=LOCAL. blastp will be used instead");
		    $cl_method="blastp";
		  }
		
		$command="t_coffee -other_pg wublast.pl --email $EMAIL $infile -D $db -p $cl_method --outfile $outfile -o xml>/dev/null 2>$error_log";
		&safe_system ( $command);
		if (-e "$outfile.xml") {`mv $outfile.xml $outfile`;}
	      }
	    else
	      {
		if ($cl_method eq "psiblast"){$cl_method ="blastp -j5";}
		
		$command="t_coffee -other_pg blastpgp.pl --email $EMAIL $infile -d $db --outfile $outfile -p $cl_method --mode PSI-Blast>/dev/null 2>$error_log";
		&safe_system ( $command);
		
		if (-e "$outfile.xml") {`mv $outfile.xml $outfile`;}
	      }
	  }
	elsif ($SERVER eq "EBI_REST" || $SERVER eq "EBI")
	  {
	    
	    $cl_method=$method;
	    &check_configuration("EMAIL","XML::Simple", "INTERNET");
	    if ($db eq "uniprot"){$db1="uniprotkb";}
	    else {$db1=$db;}
	    

	    if ($cl_method =~/wu/)
	      {
		$cl_method=~s/wu//;
		if ( $cl_method eq "psiblast"){$cl_method="blastp";}
		
		$command="t_coffee -other_pg wublast_lwp.pl --email $EMAIL -D $db1 -p $cl_method --outfile $outfile --align 7 --stype protein $infile>/dev/null 2>error_log";
		
	      }
	    else
	      {
		if ( $cl_method =~/psiblast/){$cl_method ="blastp -j5";}
		$command="t_coffee -other_pg ncbiblast_lwp.pl --email $EMAIL -D $db1 -p $cl_method --outfile $outfile --align 7 --stype protein $infile>/dev/null 2>$error_log";
	      }
	    &safe_system ( $command,5);
	    if (-e "$outfile.out.xml") {`mv $outfile.out.xml $outfile`;}
	    elsif (-e "$outfile.xml.xml"){`mv $outfile.xml.xml $outfile`;}
	    elsif (-e "$outfile.out..xml") {`mv $outfile.out..xml $outfile`;}
	    elsif (-e "$outfile.xml..xml"){`mv $outfile.xml..xml $outfile`;}
	  }
	elsif ($SERVER eq "NCBI")
	  {
	    &check_configuration ("blastcl3","INTERNET");
	    if ($db eq "uniprot"){$cl_db="nr";}
	    else {$cl_db=$db;}
	    
	    if ( $method eq "psiblast")
	      {
		add_warning($$,$$,"PSI BLAST cannot be used with the NCBI BLAST Client. Use server=EBI Or server=LOCAL. blastp will be used instead");
		$cl_method="blastp";
	      }
	    else
	      {
		$cl_method=$method;
	      }
	    $command="blastcl3 -p $cl_method -d $cl_db -i $infile -o $outfile -m 7";
	    &safe_system ($command);
	  }
	elsif ($SERVER =~/CLIENT_(.*)/)
	  {
	    my $client=$1;
	    $command="$client -p $method -d $db -i $infile -o $outfile -m 7";
	    &safe_system ($command);
	  }
	elsif ( $SERVER eq "LOCAL_blastall")
	  {
	    &check_configuration ("blastall");
	    if ($method eq "blastp")
	      {
		$command="blastall -d $db -i $infile -o $outfile -m7 -p blastp";
	      }
	    &safe_system ($command);
	  }
	elsif ( $SERVER eq "LOCAL")
	  {
	    
	    if ($ENV{"BLAST_DB_DIR"})
	      {
		$x=$ENV{"BLAST_DB_DIR"};
		$cl_db="$x$db";
	      }
	    else
	      {
		$cl_db=$db;
	      }
	    
	    if ($method eq "blastp")
	      {
		&check_configuration("blastpgp");
		$command="blastpgp -d $cl_db -i $infile -o $outfile -m7 -j1";
	      }
	    elsif ($method eq "psiblast")
	      {
		&check_configuration("blastpgp");
		$command="blastpgp -d $cl_db -i $infile -o $outfile -m7 -j5";
	      }
	    elsif ($method eq "blastn")
	      {
		&check_configuration("blastall");
		$command="blastall -p blastn -d $cl_db -i $infile -o $outfile -m7 -W6";
	      }	
	    &safe_system ($command);
	  }
	else
	  {
	    
	    myexit(add_error (EXIT_FAILURE,$$,$$,getppid(), "BLAST_FAILURE::UnknownServer",$CL));
	  }
	

	#Check that everything went well
	
	if ( !-e $outfile || !&is_valid_blast_xml($outfile))
	  {
	    
	    if ( -e $outfile)
	      {
		add_warning ($$,$$,"Corrupted Blast Output (Run $run)");
		unlink($outfile);
	      }
	    if ( -e $error_log)
	      {
		
		my $error_msg=file2string ($error_log);
		
		if ( $error_msg =~/enter a valid email/)
		  {
		    myexit(add_error (EXIT_FAILURE,$$,$$,getppid(), "BLAST_FAILURE::Invalid_or_rejected_email::$EMAIL", "$command"));  
		  }
	      }
	    if ( $run==$BLAST_MAX_NRUNS)
	      {
	
		myexit(add_error (EXIT_FAILURE,$$,$$,getppid(), "BLAST_FAILURE::UnknownReason", "$command"));
	      }
	    else
	      {
		my $out;
		if ($SERVER eq "NCBI") {$SERVER="EBI"; }
		elsif ($SERVER eq "EBI"){$SERVER="NCBI";}
		add_warning ($$,$$,"Blast for $name failed (Run: $run out of $BLAST_MAX_NRUNS. Use $SERVER)");
		$out=&run_blast ($name, $method, $db,$infile, $outfile, $run+1);
		if ($SERVER eq "NCBI") {$SERVER="EBI"; }
		elsif ($SERVER eq "EBI"){$SERVER="NCBI";}
		return $out;
	      }
	  }
	
	&cache_file("SET",$infile,$name,$method,$db,$outfile,$SERVER);
	#system ("cp $outfile ~/Dropbox/tmp/cedric.out");
	#die;
	return $outfile;
      }
  }

sub cache_file
  {
    my ($cache_mode,$infile,$name,$method,$db, $outfile,$server)=(@_);
    my $cache_file;
    #Protect names so that they can be turned into legal filenames
    $name=&clean_file_name ($name);

    if ($db=~/\//)
      {
	$db=~/([^\/]+)$/;
	$db=$1;
      }
    $cache_file_sh="$name.$method.$db.$server.tmp";
    $cache_file="$CACHE/$name.$method.$db.$server.tmp";
    
    if ($infile ne "")
      {
	$cache_file_infile_sh="$name.$method.$db.$server.infile.tmp";
	$cache_file_infile="$CACHE/$name.$method.$db.$server.infile.tmp";
      }
    
    if ($cache_mode eq "GET")
      {
	if ($CACHE eq "" || $CACHE eq "no" || $CACHE eq "ignore"  || $CACHE eq "local" || $CACHE eq "update"){return 0;}
	elsif ( !-d $CACHE)
	  {
	    print STDERR "ERROR: Cache Dir: $CACHE Does not Exist";
	    return 0;
	  }
	else
	  {
	    if ( -e $cache_file && &fasta_file1_eq_fasta_file2($infile,$cache_file_infile)==1)
	      {
		`cp $cache_file $outfile`;
		$CACHE_STATUS="READ CACHE";
		return 1;
	      }
	  }
      }
    elsif ($cache_mode eq "SET")
      {
	if ($CACHE eq "" || $CACHE eq "no" || $CACHE eq "ignore"  || $CACHE eq "local" || $CACHE eq "update"){return 0;}
	elsif ( !-d $CACHE)
	  {
	    print STDERR "ERROR: Cache Dir: $CACHE Does not Exist";
	    return 0;
	  }
	elsif (-e $outfile)
	  {
	    `cp $outfile $cache_file`;
	    if ($cache_file_infile ne ""){ `cp $infile $cache_file_infile`;}

	    #functions for updating the cache
	    #`t_coffee -other_pg clean_cache.pl -file $cache_file_sh -dir $CACHE`;
	    #`t_coffee -other_pg clean_cache.pl -file $cache_file_infile_sh -dir $CACHE`;
	    return 1;
	  }
      }
    $CACHE_STATUS="COMPUTE CACHE";
    return 0;
  }
sub file1_eq_file2
  {
    my ($f1, $f2)=@_;
    if ( $f1 eq ""){return 1;}
    elsif ( $f2 eq ""){return 1;}
    elsif ( !-e $f1){return 0;}
    elsif ( !-e $f2){return 0;}
    elsif ($f1 eq "" || $f2 eq "" || `diff $f1 $f2` eq ""){return 1;}
    
    return 0;
  }
sub clean_file_name 
  {
    my $name=@_[0];
    
    $name=~s/[^A-Za-z1-9.-]/_/g;
    return $name;
  }
sub url2file
  {
    my ($address, $out)=(@_);
    
    if (&pg_is_installed ("wget"))
	{
	  return &safe_system ("wget $address -O$out >/dev/null 2>/dev/null");
	}
    elsif (&pg_is_installed ("curl"))
      {
	return &safe_system ("curl $address -o$out >/dev/null 2>/dev/null");
      }
    else
      {
	myexit(flus_error("neither curl nor wget are installed. Imnpossible to fectch remote file"));
	exit ($EXIT_FAILURE);
      }
  }
sub fasta_file1_eq_fasta_file2
  {
    my ($f1, $f2)=@_;
    my (%s1, %s2);
    my @names;
    %s1=read_fasta_seq ($f1);
    %s2=read_fasta_seq ($f2);

    @names=(keys (%s1));
    
    foreach $n (keys(%s1))
      {
	if ($s1{$n}{seq} ne $s2{$n}{seq}){return 0;}
      } 
    
    foreach $n (keys(%s2))
      {
	if ($s1{$n}{seq} ne $s2{$n}{seq}){return 0;}
      }
    return 1;
  }
	


sub read_template_file
{
	my $pdb_templates = @_[0];
	open (TEMP, "<$pdb_templates");
	my %temp_h;
	while (<TEMP>)
{
		$line = $_;
 		$line =~/(\S+)\s(\S+)/;
 		$temp_h{$1}= $2;
#   		print  "H: $temp{$1}\n $1 $2";
}
	close(TEMP);
	return %temp_h;
}

sub calc_rna_template
{
	my ($mode, $infile, $pdbfile, $outfile)=@_;
#  	print "$mode $infile $pdbfile $outfile\n";
# 	print "@_\n";
	my %s, %h ;
	my $result;
	my (@profiles);
	&set_temporary_dir ("set",$infile,"seq.pep");
	%s=read_fasta_seq ("seq.pep");
	
	%pdb_template_h = &read_template_file($pdbfile);
	my $pdb_chain;
	open (R, ">result.aln");


	#print stdout "\n";
	foreach $seq (keys(%s))
	{
		if ($pdb_template_h{$seq} eq "")
		{
			next;
		}
		open (F, ">seqfile");
		print (F ">$s{$seq}{name}\n$s{$seq}{seq}\n");
		close (F);
# 		print "$seq";
		$pdb_chain = $pdb_template_h{$seq};
		$lib_name="$s{$seq}{name}.rfold";
		$lib_name=&clean_file_name ($lib_name);
# 		safe_system ("t_coffee -other_pg RNAplfold2tclib.pl -in=seqfile -out=$lib_name");
		
 		safe_system ("secondary_struc.py seqfile $CACHE$pdb_chain  $lib_name");
		
		if ( !-e $lib_name)
		{
		myexit(flush_error("RNAplfold failed to compute the secondary structure of $s{$seq}{name}"));
			myexit ($EXIT_FAILURE);
		}
		else
		{
			print stdout "\tProcess: >$s{$seq}{name} _F_ $lib_name\n";
			print R ">$s{$seq}{name} _F_ $lib_name\n";
		}
		unshift (@profiles, $lib_name);
	}
	close (R);
	&set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile, @profiles);
}



sub seq2rna_pair{
	my ($mode, $pdbfile1, $pdbfile2, $method, $param, $outfile)=@_;
	
	if ($method eq "runsara.py")
	{
		open(TMP,"<$pdbfile1");
		my $count = 0;
		my $line;
		while (<TMP>)
		{
			$line = $_;
			if ($count ==1)
			{
				last;
			}
			$count += 1;
		}

		
		$chain1 = substr($line,length($line)-3,1);

		close TMP;
		open(TMP,"<$pdbfile2");
		my $count = 0;
		while (<TMP>)
		{
			$line = $_;
			if ($count ==1)
			{
				last;
			}
			$count += 1;
		}
		$chain2 = substr($line,length($line)-3,1);
		close TMP;

		$tmp_file=&vtmpnam();
	
		safe_system("runsara.py $pdbfile1 $chain1 $pdbfile2 $chain2 -s -o $tmp_file --limitation 5000 > /dev/null 2> /dev/null") == 0 or die "sara did not work $!\n";
		open(TMP,"<$tmp_file") or die "cannot open the sara tmp file:$!\n";
		open(OUT,">$outfile") or die "cannot open the $outfile file:$!\n";

		my $switch = 0;
		my $seqNum = 0;
		foreach my $line (<TMP>)
		{
			next unless ($line=~/SARAALI/);
			if ($line=~/>/)
			{
				$switch =0;
				print OUT ">seq$seqNum\n";
				$seqNum++;				
			}
			if ($switch < 2){
				$switch++;
				next;
			}
	
			if ($line =~/REMARK\s+SARAALI\s+([^\*]+)\*/)
			{
				my $string = $1;
				print OUT "$string\n";
			}
		}
		close TMP; 
		close OUT;
		unlink($tmp_file);
	}
}

sub seq2tblastx_lib
  {
    my ($mode, $infile, $outfile)=@_;
    my (%s, $method,$nseq);

    $method=$mode;
    &set_temporary_dir ("set",$infile,"infile");
    %s=read_fasta_seq("infile");
    
    
    foreach $seq (keys(%s))
      {
	$slist[$s{$seq}{order}]=$s{$seq}{seq};
	$sname[$s{$seq}{order}]=$s{$seq}{name};
	$slen[$s{$seq}{order}]=length ($s{$seq}{seq});
      }
    $nseq=$#sname+1;
    open (F, ">outfile");
    print F "! TC_LIB_FORMAT_01\n";
    print F "$nseq\n";
    for ($a=0; $a<$nseq;$a++)
      {
	print F "$sname[$a] $slen[$a]  $slist[$a]\n"
      }
    close (F);
    &safe_system ("formatdb -i infile -p F");
    &safe_system ("blastall -p tblastx -i infile -d infile -m 7 -S1>blast.output");
    
    ncbi_tblastx_xml2lib_file ("outfile", file2string ("blast.output"));
    &set_temporary_dir ("unset",$mode, $method, "outfile",$outfile);
    myexit ($EXIT_SUCCESS);
    }
sub seq2tblastpx_lib
  {
    my ($mode, $infile, $outfile)=@_;
    my (%s, $method,$nseq);
    $method=$mode;
    &set_temporary_dir ("set",$infile,"infile");
    %s=read_fasta_seq("infile");
    
    foreach $seq (keys(%s))
      {
	$slist[$s{$seq}{order}]=$s{$seq}{seq};
	$sname[$s{$seq}{order}]=$s{$seq}{name};
	$slen[$s{$seq}{order}]=length ($s{$seq}{seq});
      }
    $nseq=$#sname+1;
    open (F, ">outfile");
    print F "! TC_LIB_FORMAT_01\n";
    print F "$nseq\n";
    for ($a=0; $a<$nseq;$a++)
      {
	print F "$sname[$a] $slen[$a]  $slist[$a]\n"
      }
    close (F);
    &safe_system("t_coffee -other_pg seq_reformat -in infile -output tblastx_db1 > tblastxdb");
    &safe_system ("formatdb -i tblastxdb -p T");
    &safe_system ("blastall -p blastp -i tblastxdb -d tblastxdb -m7 >blast.output");
    ncbi_tblastpx_xml2lib_file ("outfile", file2string ("blast.output"), %s);
    &set_temporary_dir ("unset",$mode, $method, "outfile",$outfile);
    myexit ($EXIT_SUCCESS);
    }


    


######################3 PID LOCK STUFF

sub file2head
      {
	my $file = shift;
	my $size = shift;
	my $f= new FileHandle;
	my $line;
	open ($f,$file);
	read ($f,$line, $size);
	close ($f);
	return $line;
      }
sub file2tail
      {
	my $file = shift;
	my $size = shift;
	my $f= new FileHandle;
	my $line;
	
	open ($f,$file);
	seek ($f,$size*-1, 2);
	read ($f,$line, $size);
	close ($f);
	return $line;
      }


sub vtmpnam
      {
	my $r=rand(100000);
	my $f="file.$r.$$";
	while (-e $f)
	  {
	    $f=vtmpnam();
	  }
	push (@TMPFILE_LIST, $f);
	return $f;
      }

sub myexit
  {
    my $code=@_[0];
    if ($CLEAN_EXIT_STARTED==1){return;}
    else {$CLEAN_EXIT_STARTED=1;}
    ### ONLY BARE EXIT
    exit ($code);
  }
sub set_error_lock
    {
      my $name = shift;
      my $pid=$$;

      
      &lock4tc ($$,"LERROR", "LSET", "$$ -- ERROR: $name $PROGRAM\n");
      return;
    }
sub set_lock
  {
    my $pid=shift;
    my $msg= shift;
    my $p=getppid();
    &lock4tc ($pid,"LLOCK","LRESET","$p$msg\n");
  }
sub unset_lock
   {
     
    my $pid=shift;
    &lock4tc ($pid,"LLOCK","LRELEASE","");
  }
sub shift_lock
  {
    my $from=shift;
    my $to=shift;
    my $from_type=shift;
    my $to_type=shift;
    my $action=shift;
    my $msg;
    
    if (!&lock4tc($from, $from_type, "LCHECK", "")){return 0;}
    $msg=&lock4tc ($from, $from_type, "LREAD", "");
    &lock4tc ($from, $from_type,"LRELEASE", $msg);
    &lock4tc ($to, $to_type, $action, $msg);
    return;
  }
sub isshellpid
  {
    my $p=shift;
    if (!lock4tc ($p, "LLOCK", "LCHECK")){return 0;}
    else
      {
	my $c=lock4tc($p, "LLOCK", "LREAD");
	if ( $c=~/-SHELL-/){return 1;}
      }
    return 0;
  }
sub isrootpid
  {
    if(lock4tc (getppid(), "LLOCK", "LCHECK")){return 0;}
    else {return 1;}
  }
sub lock4tc
	{
	  my ($pid,$type,$action,$value)=@_;
	  my $fname;
	  my $host=hostname;
	  
	  if ($type eq "LLOCK"){$fname="$LOCKDIR/.$pid.$host.lock4tcoffee";}
	  elsif ( $type eq "LERROR"){ $fname="$LOCKDIR/.$pid.$host.error4tcoffee";}
	  elsif ( $type eq "LWARNING"){ $fname="$LOCKDIR/.$pid.$host.warning4tcoffee";}
	  
	  if ($debug_lock)
	    {
	      print STDERR "\n\t---lock4tc(tcg): $action => $fname =>$value (RD: $LOCKDIR)\n";
	    }

	  if    ($action eq "LCHECK") {return -e $fname;}
	  elsif ($action eq "LREAD"){return file2string($fname);}
	  elsif ($action eq "LSET") {return string2file ($value, $fname, ">>");}
	  elsif ($action eq "LRESET") {return string2file ($value, $fname, ">");}
	  elsif ($action eq "LRELEASE") 
	    {
	      if ( $debug_lock)
		{
		  my $g=new FileHandle;
		  open ($g, ">>$fname");
		  print $g "\nDestroyed by $$\n";
		  close ($g);
		  safe_system ("mv $fname $fname.old");
		}
	      else
		{
		  unlink ($fname);
		}
	    }
	  return "";
	}
	
sub file2string
	{
	  my $file=@_[0];
	  my $f=new FileHandle;
	  my $r;
	  open ($f, "$file");
	  while (<$f>){$r.=$_;}
	  close ($f);
	  return $r;
	}
sub string2file 
    {
    my ($s,$file,$mode)=@_;
    my $f=new FileHandle;
    
    open ($f, "$mode$file");
    print $f  "$s";
    close ($f);
  }

BEGIN
    {
      srand;
    
      $SIG{'SIGUP'}='signal_cleanup';
      $SIG{'SIGINT'}='signal_cleanup';
      $SIG{'SIGQUIT'}='signal_cleanup';
      $SIG{'SIGILL'}='signal_cleanup';
      $SIG{'SIGTRAP'}='signal_cleanup';
      $SIG{'SIGABRT'}='signal_cleanup';
      $SIG{'SIGEMT'}='signal_cleanup';
      $SIG{'SIGFPE'}='signal_cleanup';
      
      $SIG{'SIGKILL'}='signal_cleanup';
      $SIG{'SIGPIPE'}='signal_cleanup';
      $SIG{'SIGSTOP'}='signal_cleanup';
      $SIG{'SIGTTIN'}='signal_cleanup';
      $SIG{'SIGXFSZ'}='signal_cleanup';
      $SIG{'SIGINFO'}='signal_cleanup';
      
      $SIG{'SIGBUS'}='signal_cleanup';
      $SIG{'SIGALRM'}='signal_cleanup';
      $SIG{'SIGTSTP'}='signal_cleanup';
      $SIG{'SIGTTOU'}='signal_cleanup';
      $SIG{'SIGVTALRM'}='signal_cleanup';
      $SIG{'SIGUSR1'}='signal_cleanup';


      $SIG{'SIGSEGV'}='signal_cleanup';
      $SIG{'SIGTERM'}='signal_cleanup';
      $SIG{'SIGCONT'}='signal_cleanup';
      $SIG{'SIGIO'}='signal_cleanup';
      $SIG{'SIGPROF'}='signal_cleanup';
      $SIG{'SIGUSR2'}='signal_cleanup';

      $SIG{'SIGSYS'}='signal_cleanup';
      $SIG{'SIGURG'}='signal_cleanup';
      $SIG{'SIGCHLD'}='signal_cleanup';
      $SIG{'SIGXCPU'}='signal_cleanup';
      $SIG{'SIGWINCH'}='signal_cleanup';
      
      $SIG{'INT'}='signal_cleanup';
      $SIG{'TERM'}='signal_cleanup';
      $SIG{'KILL'}='signal_cleanup';
      $SIG{'QUIT'}='signal_cleanup';
      
      our $debug_lock=$ENV{"DEBUG_LOCK"};
      
      
      
      
      foreach my $a (@ARGV){$CL.=" $a";}
      if ( $debug_lock ){print STDERR "\n\n\n********** START PG: $PROGRAM *************\n";}
      if ( $debug_lock ){print STDERR "\n\n\n**********(tcg) LOCKDIR: $LOCKDIR $$ *************\n";}
      if ( $debug_lock ){print STDERR "\n --- $$ -- $CL\n";}
      
	     
      
      
    }
sub flush_error
  {
    my $msg=shift;
    return add_error ($EXIT_FAILURE,$$, $$,getppid(), $msg, $CL);
  }
sub add_error 
  {
    my $code=shift;
    my $rpid=shift;
    my $pid=shift;
    my $ppid=shift;
    my $type=shift;
    my $com=shift;
    
    $ERROR_DONE=1;
    lock4tc ($rpid, "LERROR","LSET","$pid -- ERROR: $type\n");
    lock4tc ($$, "LERROR","LSET", "$pid -- COM: $com\n");
    lock4tc ($$, "LERROR","LSET", "$pid -- STACK: $ppid -> $pid\n");
   
    return $code;
  }
sub add_warning 
  {
    my $rpid=shift;
    my $pid =shift;
    my $command=shift;
    my $msg="$$ -- WARNING: $command\n";
    print STDERR "$msg";
    lock4tc ($$, "LWARNING", "LSET", $msg);
  }

sub signal_cleanup
  {
    print dtderr "\n**** $$ (tcg) was killed\n";
    &cleanup;
    exit ($EXIT_FAILURE);
  }
sub clean_dir
  {
    my $dir=@_[0];
    if ( !-d $dir){return ;}
    elsif (!($dir=~/tmp/)){return ;}#safety check 1
    elsif (($dir=~/\*/)){return ;}#safety check 2
    else
      {
	`rm -rf $dir`;
      }
    return;
  }
sub cleanup
  {
    #print stderr "\n----tc: $$ Kills $PIDCHILD\n";
    #kill (SIGTERM,$PIDCHILD);
    my $p=getppid();
    $CLEAN_EXIT_STARTED=1;
    
    
    
    if (&lock4tc($$,"LERROR", "LCHECK", ""))
      {
	my $ppid=getppid();
	if (!$ERROR_DONE) 
	  {
	    &lock4tc($$,"LERROR", "LSET", "$$ -- STACK: $p -> $$\n");
	    &lock4tc($$,"LERROR", "LSET", "$$ -- COM: $CL\n");
	  }
      }
    my $warning=&lock4tc($$, "LWARNING", "LREAD", "");
    my $error=&lock4tc($$,  "LERROR", "LREAD", "");
    #release error and warning lock if root
    
    if (isrootpid() && ($warning || $error) )
      {
	
	print STDERR "**************** Summary *************\n$error\n$warning\n";

	&lock4tc($$,"LERROR","RELEASE","");
	&lock4tc($$,"LWARNING","RELEASE","");
      } 
    
    
    foreach my $f (@TMPFILE_LIST)
      {
	if (-e $f){unlink ($f);} 
      }
    foreach my $d (@TMPDIR_LIST)
      {
	clean_dir ($d);
      }
    #No More Lock Release
    #&lock4tc($$,"LLOCK","LRELEASE",""); #release lock 

    if ( $debug_lock ){print STDERR "\n\n\n********** END PG: $PROGRAM ($$) *************\n";}
    if ( $debug_lock ){print STDERR "\n\n\n**********(tcg) LOCKDIR: $LOCKDIR $$ *************\n";}
  }
END 
  {
    
    &cleanup();
  }
   
sub blast_com2new_blast_com
    {
      my $com=shift;
      if ($ENV{"NCBI_BLAST_4_TCOFFEE"} eq "OLD"){return $com;}
      elsif (!&pg_is_installed("legacy_blast.pl")){return $com;}
      else 
	{
	  if ($com=~/formatdb/)
	    {
	      $com=~s/formatdb/makeblastdb/;
	      $com=~s/\-i/\-in/;
	      if ($com =~/pF/){$com=~s/\-pF/\-dbtype nucl/;}
	      if ($com =~/p F/){$com=~s/\-p F/\-dbtype nucl/;}
	      $com="$com -logfile /dev/null";
	      return $com;
	    }
	  elsif ($com =~/^blastn/){return $com;}
	  elsif (&is_blast_package($com))
	    {
	      my $path;
	      
	      if ( $ENV{"NCBI_BIN_4_TCOFFEE"}){$path=$ENV{"NCBI_BLAST_4_TCOFFEE"};}
	      else
		{
		  $path=`which legacy_blast.pl`;
		  $path=~s/\/legacy_blast\.pl//;
		  chomp ($path);
		}
	      $path="--path $path";
	      if ( $com=~/\>\>/){$com=~s/\>\>/ $path \>\>/;}
	      elsif ( $com=~/\>/){$com=~s/\>/ $path \>/;}
	      else {$com.=" $path";}
	      $com="legacy_blast.pl $com";
	      
	      return $com;
	    }
	}
    }
sub safe_system 
{
  my $com=shift;
  my $ntry=shift;
  my $ctry=shift;
  my $pid;
  my $status;
  my $ppid=getppid();
  if ($com eq ""){return 1;}
  
  if ( ($com=~/^blast/) ||($com=~/^formatdb/)){$com=&blast_com2new_blast_com($com);} 

  if (($pid = fork ()) < 0){return (-1);}
  if ($pid == 0)
    {
      set_lock($$, " -SHELL- $com (tcg)");
      if( $debug_generic_method ) { printf "~ exec: %s\n", $com; } 
      exec ($com);
      if( $debug_generic_method ) { printf "~ exitcode: %s\n", $?; } 
    }
  else
    {
      lock4tc ($$, "LLOCK", "LSET", "$pid\n");#update parent
      $PIDCHILD=$pid;
    }
  if ($debug_lock){printf STDERR "\n\t .... safe_system (fasta_seq2hmm)  p: $$ c: $pid COM: $com\n";}

  waitpid ($pid,WTERMSIG);

  shift_lock ($pid,$$, "LWARNING","LWARNING", "LSET");

  if ($? == $EXIT_FAILURE || lock4tc($pid, "LERROR", "LCHECK", ""))
    {
      if ($ntry && $ctry <$ntry)
	{

	  add_warning ($$,$$,"$com failed [retry: $ctry out of $ntry]");
	  lock4tc ($pid, "LRELEASE", "LERROR", "");
	  #if ($com=~/EBI/){$com=~s/EBI/NCBI/;}
	  #elsif ($com=~/NCBI/){$com=~s/NCBI/EBI/;}
	  
	  return safe_system ($com, $ntry, ++$ctry);
	}
      elsif ($ntry == -1)
	{
	  if (!shift_lock ($pid, $$, "LERROR", "LWARNING", "LSET"))
	    {
	      add_warning ($$,$$,"$com failed");
	    }
	  else
	    {
	      lock4tc ($pid, "LRELEASE", "LERROR", "");
	    }
	  return $?;}
      else
	{
	  if (!shift_lock ($pid,$$, "LERROR","LERROR", "LSET"))
	    {
	      myexit(add_error ($EXIT_FAILURE,$$,$pid,getppid(), "UNSPECIFIED system", $com));
	    }
	}
    }
  return $?;
}

sub check_configuration 
    {
      my @l=@_;
      my $v;
      foreach my $p (@l)
	{
	  
	  if   ( $p eq "EMAIL")
	    { 
	      if ( !($EMAIL=~/@/))
		{
		add_warning($$,$$,"Could Not Use EMAIL");
		myexit(add_error ($EXIT_FAILURE,$$,$$,getppid(),"EMAIL","$CL"));
	      }
	    }
	  elsif( $p eq "INTERNET")
	    {
	      if ( !&check_internet_connection())
		{
		  myexit(add_error ($EXIT_FAILURE,$$,$$,getppid(),"INTERNET","$CL"));
		}
	    }
	  elsif( $p eq "wget")
	    {
	      if (!&pg_is_installed ("wget") && !&pg_is_installed ("curl"))
		{
		  myexit(add_error ($EXIT_FAILURE,$$,$$,getppid(),"PG_NOT_INSTALLED:wget","$CL"));
		}
	    }
	  elsif( !(&pg_is_installed ($p)))
	    {
	      myexit(add_error ($EXIT_FAILURE,$$,$$,getppid(),"PG_NOT_INSTALLED:$p","$CL"));
	    }
	}
      return 1;
    }

