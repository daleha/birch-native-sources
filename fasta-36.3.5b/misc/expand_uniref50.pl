#!/usr/bin/perl -w

## usage - expand_uniref50.pl uref50acc.file > up_fasta.file
#

# (1) take a list of uniref50 accessions and uses uniref50link to get
#     the associated uniprot accessions
# (2) take the uniprot accessions and produce a fasta library file
#     from them

use strict;
use DBI;

my $db = 'uniprot';
my $user = 'wrpguest';
my $password = 'bmg=uva';

my $dbh = DBI->connect("dbi:mysql:host=wrpxdb2.bioch.virginia.edu:$db",
		       $user, $password,
                       { RaiseError => 1, AutoCommit => 1}
                      ) or die $DBI::errstr;

my %up_sth = ( 
    ur50_to_upacc => "SELECT uniprot_acc FROM  uniref50link WHERE uniref50_acc=?",
    upacc_to_seq => "SELECT * FROM  trFull WHERE acc=?",
    );

for my $sth (keys(%up_sth)) {
  $up_sth{$sth} = $dbh->prepare($up_sth{$sth});
}

my %acc_uniq = ();

while (my $line = <>) {
  next if ($line =~ m/^UniRef50_UPI/);	# _UPI accessions are not in sp-trembl
  chomp($line);
  my ($up_acc, $e_val) = split(/\t/,$line);
  processLine($up_acc,$up_sth{ur50_to_upacc});
}

for my $up_acc ( keys %acc_uniq ) {

  $up_sth{upacc_to_seq}->execute($up_acc);
  while (my $row_href = $up_sth{upacc_to_seq}->fetchrow_hashref ) {
    print ">up|". $row_href->{acc} . "|". $row_href->{name} . " (uref50|$acc_uniq{$up_acc}) " .
      $row_href->{description}. "\n";
    print $row_href->{seq} . "\n";
  }
  $up_sth{upacc_to_seq}->finish();
}

$dbh->disconnect();

sub processLine{
  my ($id,$sth)=@_;

  $id=~ s/UniRef50_//;
  my $result = $sth->execute($id);

  while (my ($acc) = $sth->fetchrow_array()) {
    $acc_uniq{$acc} = $id unless $acc_uniq{$acc};
  }
  $sth->finish();
}
