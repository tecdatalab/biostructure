#!/usr/bin/perl

if(@ARGV<1){
 print "$0 [OUTput from MAINMAST -Tree or -Grapgh] (option: [Pbject Name])\n";
 print "Generate Script for pymol -u [script]\n";
 exit;
}

my ($file,$name)=@ARGV;

if($name eq ""){
 $name="TRACE";
}

print "set connect_mode=1\n";
print "load $file, $name,0,pdb\n";
my @B=firstfile($file);

my $maxres=0;
foreach my $k (@B){
 if($$k[0] eq "BOND"){
   printf("bond resi %d and %s, resi %d and %s\n",$$k[1],$name,$$k[2],$name);
 }elsif($$k[0] eq "bond"){ #old format
  printf("bond resi %d and %s, resi %d and %s\n",$$k[2],$name,$$k[6],$name);
 }
}

print "show sticks, $name\n";
print "set connect_mode=0\n";

sub firstfile{
my $cnt=0;
my @A;
open(IN,$_[0]) or die;
while(<IN>){
  next if(/^#/);
 chomp;
 my $item;
 @{$item}=split(/[\s\t]+/,$_);
 push @A, $item
}
close(IN);
return @A;
}
