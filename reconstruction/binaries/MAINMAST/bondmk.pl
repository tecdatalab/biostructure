#!/usr/bin/perl

if(@ARGV<1){
 print "$0 [OUTput file from MAINMAST] (option: [Pbject Name])\n";
 print "Generate Script for pymol -u [script]\n";
 exit;
}

my ($file,$name)=@ARGV;

if($name eq ""){
 $name="TRACE";
}

my @B=firstfile($file);

my $maxres=0;
foreach my $k (@B){
 if($$k[0] eq "ATOM"){
  if($$k[1]>$maxres){
   $maxres=$$k[1];
  }
 }
}

print "set connect_mode=1\n";
print "load $file, $name,0,pdb\n";

#bond resi   313 and tmp, resi   327 and tmp
for($i=1;$i<$maxres;$i++){
 printf("bond resi %d and %s, resi %d and %s\n",$i,$name,$i+1,$name);
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
