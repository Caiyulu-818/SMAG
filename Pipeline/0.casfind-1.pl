$filein="$file";
$fileout="$file.3kb.fa";
open(FASTA,"$ARGV[0]");
open(FA,">$ARGV[1]");
$s=0;
while($line=<FASTA>)
{  chomp $line;
  if($line=~/^>/)
   { if($s>0)
     { if(length($seq)>3000){print FA "$pre\n$seq\n";} }
       $pre=$line;  $seq="";
   }
    else{$seq="$seq$line";}
        $s++;
}
if(length($seq)>3000)
{print FA "$pre\n$seq\n";}
close FASTA; close FA;
open(FF,"$ARGV[1]");
open(AAA,">$ARGV[2]");
$n=0;  $ss=0;
while($myline=<FF>)
{  chomp $myline;  $n++;  
   if($myline=~/^>/)
   { 
     if(length($seq)>3000  && $ss>0){print AAA "$head\n$seq\n";   }   
         $head=$myline;  $seq="";   $ss=1;
if($n>10000)
{ $n=0;  
          open(AAA,">$ARGV[2]");
         }
   }
   else{$seq="$seq$myline";}
}
          close AAA;
