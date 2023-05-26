$file=$ARGV[0];
$filein="$file";
$fileout="$file.loc";
open(A,"$filein");
      open(B,">$fileout");
      while($line=<A>)
      { chomp $line;  
        @bb=split/\t/,$line;
        if($bb[2] eq "CDS")
        {  $scaffold="$number\_$bb[0]";
           $ll=$bb[8];  $ll=~s/;/\t/g;  $ll=~s/\_/\t/g;
           @b=split/\t/,$ll; 
           $gene="$scaffold\_$b[1]";
           print B "$scaffold\t$bb[3]\t$bb[4]\t$bb[6]\t$gene\n";
           $loc{$gene}="$scaffold\t$bb[3]\t$bb[4]\t$bb[6]";
        }
      }
      close A; close B;
