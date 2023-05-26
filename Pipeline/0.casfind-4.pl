$file=$ARGV[4];
$number="$file\_metacontig\_";
open(A,"$ARGV[0]");
while($line=<A>)
{  chomp $line;
   @bb=split/\s+/,$line;  
   $line=~s/\_metacontig\_\_/\t/g;  
   @ar=split/\s+/,$line;
   $good{$bb[4]}=1;    
   $cpscaffold{$ar[1]}=1; #print $ar[1],"\n";
}
close A;

open(B,"$ARGV[1]");
      open(C,">$ARGV[2]");
      open(D,">$ARGV[3]");
      while($line=<B>)
      { chomp $line;
        $line=~s/>/>$number\_/;
        if($line=~/^>/)
         { $ll=$line; $ll=~s/>//;  
           @bb=split/\s+/,$ll; #print $bb[0],"\n";
           print C ">",$bb[0],"\t",$loc{$bb[0]},"\n";
           if(exists($good{$bb[0]})){$s=1; print D ">",$bb[0],"\t",$loc{$bb[0]},"\n";}
           else{$s=0;}
         }
         else{print C $line,"\n";
              if($s>0){print D $line,"\n";}
             } 
      }
      close B; close C; close D;

