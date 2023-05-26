$file=$ARGV[0];
$file0=$ARGV[1];
$filein="$file";
$fileout="$file.loc";
$number="$file0\_metacontig\_";
open(A,"$filein");
      open(B,">$fileout");
      while($line=<A>)
      { chomp $line;
       if($line=~/^SUMMARY BY POSITION/){$s=1;}
       if($line=~/^Help on reading this report/){$s=0;}
       if($s>0)
        {  
          @bb=split/\s+/,$line;  $len=@bb;
          if($line=~/^>/){ @ar=split/\s+/,$line; $scaffold=$ar[0];  $scaffold=~s/>//;}
          if($len>7 && $bb[1]>0 )
             {  
              if($start1<0){$start1=1;}
              if($pre eq "=")
                { $start=$bb[$len-6];  
                  $end=$bb[$len-6]+$bb[$len-5]+1;
                  $start1=$start-10000;
                  $end1=$end+10000;
                  if($start1<1){$start1=1;}
                  print B "$number\_",$scaffold,"\t",$start1,"\t",$end1,"\t",$start,"\t",$end,"\t",$bb[$len-1],"\t",$bb[$len-4],"\t",$bb[$len-3],"\t",$bb[$len-2],"\t",$bb[1],"\n";
                }
            else{ $start=$bb[$len-7];
                  $end=$bb[$len-7]+$bb[$len-6]+1;
                  $start1=$start-10000;
                  $end1=$end+10000;
                  if($start1<1){$start1=1;}
                 print B "$number\_",$scaffold,"\t",$start1,"\t",$end1,"\t",$start,"\t",$end,"\t",$bb[$len-1],"\t",$bb[$len-5],"\t",$bb[$len-4],"\t",$bb[$len-3],"\t",$bb[1],"\n";
                }
             }
          @pp=split//,$line;  $pre=$pp[0];  
        }
      }
      close A;   close B; 
