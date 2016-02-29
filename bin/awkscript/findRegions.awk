#!/usr/bin/awk -f
BEGIN{ f=0; }
{ 
 if(f==1) { 
   if($3<th){f=0; print $2;} 
 } else { 
   if($3>th){f=1; printf("%s\t%s\t",$1,$2);}  
 } 
}
