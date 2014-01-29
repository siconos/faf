#script col_counter.gp
col_count=1
good_data=1
while (good_data){
   stats "$0" u (valid(col_count))
   if ( STATS_max ){
      col_count = col_count+1
   } else {
      col_count = col_count-1
      good_data = 0
   }
}