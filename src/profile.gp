infile='profile.dat'
print infile

basename="profile"

term_choice_tikz=0
if (term_choice_tikz == 1)
{
print "term_choice_tikz";
set term tikz standalone monochrome  size 5in,3in font '\small\sf'
extension = '.tex';
set output basename.extension;
};


#call "col_counter.gp" 'profile.dat'
#print col_count   #number of columns is stored in col_count.

col_count=10
set xrange [1:]
set yrange [-.01:1.]
set key right bottom

plot for [i=2:col_count] infile using 1:i t columnhead w l
