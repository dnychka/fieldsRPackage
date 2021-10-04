# check that all files have a header of 20 lines.



newHeader
cat  << EOF > newHeader.txt
#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2021 Colorado School of Mines
# 1500 Illinois St., Golden, CO 80401
# Contact: Douglas Nychka,  douglasnychka@gmail.edu,
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2
##END HEADER
EOF

foreach i ( fields/R/*.R )
 echo $i
 sed -n '20'p $i >> lastLine.txt
end


foreach i ( *.R )
 echo $i
 sed '1,20d' $i > workR
 cat  newHeader.txt workR > $i
 rm workR
end


foreach i ( fields/R/*.R )
 sed  '21'p $i 
end


### man files

sed 's/^/%/' newHeader.txt > newDocHeader.txt


foreach i ( fields/man/*.Rd )
 echo $i
 sed -n '20'p $i >> lastLine.txt
end

foreach i ( *.Rd )
 echo $i
 sed '1,20d' $i > workR
 cat  newDocHeader.txt workR > $i
 rm workR
end

#### test files

rm lastLine.txt
foreach i ( fields/tests/*.R )
 echo $i
 sed -n '16'p $i >> lastLine.txt
end


foreach i ( *.R )
 echo $i
 sed '1,16d' $i > workR
 cat  newHeader.txt workR > $i
 rm workR
end


