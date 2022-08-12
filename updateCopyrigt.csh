

foreach i (*.R )
 echo $i
 cp $i workFile
 sed  's/Copyright (C) 2021/Copyright (C) 2022/'  < workFile > $i
 rm workFile
end
