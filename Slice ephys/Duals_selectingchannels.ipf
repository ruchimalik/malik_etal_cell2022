#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//written by Ruchi Malik to select individual channel data from dual recordings

Function Deleterows (data, every)

Wave data
variable every

variable i
Variable iteration = floor(Dimsize(data,0)/every)
Variable overhead = Dimsize(data,0)- iteration*every
//for (i=iteration; i>=0; i-=1)
for (i=0; i< iteration; i+=1)
Deletepoints i, every-1, data
endfor
if(overhead)
Deletepoints i, overhead, data
//SetScale/P x 0,(2*fADCSampleInterval/1000),"" data
endif
End

//-----------
function sparse(data, every)
//NVAR fADCSampleInterval
    wave data
    variable every
   
    extract data, data2, mod(p,every)==0
Duplicate/O data2, data; KillWaves data
end

//------------------------
function channelselect(data, every)

wave data
variable every
variable x
x = Dimsize(data,0)
variable newdim = x/2
make /O/N=(newdim) data_new
variable i
variable j
 // arrays ready. set i = 0 for selecting odd values, set i = 1 for selecting even values
 for( i=every, j=0; i<x; i +=2, j++)
   data_new[j] = data[i];

endfor

Duplicate/O data_new, data; //KillWaves data
end
