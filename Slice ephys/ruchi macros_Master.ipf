#pragma rtGlobals=1		// Use modern global access method.

//#include <Waves Average>
 Menu "Ruchi macros"
 " KillAll/9",  KillAll()
 "change_dir/5", change_dir()
 "remove300",remove300()
 "Close all graphs/3", CloseAllGraphs()
 "extractsignal",extractsignal()
 //"ZAP", ZAP()
 "Averagetraces/O", Averagetraces()
//"Load-traces/9",loadtrace()
 "Plot-traces",plottrace()
  "Adjustbaseline",adjustbaseline()
 // "measureresjump", measureresjmp()
  // "measureresjmpCC", measureresjmpCC()
  "EPSCslope", EPSCslope()
// "EPSPslope", EPSPslope()
   //" AMPA",  AMPA()
 //  "NMDA", NMDA()
   "tau", tau()  
   "line", line()
" NMDA_decay", NMDA_decay()  
//"Threshold", Threshold()
// "PlotIV", PLOTIV()
   "AP_Threshold", AP_Threshold()
   "AP_Halfwidth", AP_Halfwidth()
   "Firing Output", Firing_output()
   "Spike Duration", Spike_Daration()
   "extractCh1/1", extractCh1()
	 "extractCh2/2", extractCh2()
// "AHP amplitude", AHP_amplitude()
   //"AHP duration", AHP_duration()
   //"ADP amplitude", ADP_amplitude()
  // "ADP duration", ADP_duration()
   "measureIh", measureIh()    
   //"Fit_line", Fit_line()
  "spikeaccomodation", spikeaccomodation()
  "measureISI", measureISI()
   "measureVm",measureVm()
  "Number_Spikes",Number_Spikes()
"Number_Spikes_opto",Number_Spikes_opto()  
   //"EPSPsummation", EPSPsummation()
//  "Paired_Pulse_Ratio", Paired_Pulse_Ratio()
   //" AHP_area",  AHP_area()
   //"averagethetraces" ,averagethetraces()  
   // "EPSP_Latency",  EPSP_Latency()
  // "AMPA_peak",  AMPA_peak()
  // "AMPA_area", AMPA_area()
  //   "NMDA_peak", NMDA_peak()
 //  "NMDA_area", NMDA_area()
 //" NMDA_MK801",   NMDA_MK801()
  // "Currents", Currents()
   //"Variance", Variance()
  // "parab", parab()  
//   "TBS", TBS()
 //  "TBS_5", TBS_5()
  // "EPSPamp", EPSPamp()
  // "line", line()  
 //  "EPSPampfull", EPSPampfull()
//  "E_S_pot", E_S_pot()
 //  "EPSCampfull", EPSCampfull()
  // "Ih_plots",Ih_plots()
   " Zero_Sub/0",  Zero_Sub()
  // "Res_freq", Res_freq()
//   " resonance",  resonance()
  // "measureinputresCC", measureinputresCC()
   //"measVm", measVm()
   
   "inputres_VI", inputres_VI()
   "AP_props", AP_props()
  // "sAHP_mAHP", sAHP_mAHP()
  // " EPSPprops",  EPSPprops()
  // "TBP_area", TBP_area()
   //"backAP_props",backAP_props()
  // "Paired_Pulse_Ratio_slope",Paired_Pulse_Ratio_slope()
 //  "Gradient", Gradient()
  // "PlateauHW", PlateauHW()
  // "PlateauHW_filtered", PlateauHW_filtered()
//   "PPR_slope", PPR_slope()
   " single_EPSC", single_EPSC()
   "single_IPSC", single_IPSC()
   " IPSCtrains",  IPSCtrains()
   " IPSCtimecourse",  IPSCtimecourse()
   "rserVC", rserVC()
   "single_EPSPIPSP/E", single_EPSPIPSP()
   "LRG_ThetaGamma/P", LRG_ThetaGamma()
   "TG_spiking", TG_spiking()
   "normalizePSC", normalizePSC()
   "smooth_waves/S", smooth_waves()
   "runall/R", runall()
 End
 
 function change_dir()
 
NewPath/O Desktop "RuchiHD480:Users:Ruchi:Datafolder:"
execute "SetDataFolder root: "
	NewPath/O/M="Select Data folder" desktop  // set the path variable to the data
	End
	
 //------------------
 Function remove300()
//variable gain = 2
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
variable x  //=gain
prompt x," Enter the gain"
doprompt  "Simply enter thegain", x
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
showinfo /W = $agraph
NVAR slope = root:k1
variable average
average= mean (twave, 0,Pi)

twave= twave /x
endfor


end
//-----------------------
Function doublegraph()

string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 

variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	String name = wavename(agraph,k,1)
	String newname = name+"c2"
	Duplicate twave, $newname

endfor

Display 


end


//----------------
function extractsignal()
NVAR fADCSampleInterval
 string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 

 variable k
for(k=0;k<(numtrace(waves));k = k+1)
String atrace = stringfromlist(k,waves,";")
wave twave = waverefindexed(agraph,k, 1)
//String newwave1 = atrace +  ">" + num2str(k+1)
////String newwave2 = atrace +  ">" +"d2"
//variable x = (numpnts(twave))
//variable y= x/2
//make/O/N= (y)  $newwave1
////make/O/N=(x/2)  $newwave2
//variable i
//
//do 
//newwave1[0, (i-1)]= atrace [2*i]
// 
//while(i<(y))
//		//newwave2[0, ((i/2)-1)]= atrace[2*i+1]
//
 variable i
 for (i=0; i<(numpnts(twave)); i=i+1)
 
 wavestats /Z /R=(0,pi) twave
 
variable  x = V_min +5
 variable y = V_max -20
do   
   
 if(twave[i]< x) //&& (twave[i]> y)
deletepoints i, 1, twave
else
i+=1
endif

while(i<numpnts(twave))
endfor

Smooth/M=0 70, twave
SetScale/P x 0,(2*fADCSampleInterval/1000),"" twave

endfor
 remove300()
 
 end

// if  ((floor(i/2))!=(abs(i/2)) )
//deletepoints i, 1, twave
//else
// i+=1
// endif
// while(i<numpnts(twave))
//---------------------------
function extractCh1()
variable chann
chann = 0
NVAR fADCSampleInterval
 string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 

//prompt chann,"Enter the channel to be exracted"
//doprompt  "channel number 0 or 1", chann


 variable k
// if (chann == 0)
for(k=0;k<(numtrace(waves));k = k+1)
String atrace = stringfromlist(k,waves,";")
wave twave = waverefindexed(agraph,k, 1) 
channelselect(twave,chann)
SetScale/P x 0,(2*fADCSampleInterval/1000),"" twave
endfor


//else 
//for(k=0;k<(numtrace(waves));k = k+1)
// atrace = stringfromlist(k,waves,";")
//wave twave = waverefindexed(agraph,k, 1) 
//channelselect(twave,1)
//SetScale/P x 0,(2*fADCSampleInterval/1000),"" twave
//endfor
//endif
end
 
 
 //---------------
 function extractCh2()
variable chann
chann = 1
NVAR fADCSampleInterval
 string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 

//prompt chann,"Enter the channel to be exracted"
//doprompt  "channel number 0 or 1", chann


 variable k
// if (chann == 0)
for(k=0;k<(numtrace(waves));k = k+1)
String atrace = stringfromlist(k,waves,";")
wave twave = waverefindexed(agraph,k, 1) 
channelselect(twave,chann)
SetScale/P x 0,(2*fADCSampleInterval/1000),"" twave
endfor


//else 
//for(k=0;k<(numtrace(waves));k = k+1)
// atrace = stringfromlist(k,waves,";")
//wave twave = waverefindexed(agraph,k, 1) 
//channelselect(twave,1)
//SetScale/P x 0,(2*fADCSampleInterval/1000),"" twave
//endfor
//endif
end
 
 //----------------
 
 
function samplingrate()
Variable/G S_freq
Variable freq
Prompt freq, " Enter the sampling frequency"
Doprompt  " Enter here",freq
S_freq = freq
end.

Function CloseAllGraphs()
	String name
	do
		name = WinName(0,1)		// name of the front graph
		if (strlen(name) == 0)
			break			// all done
		endif
		DoWindow/K $name // Close the graph
	while(1)
End

 //---------------------------------------------------------------Macros-------------------------------------------------------//
// Function loadtrace()
// newdatafolder /O/S root:raw_data
// Setdatafolder root
// newdatafolder /O/S root:Cells
//
//      newdatafolder /O/S root:cells:newcell
//	Print "Choose file to load"
//	LoadWave/T/P=EP_Data/O
//     //  Print WaveList("*",";","") 
//      // String new=wavelist("*",";","")
//       String newwavename=Stringfromlist(0, new)
//       String newfolder=newwavename[1,15]
//   
//       String newfolder_new = newwavename[1,15] + "new"
//       renamedatafolder root:cells:newcell, $newfolder
//      // newdatafolder /O/S root:raw_data:$newfolder
//       setdatafolder root:cells:$newfolder
//       Duplicatedatafolder root:cells:$newfolder, root:raw_data:$newfolder
//      setdatafolder root:cells:$newfolder
//      variable choice=1
//prompt choice, "Enter 1 if baseline need to load more"
//doprompt  "Your choice is", choice
//if (choice==1)
//loadtrace()
//endif
//      //Duplicatedatafolder root:cells:$newfolder,root:cells:$newfolder_new
//       
//       End






Function KillAll()
	CloseAllGraphs()
	//setdatafolder root:compiled
	KillWaves/A/Z
	//KillVariables/A/Z
	//KillStrings/A/Z
	//setdatafolder root:
End





// // -------------------------------------------------------------------------------------------------------------------------------------------- //    
//       Function plottrace()
//       String list1= WaveList("*",";","") 
//       //String Wave1
//      Variable index1=0
//do
//// Get the next wave name
//String waveA
//WaveA = StringFromList(index1, list1)
//if (strlen(WaveA) == 0)
//break // Ran out of waves
//endif
//if (index1 == 0) // Is this the first wave?
//  String newgraphname=Stringfromlist(0, list1)
//       String newgraph1=newgraphname[9,12]
//       String newgraph2=newgraphname[1,8]
//       String a= a
//       String newgraph= a + newgraph1 + newgraph2
//Display /N=$newgraph, $WaveA
//else
//AppendToGraph $WaveA
//SetAxis/A/N=1 left
//endif
//index1 += 1
//while (1) // Loop until break above
//End
//


//--------------------------------------------------------------------------------------------------------------------------------------------------------------
function numtrace(w)                                                                             //  returns the number of traces in graph
string w
string traces
variable count=0;
variable i = 0;
do
	traces = stringfromlist(i,w)
	count= count+1;
	i = i+1;
while( strlen(traces) != 0)
return (count -1) 


end

//-------------------------------------------------------------------------------------------------------------------------------------------------------
Function UserCursorAdj(grfName)                                                 // generates cursors for user-friendly data selections
String grfName
DoWindow/F $grfName // Bring graph to front
if (V_Flag == 0) // Verify that graph exists
	Abort "UserCursorAdjust: No such graph."
	return -1
endif

NewPanel/K=2 /W=(139,341,382,432) as "Pause for Cursor"
DoWindow/C tmp_PauseforCursor // Set to an unlikely name
AutoPositionWindow/E/M=1/R=$grfName // Put panel near the graph
//Cursor /M /A=1 /K /L=1 /H=1 /C=(0,0,65000) /S=0 A  
//Cursor /M /A=1 /K /L=1 /H=1 /C=(0,0,65000) /S=0 B  
DrawText 21,20,"Adjust the cursors and then"
DrawText 21,40,"Click Continue."
Button button0,pos={80,58},size={92,20},title="Continue"
Button button0,proc=UserCursorAdj_ContButtonProc
PauseForUser tmp_PauseforCursor,$grfName
return 0
End
//--------------------------------------------------------------------------------------------------------------------------------------------------------------
Function UserCursorAdj_ContButtonProc(ctrlName) : ButtonControl
String ctrlName
DoWindow/K tmp_PauseforCursor // Kill self
End
/////_________________
Function Averagetraces()
//remove300()
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1);


			variable i=0
			//string agraph= winname(0,1);
			//SVAR waves = tracenamelist(agraph, ";", 1)
			print "ok"
			Variable numTraces = ItemsInList(waves,";")
			//print waves
			print numTraces, "traces on graph"
			Duplicate/O 	 $(StringFromList(0, waves)) avg_data1
			i=1
			do
	 			WAVE w = $(StringFromList(i, waves))
	  			avg_data1  =  avg_data1 +    w
				i += 1
				//print i
			while ( i <numTraces )
			print "num: ", numTraces
			avg_data1 = avg_data1/numTraces
			//SetDataFolder c1folder
			String tname=	"graph_avg"
			Duplicate/O root:avg_data1	$tname
			
			
			//Display/K=0  root:graph_avg
			//appendtograph  	$tname
				Display/K=0 $tname
			
			//ModifyGraph lsize($tname)=2,rgb($tname)=(0,0,65535)		
				end	
//--------------------------------------------------------------------------------------------------------------------------------------------------------------
Function measureresjmp()
//measures the series resistance jump specified by cursors for both baseline and the jump region,
// Then asks the user for the current -inj 
// Asks if any averages has to be done and then plots a graph
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 

make /O/D/N= (numtrace(waves)) res,Ichange
showinfo /W = $agraph
NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
Print " Place the cursors inoder to measure the resistance jump"
usercursoradj(agraph);
variable k,voltinj = 1E-2
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	Ichange[k] = V_max-V_min;
	print " Maximum of the jump = ",V_max
	Print "Minimum of the jump = ", V_min
	print " The current change = ",Ichange[k]
endfor
prompt voltinj,"Enter the Hyperpolarizing voltage injected"
doprompt  "Voltage(mV)", voltinj 
res = (voltinj/ichange) 


  display /N=ser_res, res


//killwaves slopes 


 setscale /P   x, 0, 0.3, res

 
 Variable avgres= mean(res,0,Pi)
print avgres


make  /O/D/N= (numtrace(waves)) newres 
newres=( res*0) + avgres
print newres
 setscale /P   x, 0, 0.3, newres

Appendtograph /W=ser_res, newres
SetAxis left 10000000,30000000 
Label left "\\Z14\\K(0,15872,65280)Series Resistance"
Label bottom "\\Z14\\K(0,15872,65280)Time (min)"
ModifyGraph nticks=5
ModifyGraph mode(res)=2,lsize(res)=6
ModifyGraph rgb(res)=(0,15872,65280)
ModifyGraph rgb=(0,15872,65280),lstyle(newres)=2
// setscale /P   x, 0, nuslopes;
//Variable avgres1
//avgres1= mean(res,0,Pi)

//make  /O/D/N= (numtrace(waves)) 
//print avgres
end
//killwaves Vchange

//------------------------------------------------------------------------------------------------------------------------------------------------------------
Function adjustbaseline2(x, y)

//new function created for single EPSP IPSP
variable x
variable y
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) slopes
showinfo /W = $agraph
Print " Place the cursors inoder to adjust baseline"
//usercursoradj(agraph);
variable k,firstavg,nextavg,nextavgch
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	//wavestats /Q/Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	wavestats /Q/Z /R=(x,y) twave
	if(k==0)
		firstavg = V_avg;
	else 
		nextavgch = V_avg- firstavg
		twave = twave-nextavgch
	endif
endfor
end





//------------------------------------------------------------------------------------------------------------------------------------------------------------
Function adjustbaseline()
variable x
variable y
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) slopes
showinfo /W = $agraph
Print " Place the cursors inoder to adjust baseline"
usercursoradj(agraph);
variable k,firstavg,nextavg,nextavgch
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	wavestats /Q/Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	//wavestats /Q/Z /R=(x,y) twave
	if(k==0)
		firstavg = V_avg;
	else 
		nextavgch = V_avg- firstavg
		twave = twave-nextavgch
	endif
endfor
end
//-------------------------------------------------------------------------------------------------------------------------
Function EPSCslope()
//for LTP expts only, select 9 traces for pre n 90 for post
// measure the slope in 10-90 % of the rising phase of the EPSP
string agraph= winname(0,1),slopename,atrace

string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) slopes
showinfo /W = $agraph

NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
usercursoradj(agraph);
variable k,x80percentpeak,Ichange,x20percentpeak
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	atrace = stringfromlist(k,waves,";")
	wavestats /Q/Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	Ichange = V_max-V_min;
	//Print " the min voltage  = ",V_min
	//Print " the max voltage =", V_max
	//Print "The voltage chage = ",vchange
	//print "start cursor",xcsr(A,agraph)
	//print "endt cursor",xcsr(B,agraph)
	//x80percentpeak = findepsp(twave,(V_min+(Vchange*0.8)),xcsr(A,agraph),xcsr (B, agraph))
	//x20percentpeak = findepsp(twave,(V_min+(Vchange*0.2)),xcsr(A,agraph),xcsr (B, agraph))
	//print V_max,x90percentpeak
	curvefit  /N/Q   line  twave(xcsr(A,agraph),xcsr (B, agraph)) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	slopes[k] = slope;
endfor
InsertPoints 9,4, slopes
string tempgraph = "temp1"
display /n= $tempgraph slopes

string slopename1
prompt slopename1," Enter the name of the slope wave now"
doprompt  "Simply enter the name madi", slopename1
make /O/D/N= (numtrace(waves)) $slopename1
duplicate /O slopes $slopename1
ModifyGraph mode=2,lsize=6

killwindow  $tempgraph
//killwaves slopes 
display $slopename1 

Variable avg= mean($slopename1,0,8)
print avg    
variable newavg=100/avg


make /O/D/N= (numtrace(waves)) nuslopes
for(k=0;k<(numtrace(waves));k = k+1)


	nuslopes=slopes*newavg 

	endfor
 
 display /n=$slopename1 nuslopes
 setscale /P   x, 0, 0.3, nuslopes; 
 ModifyGraph mode=2,lsize=6
 SetAxis left 0.1,300 
 SetAxis bottom 0,32 

 Label left "\\Z14 % EPSC Slope"
Label bottom "\\Z14 Time (min)"
end.









//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
Function EPSPslope()
//for Current clamp LTP expts only, select 15 traces for pre n 90 for post
// measure the slope in 10-90 % of the rising phase of the EPSP
string agraph= winname(0,1),slopename,atrace

string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) slopes
showinfo /W = $agraph
variable choice=1
prompt choice, "Enter 1 if baseline need to be adjusted"
doprompt  "Your choice is", choice
if (choice==1)
adjustbaseline()
endif

NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
usercursoradj(agraph);
variable k,x80percentpeak,Vchange,x20percentpeak
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	atrace = stringfromlist(k,waves,";")
	wavestats /Q/Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	Vchange = V_max-V_min;
	//Print " the min voltage  = ",V_min
	//Print " the max voltage =", V_max
	//Print "The voltage chage = ",vchange
	//print "start cursor",xcsr(A,agraph)
	//print "endt cursor",xcsr(B,agraph)
	//x80percentpeak = findepsp(twave,(V_min+(Vchange*0.8)),xcsr(A,agraph),xcsr (B, agraph))
	//x20percentpeak = findepsp(twave,(V_min+(Vchange*0.2)),xcsr(A,agraph),xcsr (B, agraph))
	//print V_max,x90percentpeak
	curvefit  /N/Q   line  twave(xcsr(A,agraph),xcsr (B, agraph)) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	slopes[k] = slope;
endfor
//insert a single point showing the 100 hz LTP data point
variable base
prompt base,"Enter the number of traces in baseline"
doprompt "enter the name madi", base
InsertPoints base,3, slopes
string tempgraph = "temp1"
display /n= $tempgraph slopes

string slopename1
prompt slopename1," Enter the name of the slope wave now"
doprompt  "Simply enter the name madi", slopename1
make /O/D/N= (numtrace(waves)) $slopename1
duplicate /O slopes $slopename1
ModifyGraph mode=2,lsize=6

killwindow  $tempgraph
//killwaves slopes 
display $slopename1 
//calculate baseline mean
variable trace
trace= base-1
Variable avg= mean($slopename1,0,trace)
print avg    
variable newavg=100/avg // calculate % average of baseline

//make new wave by dividing all original slope values by the %average of the baseline
make /O/D/N= (numtrace(waves)) nuslopes
for(k=0;k<(numtrace(waves));k = k+1)


	nuslopes=slopes*newavg 

	endfor
 //set axis, labels etc
 display /n=$slopename1 nuslopes
 setscale /P   x, 0, 0.3, nuslopes; 
 ModifyGraph mode=2,lsize=6
 SetAxis left 0.1,300 
 SetAxis bottom 0,32 

 Label left "\\Z14 % EPSP Slope"
Label bottom "\\Z14 Time (min)"
end.
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



//measure series resistance in current clamp mode
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Function measureresjmpCC()
//measures the series resistance jump specified by cursors for both baseline and the jump region,
// Then asks the user for the current -inj 
// Asks if any averages has to be done and then plots a graph
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) res,Vchange
showinfo /W = $agraph
variable choice=1
prompt choice, "Enter 1 if baseline need to be adjusted"
doprompt  "Your choice is", choice
if (choice==1)
adjustbaseline()
endif
NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
Print " Place the cursors in oder to measure the resistance jump"
usercursoradj(agraph);
variable k,currinj = 1E-12
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	Vchange[k] = V_max-V_min;
	print " Maximum of the jump = ",V_max
	Print "Minimum of the jump = ", V_min
	print " The voltage change = ",Vchange[k]
endfor
prompt CurrInj,"Enter the Hyperpolarizing current injected"
doprompt  "Current(pA)", CurrInj 
res = (Vchange/CurrInj) 
display res


 Variable avgres=(mean(res,0,Pi))
print avgres
Variable avg
avg=100/avgres
print avg
make /O/D/N=(numtrace (waves)) newres

for(k=0;k<(numtrace(waves));k = k+1)


	newres=res*avg

	endfor

display res
 setscale /P   x, 0, 0.3, res

Appendtograph /W=ser_res, res
SetAxis left 5e6,35e6
Label left "\\Z14\\K(0,15872,65280) Series Resistance"
Label bottom "\\Z14\\K(0,15872,65280)Time (min)"
ModifyGraph nticks=5
ModifyGraph mode(res)=2,lsize(res)=6
ModifyGraph rgb(res)=(0,15872,65280)
ModifyGraph rgb=(0,15872,65280),lstyle(res)=2


SetAxis bottom 0,32 
// setscale /P   x, 0, nuslopes;
//Variable avgres1
//avgres1= mean(res,0,Pi)

//make  /O/D/N= (numtrace(waves)) 
//print avgres
end
//killwaves Vchange

 


//-----------------------------------------------------------------------------------------------------------------------------------------------//

Function AMPA()
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) peak
showinfo /W = $agraph
NVAR slope = root:k1
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Peak[k] = V_min
endfor
Variable Average_AMPA
Average_AMPA= mean(peak,0,Pi)
print average_AMPA

display /N=AMPA, peak

end
//---------------------------------------------------------------------------------------
Function NMDA()
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) peak
showinfo /W = $agraph
NVAR slope = root:k1
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Peak[k] = V_max
endfor
Variable Average_NMDA
Average_NMDA= mean(peak,0,Pi)
print average_NMDA

display /N=NMDA, peak

end

//------------------------------------------------------------

Function sAHP_mAHP()

string agraph= winname(0,1),atrace
string waves = tracenamelist(agraph,";",1); 
showinfo /W = $agraph
//zero_sub()
//usercursoradj(agraph);
wave twave = waverefindexed(agraph,0, 1)
Wavestats /Z/R=(0,60) twave
variable b= V_avg
//usercursoradj(agraph);
wave twave = waverefindexed(agraph,0, 1)
Wavestats /Z/R=(900,1100) twave

Variable K
K= V_min

wave twave = waverefindexed(agraph,0, 1)
Wavestats /Z/R=(1100,1200) twave

Variable m
m= V_min



variable sl= (-1)*(k-b)
variable med= (-1) * (m-b)

Printf  "The AHP of this cell is %g ", sl
Printf  "The AHP of this cell is %g ", med

END



//-------------------------------------------------------------------------------------------------------------
Function AHP_amplitude()

string agraph= winname(0,1),atrace
string waves = tracenamelist(agraph,";",1); 
showinfo /W = $agraph
//zero_sub()
usercursoradj(agraph);
wave twave = waverefindexed(agraph,0, 1)
Wavestats /Z/R=(0,40) twave
variable m= V_avg
usercursoradj(agraph);
wave twave = waverefindexed(agraph,0, 1)
Wavestats /Z/R=(100,300) twave

Variable K
K= V_min

variable n= (-1)*(k-m)

Printf  "The AHP of this cell is %g ", n
END

//-------------------------------------------------------------------------------------------------------------
Function ADP_amplitude()

string agraph= winname(0,1),atrace
string waves = tracenamelist(agraph,";",1); 
showinfo /W = $agraph
usercursoradj(agraph);
wave twave = waverefindexed(agraph,0, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Variable K
K= V_max*(-1000)
Printf  "The ADP of this cell is %g ", K
END

//---------------------------------------------------------------------------------------------------------------------

Function PlotIV()

string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) Vchange
showinfo /W = $agraph
NVAR slope = root:k1
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
 
 if  (V_min*(-1)>V_max)
 
 Vchange[K]=V_min
 elseif (V_min*(-1)<V_max)
 Vchange[K]=V_max
 endif
endfor
make /O/D/N=(numtrace(waves))current
make /O/D/N=(numtrace(waves))new
variable x
variable y

for (x=0;x<(numtrace(waves));x = x+1)
y= -0.0000000001
current[x] = -50*x
new[x] = (current[x] + y)
endfor
string voltage
Display /n=voltage Vchange vs new
ModifyGraph mode=3,marker=8
make /O/D/N= 11 Vchange_new, Vchange_final

for(k=0;k<(numtrace(waves));k = k+1)
Vchange_new[k]= (Vchange[K]) +(Vchange [K+11])+(Vchange[k+22])
Vchange_final [k]= (Vchange_new[k])/3
endfor

variable V
V = 0.1
make /O/D/N =11 vol
variable g
for (g=0;g<11;g=g+1)
vol [g]= (-0.5)+ (g*V)
endfor

Display Vchange_final vs vol
CurveFit line  Vchange_final /X=vol /D


//curvefit  /N /Q line  Vchange(,6) /D

ModifyGraph mode(Vchange_final)=3,marker(Vchange_final)=19;DelayUpdate
ModifyGraph lsize(fit_Vchange_final)=2,rgb(fit_Vchange_final)=(4352,4352,4352)
ModifyGraph lstyle(fit_Vchange_final)=7
ModifyGraph zero=1
End
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Function inputres_VI()

//remove300()

string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) Vchange, V_initial, V_stable, new_Vchange
showinfo /W = $agraph
NVAR slope = root:k1
//usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
//Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph))twave
wavestats /Z /R=(20,100) twave
 V_initial[k]= V_avg
 endfor
 
 //usercursoradj(agraph);
variable l
for(l=0;l<(numtrace(waves));l = l+1)
	wave twave = waverefindexed(agraph,l, 1)
//Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph))twave
wavestats /Z /R=(760,890) twave
 V_stable[l]= V_avg
endfor

variable m
for(m=0;m<(numtrace(waves));m = m+1)
	wave twave = waverefindexed(agraph,m, 1)

Vchange[m]= V_initial[m]-V_stable[m]
endfor

for(m=0;m<(numtrace(waves)); m=m+1)
	wave twave = waverefindexed(agraph,m, 1)

new_Vchange[m]= (-1)* Vchange[m]
endfor

variable V
V = 10
string num
prompt num," Enter the # of traces"
doprompt  "Simply enter the number", num
make /O/D/N= (str2num(num)) vol
variable g
for (g=0;g<str2num(num);g=g+1)
vol [g]= (-1)*((50)- (g*V))
endfor

Display new_Vchange vs vol
CurveFit/X=1/NTHR=0 line  new_Vchange /X=vol /D 
ModifyGraph mode(new_Vchange)=3,rgb(fit_new_Vchange)=(4369,4369,4369)
ModifyGraph marker(new_Vchange)=19


//curvefit  /N /Q line  Vchange(,6) /D

ModifyGraph mode(new_Vchange)=3,marker(new_Vchange)=19;DelayUpdate
ModifyGraph lsize(fit_new_Vchange)=2,rgb(fit_new_Vchange)=(4352,4352,4352)
ModifyGraph lstyle(fit_new_Vchange)=7
ModifyGraph zero=1

End

//----------------------------------------------------------
Function Fit_line()
string newgraph
usercursoradj(newgraph);

	curvefit  /N /Q line  Vchange(xcsr (A, newgraph), xcsr (B, newgraph)) /D
	
End
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
Function countspikefrequency(trace)
// This function will read a single trace and then returns spike frequency for each spike in the trace

wave trace
variable choice
variable n = countAPnum (trace)
make /O /D /N= (n) peaktime, peakwidth, peakamp, Intspint
make /O/D/N=1 ratio
make /O/D/N= (n)  firingfreq //, pamp = NaN, pwidth = NaN
variable i = 0 , xpos1 =110, j;


do
	findpeak /Q /B=0.01 /M=5 /R=(xpos1, 920) trace
	//findpeak /Q /B=1 /M=0.02 /R=(xpos1, 0.8) trace
	if (V_flag != 0)  
	
		break
	endif
	
	
	peaktime[i] = V_peakLoc
	peakamp [i] = V_peakval
	peakwidth [i] = V_peakwidth
	i = i+1;
	
	//peakamp[i]= V_peakVal
	//peakwidth[i]=V_peakWidth
	
	xpos1 = V_peakLoc + 1;
while(1)

if (i>1) 
	for (j=0;j<=i-1;j=j+1)
	Intspint[j]= (peaktime[j+1]-peaktime[j])
		firingfreq[j] = 100/(peaktime[j+1]-peaktime[j])
		DeletePoints n,1, firingfreq
		wavestats firingfreq
		//variable avgfiring= V_avg*1000
		
		
		
		
		//pamp[j]= peakamp[i]
		//pwidth[j]=peakwidth[i]
	endfor
else
	firingfreq=0
	print " Sorry there are no firings in this trace" 
//endif
//if (choice==0)

	
	//return (pamp)
	//return(peakwidth)
//else
	//return (i)


return firingfreq [j] 	
endif

//ratio = firingfreq[i]/firingfreq[1]
//killwaves peaktime 
//killwaves firingfreq
//killwaves avgfiring
variable CV_ISI= (sqrt(variance (Intspint)))*100/(mean (Intspint))
variable CV_amp= (sqrt(variance (peakamp)))*100/(mean (peakamp))
variable CV_width= (sqrt(variance (peakwidth)))*100/(mean (peakwidth))
make /o/D/N=(6) CoV
CoV[0] = CV_ISI
CoV[1]= CV_amp
CoV[2]= CV_width
wavestats firingfreq
variable x = V_npnts
print x, "muber of point"
variable num = x-1
variable SFA_f= firingfreq[0]/firingfreq[num]
wavestats peakamp
variable y = V_npnts
variable SFA_a= peakamp[0]/peakamp[y]
wavestats peakwidth
variable z = V_npnts
variable SFA_hw= peakwidth[z]/peakwidth[0]
make /o/D/N=(3) SFA
CoV[3] = SFA_f
CoV[4]=SFA_a
CoV[5]= SFA_hw
end
//------------------------------------------------------------------------------------------------------------------------------------------------------------
function spikeaccomodation()
variable index = 0,first =0//,maxap = 0;	
string agraph= winname(0,1);
//SVAR c1folder = root:c1folder
//setdatafolder c1folder
string atrace, aspikefreq
make /O /D /N= 1 peaks

do
	atrace = wavename(agraph,index,1)
	
	if (waveexists($atrace) == 0)
		break
	endif
	
	aspikefreq = atrace +  ">" + num2str (index+1) 
	peaks =   countspikefrequency($atrace)

	duplicate /O peaks $aspikefreq
	
	
	//if (countspikefrequency($atrace,1) >= maxap)
		//maxap = countspikefrequency($atrace,1)
		//maxname = aspikefreq
	//endif
	
//		if ((peaks[1] != 0) && (first == 0))
//		         first = 1
//			display $aspikefreq
//			print "displaying...",aspikefreq
//			DrawText 0.332157618631268,1.20175222997122,"Spike number"
//			SetDrawEnv textrot= 90
//			DrawText -0.117753623188405,1.13382899628252,"Instantaneous frequency (Hz)"
//		elseif (peaks[1] != 0)
//			appendtograph /W= graph1 $aspikefreq
//			print "adding to display...",aspikefreq
//		endif
	index = index + 1
while(1)
//print maxname
//string ap6 = " MaxAP" + "__" + maxname
//rename $maxname $ap6 
//appendtograph /w = newgraph $maxname
//display /n = newgraph $ap6
//setaxis left 0,120
//setaxis bottom 0,20
//DrawText 0.332157618631268,1.20175222997122,"Spike number"
//SetDrawEnv textrot= 90
//DrawText -0.117753623188405,1.13382899628252,"Instantaneous frequency (Hz)"
//killwaves peaks;
//killwindow graph1
//killwindow graph2
setdatafolder root:
end
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
Function measureVm()
NVAR freq = s_freq
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) Vm
showinfo /W = $agraph
NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
Print " Indicate the window by placing the cursors to calculate average Vm"
//usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	wavestats /Z /R=(1,2300) twave
	Vm[k] = V_avg;
endfor
display Vm
end
//--------------------------------------------------------------------------------------------------------------------------------------------
Function ISI(trace, choice)
// This function will read a single trace and then returns spike latency for each spike in the trace

wave trace
variable choice
variable n = countAP (trace)
make /O /D /N= (n) peaktime, firingfreq,Intspint

variable i = 0 , xpos1 =0.1,,j;
do
	findpeak /Q /B=0.01 /M=5 /R=(xpos1, pnt2x(trace, numpnts(trace))) trace
	//findpeak /Q /B=1 /M=0.02 /R=(xpos1, 0.8) trace
	if (V_flag != 0)  
		break
	endif
	
	peaktime[i] = V_peakLoc 
	i = i+1;
	//peakamp[i]= V_peakVal
	//peakwidth[i]=V_peakWidth
	
	xpos1 = V_peakLoc+1;
while (1)

//if (i>1) 
	for (j=0;j<=i-1;j=j+1)
	Intspint[j]= (peaktime[j+1]-peaktime[j])
		firingfreq[j] = 100/(peaktime[j+1]-peaktime[j])
		wavestats firingfreq
		//variable avgfiring= V_avg*1000
		
		
		
		
		//pamp[j]= peakamp[i]
		//pwidth[j]=peakwidth[i]
	endfor
//else
	//firingfreq=0
//make /O /D /N= (25)  firingfreq = NaN//, pamp = NaN, pwidth = NaN
//if (i>1) 
	//for (j=1;j<=i-1;j=j+1)
		//firingfreq[j] = peaktime[j]
		
	//endfor
//else
	//firingfreq=0
	//print " Sorry there are no firings in this trace"
//endif
//if (choice==0)
	return (peaktime[i])
	//return (pamp)
	//return(peakwidth)
//else
	//return (i)
//endif
variable CV_ISI= (sqrt(variance (Intspint)))/(mean (Intspint,0,pi))
print CV_ISI
killwaves peaktime
killwaves firingfreq
//endif
end
//------------------------------------------------------------------------------------------------------------------------------------------------------------
function measureISI()
variable index = 0,first =0//,maxap = 0;
string agraph= winname(0,1);

//SVAR c1folder = root:c1folder
//setdatafolder c1folder

string atrace, aspikefreq//, maxname
make /O /D /N= (150)  peaks
do
	atrace = wavename(agraph,index,1)
	if (waveexists($atrace) == 0)
		break
	endif
	
	aspikefreq = atrace +  ">" + num2str(index+1) 
	peaks =  ISI($atrace,0)
	duplicate /O peaks $aspikefreq
	//if (countspikefrequency($atrace,1) >= maxap)
		//maxap = countspikefrequency($atrace,1)
		//maxname = aspikefreq
	//endif
	
//		if ((peaks[1] != 0) && (first == 0))
//		         first = 1
//			display $aspikefreq
//			print "displaying...",aspikefreq
//			DrawText 0.332157618631268,1.20175222997122,"Spike number"
//			SetDrawEnv textrot= 90
//			DrawText -0.117753623188405,1.13382899628252,"Instantaneous frequency (Hz)"
//		elseif (peaks[1] != 0)
//			appendtograph /W= graph1 $aspikefreq
//			print "adding to display...",aspikefreq
//		endif
//	index = index +1
while(1)
////print maxname
////string ap6 = " MaxAP" + "__" + maxname
////rename $maxname $ap6 
////appendtograph /w = newgraph $maxname
////display /n = newgraph $ap6
////setaxis left 0,120
////setaxis bottom 0,20
////DrawText 0.332157618631268,1.20175222997122,"Spike number"
////SetDrawEnv textrot= 90
////DrawText -0.117753623188405,1.13382899628252,"Instantaneous frequency (Hz)"
////killwaves peaks;
////killwindow graph1
////killwindow graph2
setdatafolder root:
end



//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
Function measureIh()
NVAR freq = s_freq
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
SVAR c1folder = root:c1folder
//setdatafolder c1folder
//Zero_Sub()
//adjustbaseline()
Zero_Sub()
make /O/D/N= (numtrace(waves)) SSV
make /O/D/N= (numtrace(waves)) sag, Reb, timeofsag, timeofreb, TOP1, TOP2, TOP3, TOP4
showinfo /W = $agraph
NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
Print " Indicate the window by placing the cursors to calculate average Vm"
//usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	//wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	wavestats /Z /R=(140,200) twave
	//Sag[k] = V_min;
	Timeofsag [k]=V_minLoc
TOP1 [k]=Timeofsag [K]-0.0005
TOP2 [k]=Timeofsag [k]+0.0005

Wavestats /Z/R=(TOP1[K],TOP2 [k]) twave
Sag[k] = V_avg


endfor
//usercursoradj(agraph);
variable m
for(m=0;m<(numtrace(waves));m = m+1)
	wave twave = waverefindexed(agraph,m, 1)
	//wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	wavestats /Z /R=(800,890) twave
	SSV[m] = V_avg;
endfor
//usercursoradj(agraph);

variable n
for(n=0;n<(numtrace(waves));n = n+1)
	wave twave = waverefindexed(agraph,n, 1)
	//wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	wavestats /Z /R=(910,1010) twave
	//Reb[n] = V_max;
	Timeofreb [n]=V_maxLoc
TOP3[n]=Timeofreb [n]-0.0005
TOP4 [n]=Timeofreb [n]+0.0005

Wavestats /Z/R=(TOP3[n],TOP4 [n]) twave
Reb[n] = V_avg

endfor
make /O/D/N= (numtrace(waves)) sagR
make /O/D/N= (numtrace(waves)) sagV
make /O/D/N= (numtrace(waves)) Rebound

variable j
for (j=0;j<(numtrace(waves));j = j+1)
SagR[j]=Sag[j]/SSv[j]
SagV[j]=sag[j]-SSV[j]
Rebound[j]= Reb[j]
display rebound vs SSV

CurveFit/X=1/NTHR=0 line  rebound /X=SSV/D 

endfor
Variable Sag_ratio
Sag_ratio= mean(sagR,0,Pi)
//Variable Sag_voltage
//Sag_voltage=mean(sagV,0,Pi)

//ariable V
//V = 50
//make /O/D/N =11 vol
//variable g
//for (g=0;g<11;g=g+1)
//vol [g]= (-1)*((50))- (g*V)
//endfor
//Display SagV vs vol
//SetAxis bottom 0,-500
ModifyGraph mode=4,marker=19,msize=5
//Display Rebound vs vol
//SetAxis bottom 0,-500
//ModifyGraph mode=4,marker=19,msize=5

Print " The average sag ratio for this cell is", Sag_ratio
//Print " The average sag voltage for this cell is", Sag_voltage,"mV"

//setdatafolder root:

end


//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function EPSPsummation()
NVAR freq = s_freq
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
Zero_Sub()
SetAxis bottom 50,350
//adjustbaseline()
make /O/D/N= (numtrace(waves)) SSV,Sag, sumratio

showinfo /W = $agraph
NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
Print " Indicate the window by placing the cursors to calculate average Vm"
usercursoradj(agraph);

variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	//wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	wavestats /Z /R=(60,90) twave
	SSV[k] = V_max;
endfor
usercursoradj(agraph);
variable m
for(m=0;m<(numtrace(waves));m = m+1)
	wave twave = waverefindexed(agraph,m, 1)
	//wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	wavestats /Z /R=(260,290) twave
	
	Sag[m] = V_max;
endfor


variable j
for (j=0;j<(numtrace(waves));j = j+1)
Sumratio[j]=sag[j]/SSV[j]

endfor
Variable Sum_ratio
Sum_ratio= mean(sumratio,0,Pi)



Print " The average sumratio for this cell is", Sum_ratio

end
//-------------------
Function countAPnum(trace)
// This function will read a single trace and then returns number of APS for each spike in the trace
wave trace


make /O /D /N= (200)  peaktime
variable i = 0 , xpos1 =0.01,,j;
do
	findpeak /Q /B=0.0001/M=5/R=(xpos1, pnt2x(trace, numpnts(trace))) trace
	//findpeak /Q /B=1 /M=0.02 /R=(xpos1, 0.7) trace
	if (V_flag != 0)  
		break
	endif
	i = i+1;
	peaktime[i] = V_peakLoc
	
	xpos1 = V_peakLoc+1;
while(1)
variable AP=(i)
variable freq = Ap/0.8
return(AP)
end.


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Function countAP(trace)
// This function will read a single trace and then returns number of APS for each spike in the trace
wave trace


make /O /D /N= (200)  peaktime
variable i = 0 , xpos1 =0.01,,j;
do
	findpeak /Q /B=0.0001/M=5/R=(xpos1, pnt2x(trace, numpnts(trace))) trace
	//findpeak /Q /B=1 /M=0.02 /R=(xpos1, 0.7) trace
	if (V_flag != 0)  
		break
	endif
	i = i+1;
	peaktime[i] = V_peakLoc
	
	xpos1 = V_peakLoc+1;
while(1)
variable AP=(i)
variable freq = Ap/0.8
return(freq)
end.
//------------------------------------------------------------------------------------------------------------------------------------------------------------
Function Number_spikes()
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) freq, APnum

variable k
for(k=0;k<(numtrace(waves));k = k+1)

	wave twave = waverefindexed(agraph,k, 1)
make /O /D /N= (200)  peaktime
variable i = 0 , xpos1 =0.01,,j;
do
	findpeak /Q /B=0.0001/M=5/R=(xpos1, pnt2x(trace, numpnts(trace))) twave
	if (V_flag != 0)  
		break
	endif
	i = i+1;
	peaktime[i] = V_peakLoc
	
	xpos1 = V_peakLoc+1;
while(1)
variable AP=(i)
variable freq1 = Ap/0.8
Freq [k] = freq1
APnum [k] = AP
endfor

make /o/n=20 currentinj
currentinj [0] = 0
variable m 
for(m=1;m<(20);m = m+1)
currentinj [m] = currentinj [m-1] + 25

endfor
string newgraph
display /N=newgraph, freq vs currentinj 
//DoWindow/F newgraph // Bring graph to front
//usercursoradj(newgraph);
//curvefit  /N/Q   line  freq(xcsr (A, newgraph), xcsr (B, newgraph)) /D

//variable index = 0,first =0,maxap = 0;
//string agraph= winname(0,1);
//string atrace, spikenumber
//SVAR c1folder = root:c1folder
//setdatafolder c1folder
//
//make /O /D /N= (1)  peaks
//do
//	atrace = wavename(agraph,index,1)
//	if (waveexists($atrace) == 0)
//		break
//	endif
//	
//	spikenumber = "b" + atrace +  ">" + num2str(index+1) 
//	peaks =  countAP($atrace)
//	duplicate /O peaks $spikenumber
//	
//		if ((peaks[1] != 0) && (first == 0))
//		        first = 1
//			display $spikenumber 
//			print "displaying...",spikenumber
//			
//		elseif (peaks[1] != 0)
//			appendtograph /W= graph1 $spikenumber
//			print "adding to display...",spikenumber
//		endif
//	index = index +1
//while(1)
// print "peaks", peaks
//
////killwindow graph1
////killwindow graph2
setdatafolder root:
end
	
//-------------------------------------------
Function Number_spikes_opto()
//remove300(2)
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) freq, APnum, freq2, APnum2

variable k
for(k=0;k<(numtrace(waves));k = k+1)

	wave twave = waverefindexed(agraph,k, 1)
make /O /D /N= (200)  peaktime
variable i = 0 , xpos1 =1000,,j;
do
	//findpeak /Q /B=0.0001/M=5/R=(xpos1, pnt2x(trace, numpnts(trace))) twave
	findpeak /Q /B=0.0001/M=5/R=(xpos1, 2000) twave
	if (V_flag != 0)  
		break
	endif
	i = i+1;
	peaktime[i] = V_peakLoc
	
	xpos1 = V_peakLoc+1;
while(1)
variable AP=(i)
variable freq1 = Ap/0.8
Freq [k] = freq1
APnum [k] = AP
endfor

//code added for Karuna's paper analysis on 11/15/18

for(k=0;k<(numtrace(waves));k = k+1)

	wave twave = waverefindexed(agraph,k, 1)
make /O /D /N= (200)  peaktime
variable q = 0 , xpos2 =3000,,p;
do
	//findpeak /Q /B=0.0001/M=5/R=(xpos1, pnt2x(trace, numpnts(trace))) twave
	findpeak /Q /B=0.0001/M=5/R=(xpos2, 4000) twave
	if (V_flag != 0)  
		break
	endif
	q = q+1;
	peaktime[q] = V_peakLoc
	
	xpos2 = V_peakLoc+1;
while(1)
variable AP2=(q)
variable freq_2 = Ap/0.8
Freq2 [k] = freq_2
APnum2 [k] = AP2
endfor




make /o/n=20 currentinj
currentinj [0] = 0
variable m 
for(m=1;m<(20);m = m+1)
currentinj [m] = currentinj [m-1] + 25

endfor
string newgraph
display /N=newgraph, freq vs currentinj 
//DoWindow/F newgraph // Bring graph to front
//usercursoradj(newgraph);
//curvefit  /N/Q   line  freq(xcsr (A, newgraph), xcsr (B, newgraph)) /D

//variable index = 0,first =0,maxap = 0;
//string agraph= winname(0,1);
//string atrace, spikenumber
//SVAR c1folder = root:c1folder
//setdatafolder c1folder
//
//make /O /D /N= (1)  peaks
//do
//	atrace = wavename(agraph,index,1)
//	if (waveexists($atrace) == 0)
//		break
//	endif
//	
//	spikenumber = "b" + atrace +  ">" + num2str(index+1) 
//	peaks =  countAP($atrace)
//	duplicate /O peaks $spikenumber
//	
//		if ((peaks[1] != 0) && (first == 0))
//		        first = 1
//			display $spikenumber 
//			print "displaying...",spikenumber
//			
//		elseif (peaks[1] != 0)
//			appendtograph /W= graph1 $spikenumber
//			print "adding to display...",spikenumber
//		endif
//	index = index +1
//while(1)
// print "peaks", peaks
//
////killwindow graph1
////killwindow graph2
setdatafolder root:
end
	


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Function Paired_Pulse_Ratio()
NVAR freq = s_freq
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) PP1
make /O/D/N= (numtrace(waves)) PP2, baseline, PP_1, PP_2, baseline2, time1, time2, timediff

showinfo /W = $agraph
NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
Print " Indicate the window by placing the cursors to calculate average Vm"
//usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)

	wave twave = waverefindexed(agraph,k, 1)
	wavestats/Z/R= (10,30) twave
	baseline[k]= V_avg
	wavestats /Z /R=(50,70) twave
	PP_1[k] = V_max
	time1[k]= V_maxloc
	PP1[k]= PP_1[K]-baseline[k];
endfor
//usercursoradj(agraph);
variable m
for(m=0;m<(numtrace(waves));m = m+1)
	wave twave = waverefindexed(agraph,m, 1)
	wavestats/Z/R= (10,30) twave
	baseline2[m]= V_avg
	wavestats /Z /R=(70,400) twave
	PP_2[m] = V_max
	time2[m]= V_maxloc
	PP2[m]=PP_2[m]-baseline2[m]
endfor
make /O/D/N= (numtrace(waves)) PPR

variable j
for (j=0;j<(numtrace(waves));j = j+1)
PPR[j]=PP2[j]/PP1[j]
timediff[j]= time2[j]- time1[j]


endfor
Variable Paired_Pulse_Ratio, time_diff
Paired_Pulse_Ratio= mean(PPR,0,Pi)
time_diff= mean(timediff,0,Pi)

Print " The average PPR for this cell is", Paired_Pulse_Ratio
Print "the time diff is", time_diff

display PPR

end
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

function averagethetraces() 
make /O /D /N= ( numpnts(waverefindexed(winname(0,1),0,1)))  avgedtrace = 0; // make a new wave with the same points as the first wave in the graph window     
variable index=0
do
wave awave = waverefindexed(winname(0,1),index,1) 
 avgedtrace = awave+avgedtrace
	index = index+1;
while(waveexists(waverefindexed(winname(0,1),index,1)))
avgedtrace = avgedtrace/(index-1)
copyscales /P awave, avgedtrace
//setscale /P x 0,50E-06, " s ", avgedtrace
//setscale /P y 0,1, " A ", avgedtrace
print "the number of waves are ",index
end
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Function EPSP_Latency()


NVAR freq = s_freq
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) SA
make /O/D/N= (numtrace(waves)) Peak

showinfo /W = $agraph
NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
Print " Indicate the window by placing the cursors to calculate average Vm"
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	SA[k] = V_maxloc;
endfor
usercursoradj(agraph);
variable m
for(m=0;m<(numtrace(waves));m = m+1)
	wave twave = waverefindexed(agraph,m, 1)
	wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	Peak[m] = V_maxloc;
endfor
make /O/D/N= (numtrace(waves)) Latency

variable j
for (j=0;j<(numtrace(waves));j = j+1)
Latency[j]=peak[j]-SA[j]


endfor
//Variable Paired_Pulse_Ratio
//Paired_Pulse_Ratio= mean(PPR,0,Pi)

//Print " The average PPR for this cell is", Paired_Pulse_Ratio

display Latency
setscale /P   x, 0, 0.3, Latency
SetAxis left 0.005,0.025 
Label left "\\F'Arial'\\Z14\\K(39168,0,0)\\f01 EPSP Latency"
Label bottom "\\F'Arial'\\Z14\\K(39168,0,0)\\f01Time (min)"
ModifyGraph nticks=5
ModifyGraph mode(latency)=2,lsize(latency)=6
ModifyGraph rgb(latency)=(26112,0,0)
//ModifyGraph rgb=(0,15872,65280),lstyle(newres)=2

end

//-----------------------------------------------------------------------------------------------------------------------------------------------
function  Input_resistance()
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) res_in,Vchange
showinfo /W = $agraph
NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
Print " Place the cursors in oder to measure the resistance jump"
usercursoradj(agraph);
variable k,currinj = 1E-12
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	Vchange[k] = V_avg;
	//print " Maximum of the jump = ",V_max
	//Print "Minimum of the jump = ", V_min
	//print " The voltage change = ",Vchange[k]
endfor
prompt CurrInj,"Enter the Hyperpolarizing current injected"
doprompt  "Current(pA)", CurrInj 
res_in = (-Vchange/CurrInj) 


 
 Variable avgres_in=(mean(res_in,0,14))
print avgres_in
Variable avg_in
avg_in=100/avgres_in
print avg_in
make /O/D/N=(numtrace (waves)) newres_in

for(k=0;k<(numtrace(waves));k = k+1)


	newres_in=res_in*avg_in

	endfor

display newres_in
 setscale /P   x, 0, 0.3, newres_in

Appendtograph /W=ser_res, newres_in
SetAxis left -50,250
Label left "\\Z14\\K(26112,0,10240) % Input resistance"
Label bottom "\\Z14\\K(26112,0,10240)Time (min)"

ModifyGraph nticks=5
ModifyGraph mode(newres_in)=2,lsize(newres_in)=6

ModifyGraph rgb=(26112,0,10240),lstyle(newres_in)=2
SetAxis bottom 0,32 
// setscale /P   x, 0, nuslopes;
//Variable avgres1
//avgres1= mean(res,0,Pi)

//make  /O/D/N= (numtrace(waves)) 
//print avgres
end
//killwaves Vchange
//------------------------------------------------------------------------------------------------------------------------------------------

Function AMPA_peak()
variable choice=1
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) peak, timeofpeak, TOP1,TOP2, peak1
showinfo /W = $agraph
prompt choice, "Enter 1 if baseline need to be adjusted"
doprompt  "Your choice is", choice
if (choice==1)
adjustbaseline()
endif
NVAR slope = root:k1
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Peak[k]=V_min
Timeofpeak [k]=V_minLoc
TOP1 [k]=Timeofpeak [K]-0.0005
TOP2 [k]=Timeofpeak [k]+0.0005

Wavestats /Z/R=(TOP1[K],TOP2 [k]) twave
Peak1[k] = V_avg

endfor
Variable Average_AMPA
Average_AMPA= mean(peak1,0,Pi)
print average_AMPA

display /N=AMPA, peak1


end
//-----------------------------------------------------//

Function AMPA_area()
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) Area_ampa,SA, Xpos, Ypos

showinfo /W = $agraph
SetAxis bottom 0,0.2 
SetAxis left -8e-10,4e-10 
NVAR slope = root:k1
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
wave twave = waverefindexed(agraph,k, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
SA [k]= V_minloc
Xpos [k]= V_minloc + 0.0001
Ypos [k] = Xpos[k] + 0.4251
Area_ampa [k]=area (twave ,Xpos[k], Ypos[k])
endfor
Variable Average_area, final_average
Average_area= mean(Area_ampa,0,Pi)
Final_average =Average_area*(-1000000000000000)
Print "AMPA_area is", final_average, "pA*ms"
Display /N=AMPA_area, Area_ampa
end
  //------------------------------------------------------------------------------//
  
  Function NMDA_peak ()
  
variable choice=1
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) peak, timeofpeak, TOP1,TOP2, peak1
showinfo /W = $agraph
prompt choice, "Enter 1 if baseline need to be adjusted"
doprompt  "Your choice is", choice
if (choice==1)
adjustbaseline()
endif
NVAR slope = root:k1
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Peak[k]=V_max
Timeofpeak [k]=V_maxLoc
TOP1 [k]=Timeofpeak [K]-0.005
TOP2 [k]=Timeofpeak [k]+0.005

Wavestats /Z/R=(TOP1[K],TOP2 [k]) twave
Peak1[k] = V_avg

endfor
Variable Average_NMDA
Average_NMDA= mean(peak1,0,Pi)
print average_NMDA

display /N=NMDA, peak1

end
//------------------------------------------------------------------------------------------

Function NMDA_area()
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) Area_nmda, SA, Xpos, Ypos

showinfo /W = $agraph
SetAxis bottom 1.3,2.2 
SetAxis left -2e-10,4e-10 
NVAR slope = root:k1
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
wave twave = waverefindexed(agraph,k, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
SA [k]= V_minloc
Xpos [k]= V_minloc + 0.0001
Ypos [k] = Xpos[k] + 0.4251
Area_nmda [k]=area (twave ,Xpos[k], Ypos[k])
endfor
Variable Average_area, final_average
Average_area= mean(Area_nmda,0,Pi)
Final_average =Average_area*(1000000000000000)
Print "NMDA_area is", final_average, "pA*ms"
Display /N=NMDA_area, Area_NMDA
end
  
  
  //-------------------------------------------------------------------------------------------------------------------------------------
  
  Function NMDA_decay()  
    
string agraph= winname(0,1),atrace
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) taus1, taus2, Current1, current2
showinfo /W = $agraph
//SetAxis bottom 1.3,2.2 
//SetAxis left -2e-10,4e-10 
//RemoveFromGraph Area_nmda
//RemoveFromGraph SA
//RemoveFromGraph Ypos
//RemoveFromGraph Xpos
usercursoradj(agraph);
NVAR tau = root:k2
variable i,k
for(k=0;k<(numtrace(waves));k = k+1)
	atrace = stringfromlist(k,waves,";")
	wave twave = waverefindexed(agraph,k, 1)
	CurveFit   /N  dblexp  twave (xcsr (A, agraph), xcsr (B, agraph)) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	taus1[k] = 1/ K2;
	Taus2 [k]=1/k4
	Current1 [k]= k1
	Current2 [k]=k3
endfor
wavestats /R= [0, ] /Q/Z taus1
variable decay1
Decay1= V_Avg
wavestats /R= [0, ] /Q/Z taus2
Variable decay2
decay2= V_Avg
wavestats /R= [0, ] /Q/Z current1
variable curr1
curr1= V_Avg
wavestats /R= [0, ] /Q/Z current2
Variable curr2
curr2= V_Avg
printf " The average time constants of the cell is %g", Decay1
printf " The average time constants of the cell is %g", Decay2
print curr1
print curr2
variable weight
weight=(decay1* (curr1/(curr1+curr2)))+ ((curr2/(curr1+curr2))*decay2)
Print weight
display taus1, taus2
//string fits = wavelist("fit_*",";","WIN:"),trace;
//for(i=0;i<(numtrace(waves));i=i+1)
//trace = stringfromlist (i, fits)
//modifygraph rgb(trace) = (0,0,0); doupdate;
//endfor
End.

//---------------------------------------
Function tau()  
    
string agraph= winname(0,1),atrace
string waves = tracenamelist(agraph,";",1); 
//make /O/D/N= (numtrace(waves)) taus

showinfo /W = $agraph
setaxis bottom 58.9,400
//RemoveFromGraph Area_ampa
//RemoveFromGraph SA
//RemoveFromGraph Ypos
//RemoveFromGraph Xpos
//SetAxis bottom 0,0.202588
//SetAxis left -8e-10,5e-11 
 
usercursoradj(agraph);
//NVAR tau = root:
variable k
for (k=0;k<(numtrace(waves));k = k+1)
	atrace = stringfromlist(k,waves,";")
	wave twave = waverefindexed(agraph,k, 1)
	//wave dvave= Xwavename
	variable x= pcsr (A, agraph)
	variable y = pcsr (B, agraph)
	
	print x, y
	CurveFit/NTHR=0 dblexp_XOffset  twave[X,Y] /D 
	//CurveFit dblexp_XOffset, twave[x, y] /D 
	//CurveFit/M=2/W=0 dblexp_XOffset, graph_avg[1184,3085]/D
	//curvefit  /NTHR=0  line  twave/X= wave13/D
	//modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	//taus[k] = tau;
endfor
End
//wavestats /R= [0, ] /Q/Z taus
//printf " The average time constant of the cell is %g", V_avg 
//display taus
//string fits = wavelist("fit_*",";","WIN:"),trace;
//for(i=0;i<(numtrace(waves));i=i+1)
//trace = stringfromlist (i, fits)
//modifygraph rgb(trace) = (0,0,0); doupdate;
//endfor

//-----------------------------------------------------------------------------------
 
  Function NMDA_MK801()
  
variable choice=1
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) peak, timeofpeak, TOP1,TOP2, peak1
showinfo /W = $agraph
prompt choice, "Enter 1 if baseline need to be adjusted"
doprompt  "Your choice is", choice
if (choice==1)
adjustbaseline()
endif
NVAR slope = root:k1
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Peak[k]=V_max
Timeofpeak [k]=V_maxLoc
TOP1 [k]=Timeofpeak [K]-0.0025
TOP2 [k]=Timeofpeak [k]+0.0025

Wavestats /Z/R=(TOP1[K],TOP2 [k]) twave
Peak1[k] = V_avg

endfor
display /N=NMDA, peak1
ModifyGraph mode=3,marker=8

Variable avg= mean(peak1,0,0)
print avg    
variable newavg=1/avg


make /O/D/N= (numtrace(waves)) nupeak
for(k=0;k<(numtrace(waves));k = k+1)


	nupeak=peak1*newavg 

	endfor
	display /N=NMDA_decay, nupeak
	ModifyGraph mode=4,marker=8
	




end
//------------------------------------------------------------------------------------------
//Function FT()
//FFT /RP=[startPoint,endPoint ] /RX=(startX,endX ) /PAD={dim1 [, dim2, dim3, dim4 ]} /Z  [/COLS/HCC/HCS/ROWS/REAL/MAG/...



//------------------------------------------------------------------------
Function Var()

variable choice=1
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) peak, timeofpeak, TOP1,TOP2, peak1, TOP3, TOP4, Baseline, PeakV1, PeakV2, baseV1, basev2
showinfo /W = $agraph
prompt choice, "Enter 1 if baseline need to be adjusted"
doprompt  "Your choice is", choice
if (choice==1)
adjustbaseline()
endif
NVAR slope = root:k1
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Peak[k]=V_min
Timeofpeak [k]=V_minLoc
TOP1 [k]=Timeofpeak [K]-0.0005
TOP2 [k]=Timeofpeak [k]+0.0005
TOP3 [k]= timeofpeak[k] -0.016
TOP4 [k]= timeofpeak [k]-0.017

Wavestats /Z/R=(TOP1[K],TOP2 [k]) twave
Peak1[k] = V_avg
Wavestats /Z/R=(TOP4[K],TOP3 [k]) twave
baseline[k] = V_avg

endfor
Variable Average_AMPA1, average_ampa
Average_AMPA1= mean(peak1)
Average_ampa= average_AMPA1*1000000000000

Variable Average_base1, average_base
Average_base1= mean(baseline)
Average_base= average_base1*1000000000000


 for(k=0;k<(numtrace(waves));k = k+1)
PeakV1 [k] = peak1 [k]-average_ampa1
peakV2[k]= peakV1[k]*peakV1[k]
endfor
variable peak_var1, peak_variance
peak_var1= mean(peakV2)
peak_variance= (peak_var1*1000000000000000000000000)/(numtrace(waves)-1)

 for(k=0;k<(numtrace(waves));k = k+1)
baseV1 [k] = baseline[k]-average_base1
baseV2[k]= baseV1[k]*baseV1[k]
endfor
variable base_var1, baseline_variance
base_var1= mean(baseV2)
baseline_variance= (base_var1*1000000000000000000000000)/(numtrace(waves)-1)



print "Peak amplitude is ", average_AMPA, "pA"
print " Peak variance is", peak_variance, "pA2"

print "Baseline amplitude is ", average_base, "pA"
print " baseline variance is", baseline_variance, "pA2"




display /N=AMPA, peak1


end
///--------------------------------------------
Function parab()  
    
string agraph= winname(0,1),atrace
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) taus
showinfo /W = $agraph
//RemoveFromGraph Area_ampa
//RemoveFromGraph SA
//RemoveFromGraph Ypos
//RemoveFromGraph Xpos
//SetAxis bottom 0,0.202588
//SetAxis left -8e-10,5e-11 
 

NVAR tau = root:k2
variable i,k
for(k=0;k<(numtrace(waves));k = k+1)
	atrace = stringfromlist(k,waves,";")
	wave twave = waverefindexed(agraph,k, 1)
	curvefit  /X=1 poly 3  , baseline/D
	
endfor

//string fits = wavelist("fit_*",";","WIN:"),trace;
//for(i=0;i<(numtrace(waves));i=i+1)
//trace = stringfromlist (i, fits)
//modifygraph rgb(trace) = (0,0,0); doupdate;
//endfor
End.
//----------------------------------------------------------------------------
//Function delete_APs
string agraph= winname(0,1),atrace
showinfo /W = $agraph
usercursoradj(agraph);
	DeletePoints xcsr (A, agraph)
//===================================
Function TBS()
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
variable peak, P11, P12, peak1, timeofpeak
showinfo /W = $agraph
NVAR slope = root:k1
usercursoradj(agraph);
	wave twave = waverefindexed(agraph,0, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Peak=V_max
Timeofpeak =V_maxLoc
P11 =Timeofpeak -0.005
P12 =Timeofpeak +0.005
Wavestats /Z/R=(P11,P12) twave
Peak1 = V_avg

print peak1

end
//=================================================================================
Function TBS_5()
variable k
for(k=0;k<6;k = k+1)
prompt K, "Enter the value of K"
doprompt  "K is", K
if (K==1)
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
variable peak, P11, P12, peak1, timeofpeak
showinfo /W = $agraph
NVAR slope = root:k1
usercursoradj(agraph);
	wave twave = waverefindexed(agraph,0, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Peak=V_max
Timeofpeak =V_maxLoc
P11 =Timeofpeak -0.00005
P12 =Timeofpeak +0.00005
Wavestats /Z/R=(P11,P12) twave
Peak1 = V_max

print peak1

variable EPSP_1=peak1
endif
if (K==2)


showinfo /W = $agraph
NVAR slope = root:k1
usercursoradj(agraph);
	wave twave = waverefindexed(agraph,0, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Peak=V_max
Timeofpeak =V_maxLoc
P11 =Timeofpeak -0.005
P12 =Timeofpeak +0.005
Wavestats /Z/R=(P11,P12) twave
Peak1 = V_max

print peak1
variable EPSP_2=peak1
endif
if (K==3)

showinfo /W = $agraph
NVAR slope = root:k1
usercursoradj(agraph);
	wave twave = waverefindexed(agraph,0, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Peak=V_max
Timeofpeak =V_maxLoc
P11 =Timeofpeak -0.005
P12 =Timeofpeak +0.005
Wavestats /Z/R=(P11,P12) twave
Peak1 = V_max

print peak1
variable EPSP_3=peak1
endif
if (K==4)

showinfo /W = $agraph
NVAR slope = root:k1
usercursoradj(agraph);
	wave twave = waverefindexed(agraph,0, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Peak=V_max
Timeofpeak =V_maxLoc
P11 =Timeofpeak -0.005
P12 =Timeofpeak +0.005
Wavestats /Z/R=(P11,P12) twave
Peak1 = V_max

print peak1
variable EPSP_4=peak1
endif
if (K==5)

showinfo /W = $agraph
NVAR slope = root:k1
usercursoradj(agraph);
	wave twave = waverefindexed(agraph,0, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Peak=V_max
Timeofpeak =V_maxLoc
P11 =Timeofpeak -0.005
P12 =Timeofpeak +0.005
Wavestats /Z/R=(P11,P12) twave
Peak1 = V_max

print peak1
variable EPSP_5=peak1
endif
endfor

variable average_TBS= (EPSP_1+EPSP_2+EPSP_3+EPSP_4+EPSP_5)/5
Print EPSP_1, EPSP_2, EPSP_3, EPSP_4, EPSP_5, Average_TBS
end
//==================================================================================================
Function EPSPamp()
//for Current clamp LTP expts only, select 15 traces for pre n 90 for post
// measure the slope in 10-90 % of the rising phase of the EPSP
string agraph= winname(0,1),slopename,atrace

string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) peak, timeofpeak, TOP1, TOP2, peak1
showinfo /W = $agraph
variable choice=1
prompt choice, "Enter 1 if baseline need to be adjusted"
doprompt  "Your choice is", choice
if (choice==1)
adjustbaseline()
endif

NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
usercursoradj(agraph);
variable k,x80percentpeak,Vchange,x20percentpeak
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	atrace = stringfromlist(k,waves,";")
	wavestats /Q/Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	Peak[k]=V_max
Timeofpeak [k]=V_maxLoc
TOP1 [k]=Timeofpeak [K]-0.0005
TOP2 [k]=Timeofpeak [k]+0.0005

Wavestats /Z/R=(TOP1[K],TOP2 [k]) twave
Peak1[k] = V_avg
endfor
//insert a single point showing the 100 hz LTP data point
variable base
prompt base,"Enter the number of traces in baseline"
doprompt "enter the name madi", base
InsertPoints base,3, peak1
string tempgraph = "temp1"

display /n= $tempgraph peak1

string peakname
prompt peakname," Enter the name of the slope wave now"
doprompt  "Simply enter the name madi", peakname
make /O/D/N= (numtrace(waves)) $peakname
duplicate /O peak1 $peakname
ModifyGraph mode=2,lsize=6

killwindow  $tempgraph

variable trace
trace= base-1
Variable avg_baseline= mean($peakname,0,trace)
print avg_baseline    




end.
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

Function line()  
    
string agraph= winname(0,1),atrace
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) taus
showinfo /W = $agraph
//RemoveFromGraph Area_ampa
//RemoveFromGraph SA
//RemoveFromGraph Ypos
//RemoveFromGraph Xpos
//SetAxis bottom 0,0.202588
//SetAxis left -8e-10,5e-11 
 
usercursoradj(agraph);
variable x = xcsr (A, agraph)
variable y = xcsr (B, agraph)

NVAR tau = root:k2
variable i,k
for(k=0;k<(numtrace(waves));k = k+1)
	atrace = stringfromlist(k,waves,";")
	wave twave = waverefindexed(agraph,k, 1)
	//wave dwave= XWavename (agraph, atrace)
	CurveFit/M=2/W=0 line, freq[x,y]/X=currentinj[x,y]/D
	//curvefit  /NTHR=0   line  twave/X=wave13 /D
	//modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	taus[k]=tau

endfor

End.
//==========================================================================
Function EPSPampfull()
//for Current clamp LTP expts only, select 15 traces for pre n 90 for post
// measure the slope in 10-90 % of the rising phase of the EPSP
string agraph= winname(0,1),slopename,atrace

string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) peak1, peak, timeofpeak, TOP1, TOP2
showinfo /W = $agraph
variable choice=1
prompt choice, "Enter 1 if baseline need to be adjusted"
doprompt  "Your choice is", choice
if (choice==1)
adjustbaseline()
endif

NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	atrace = stringfromlist(k,waves,";")
	wavestats /Q/Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	Peak[k]=V_max
Timeofpeak [k]=V_maxLoc
TOP1 [k]=Timeofpeak [K]-0.0005
TOP2 [k]=Timeofpeak [k]+0.0005

Wavestats /Z/R=(TOP1[K],TOP2 [k]) twave
Peak1[k] = V_avg
endfor
//insert a single point showing the 100 hz LTP data point
variable base
prompt base,"Enter the number of traces in baseline"
doprompt "enter the name madi", base
InsertPoints base,3, peak1
string tempgraph = "temp1"

display /n= $tempgraph peak1

string peakname
prompt peakname," Enter the name of the slope wave now"
doprompt  "Simply enter the name madi", peakname
make /O/D/N= (numtrace(waves)) $peakname
duplicate /O peak1 $peakname
ModifyGraph mode=2,lsize=6

killwindow  $tempgraph

variable trace
trace= base-1
Variable avg_baseline= mean($peakname,0,trace)
print avg_baseline    
	

variable newavg=100/avg_baseline // calculate % average of baseline

//make new wave by dividing all original slope values by the %average of the baseline
make /O/D/N= (numtrace(waves)) nuamps
for(k=0;k<(numtrace(waves));k = k+1)


nuamps=peak1*newavg 

	endfor
 //set axis, labels etc
 display /n=$peakname nuamps
 //setscale /P   x, 0, 0.3, nuslopes; 
 ModifyGraph mode=2,lsize=6
 //SetAxis left 0.1,300 
 //SetAxis bottom 0,32 

 Label left "\\Z14 % normalised EPSP amplitude (mV)"
Label bottom "\\Z14 Time (min)"
end.
//=======================================================================================================================//
Function E_S_pot()
// measure the slope in 10-90 % of the rising phase of the EPSP


string agraph= winname(0,1),slopename,atrace
//SVAR c1folder = root:c1folder
//setdatafolder c1folder
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) slopes , baseline, peak
showinfo /W = $agraph
//SetAxis left -0.005, 0.05
//SetAxis bottom 0.01, 0.1
//variable choice=1
//prompt choice, "Enter 1 if baseline need to be adjusted"
//doprompt  "Your choice is", choice
//if (choice==1)
//adjustbaseline()
//endif

//usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	atrace = stringfromlist(k,waves,";")
	wavestats /Q/Z /R=(0,15) twave
	baseline[k]=V_avg
endfor
//SetAxis left -0.0015, 0.02
//SetAxis bottom 0.015, 0.035
NVAR slope = root:k1
//usercursoradj(agraph);

variable x80percentpeak,Vchange,x20percentpeak, xpos, ypos
variable stim=1
prompt stim, "Enter 1 if stim is at 100"
doprompt  "Your choice is", stim
if (stim==1)

for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	atrace = stringfromlist(k,waves,";")
	//wavestats /Q/Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	//Vchange = V_max-V_min;
	//xpos = xcsr(A, agraph)
	//ypos = xpos+ 2
	CurveFit/M=2/W=0 line, twave(103,105) /D
	modifygraph rgb ($("fit_"+ atrace))= (0,0,0);
	slopes[k]= slope
	endfor
else
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	atrace = stringfromlist(k,waves,";")
	//wavestats /Q/Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	//Vchange = V_max-V_min;
	//xpos = xcsr(A, agraph)
	//ypos = xpos+ 2
CurveFit/M=2/W=0 line, twave(52.5,53.5) /D
modifygraph rgb ($("fit_"+ atrace))= (0,0,0);
slopes[k]= slope
	endfor
endif
	//CurveFit/M=2/W=0 line, twave(53,55) /D
	//curvefit  /N/Q   line  twave(xpos,ypos) /D
	//modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	
//SetAxis left -0.005, 0.05
//SetAxis bottom 0.01, 0.1

//usercursoradj(agraph);

for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	atrace = stringfromlist(k,waves,";")
	wavestats /Q/Z /R=(20,200) twave
	Peak[k]=V_max
endfor

variable v
for (v=100; v>100;v=200)
findvalue /U=(100-150) peak
endfor
//setdatafolder root:

end.
//========================================================================================================================//
Function EPSCampfull()
//for voltage clamp LTP expts only, select 15 traces for pre n 90 for post
// measure the slope in 10-90 % of the rising phase of the EPSP
string agraph= winname(0,1),slopename,atrace

string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) peak1, ppeak1, timeofpeak1, TOP11, TOP12, peak2, ppeak2, timeofpeak2, TOP21,TOP22, PPR
showinfo /W = $agraph
variable choice=1
prompt choice, "Enter 1 if baseline need to be adjusted"
doprompt  "Your choice is", choice
if (choice==1)
adjustbaseline()
endif

NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	atrace = stringfromlist(k,waves,";")
	wavestats /Q/Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	Peak1[k]=V_min
Timeofpeak1[k]=V_minLoc
TOP11 [k]=Timeofpeak1 [K]-0.0005
TOP12 [k]=Timeofpeak1 [k]+0.0005

Wavestats /Z/R=(TOP11[K],TOP12 [k]) twave
PPeak1[k] = V_avg
endfor
usercursoradj(agraph);
variable m
for(m=0;m<(numtrace(waves));m = m+1)
	wave twave = waverefindexed(agraph,m, 1)
	wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
		Peak2[m]=V_min
Timeofpeak2[m]=V_minLoc
TOP21 [m]=Timeofpeak2 [m]-0.0005
TOP22 [m]=Timeofpeak2 [m]+0.0005

Wavestats /Z/R=(TOP21[m],TOP22 [m]) twave
PPeak2[m] = V_avg
endfor
make /O/D/N= (numtrace(waves)) PPR

variable j
for (j=0;j<(numtrace(waves));j = j+1)
PPR[j]=PPeak2[j]/PPeak1[j]
endfor
//insert a single point showing the 100 hz LTP data point
variable base
prompt base,"Enter the number of traces in baseline"
doprompt "enter the name madi", base
InsertPoints base,3, ppeak1
insertpoints base,3, PPR
string tempgraph = "temp1"

display /n= $tempgraph peak1

string peakname
prompt peakname," Enter the name of the slope wave now"
doprompt  "Simply enter the name madi", peakname
make /O/D/N= (numtrace(waves)) $peakname
duplicate /O peak1 $peakname
ModifyGraph mode=2,lsize=6

killwindow  $tempgraph

variable trace
trace= base-1
Variable avg_baseline= mean($peakname,0,trace)
print avg_baseline    
	

variable newavg=100/avg_baseline // calculate % average of baseline

//make new wave by dividing all original slope values by the %average of the baseline
make /O/D/N= (numtrace(waves)) nuamps
for(k=0;k<(numtrace(waves));k = k+1)


nuamps=ppeak1*newavg 

	endfor
 //set axis, labels etc
 display /n=$peakname nuamps
 //setscale /P   x, 0, 0.3, nuslopes; 
 ModifyGraph mode=2,lsize=6
 //SetAxis left 0.1,300 
 //SetAxis bottom 0,32 

 Label left "\\Z14 % normalised EPSP amplitude (mV)"
Label bottom "\\Z14 Time (min)"



string tempgraph2 = "temp"

display /n= $tempgraph peak1

string peakname2
prompt peakname2," Enter the name of the slope wave now"
doprompt  "Simply enter the name madi", peakname2
make /O/D/N= (numtrace(waves)) $peakname2
duplicate /O PPR $peakname2
ModifyGraph mode=2,lsize=6

killwindow  $tempgraph


Variable avg_baseline2= mean($peakname2,0,trace)
print avg_baseline2    
	

variable newavg2=100/avg_baseline2 // calculate % average of baseline

//make new wave by dividing all original slope values by the %average of the baseline
make /O/D/N= (numtrace(waves)) nuPPR
for(k=0;k<(numtrace(waves));k = k+1)


nuPPR=PPR*newavg2

	endfor
 //set axis, labels etc
 display /n=$peakname2 nuPPR
 //setscale /P   x, 0, 0.3, nuslopes; 
 ModifyGraph mode=2,lsize=6
 //SetAxis left 0.1,300 
 //SetAxis bottom 0,32 

 Label left "\\Z14 % normalised PPR"
Label bottom "\\Z14 Time (min)"
end.
//---------------------------------------------------------------------------------------------------------------------------------//
Function Ih_plots()
NVAR freq = s_freq
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) SSV
make /O/D/N= (numtrace(waves)) sag
make /O/D/N= (numtrace(waves)) rebound
showinfo /W = $agraph
NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
Print " Indicate the window by placing the cursors to calculate average Vm"
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	SSV[k] = V_avg;
endfor
usercursoradj(agraph);
variable m
for(m=0;m<(numtrace(waves));m = m+1)
	wave twave = waverefindexed(agraph,m, 1)
	wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	Sag[m] = V_min;
endfor
usercursoradj(agraph);
variable l
for(l=0;l<(numtrace(waves));l = l+1)
	wave twave = waverefindexed(agraph,l, 1)
	wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	rebound[l] = V_max;
endfor
make /O/D/N= (numtrace(waves)) sagR
make /O/D/N= (numtrace(waves)) sagV

variable j
for (j=0;j<(numtrace(waves));j = j+1)
SagR[j]=SSV[j]/sag[j]
SagV[j]=SSV[j]-sag[j]

endfor
Variable Sag_ratio
Sag_ratio= mean(sagR,0,Pi)
Variable Sag_voltage
Sag_voltage=mean(sagV,0,Pi)
Print " The average sag ratio for this cell is", Sag_ratio
Print " The average sag voltage for this cell is", Sag_voltage,"mV"
display SagV vs SSV
display Rebound vs SSV


end

//___________________________------------------------------------------------------------//
Function Zero_Sub()

string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
showinfo /W = $agraph
NVAR slope = root:k1
variable average
average= mean (twave, 20, 120) // changes pi to 700
twave= twave-average
endfor
end

//--------------------------------------------------------------------------------------------------------------------------------------------------------------//
Function Res_freq()
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) fourier
variable k
wave twave = waverefindexed(agraph,K,1)

FFT   twave


end
//-----------------------------------------------------------------------------------------------//

Function resonance()
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) peak
showinfo /W = $agraph
NVAR slope = root:k1
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave
Peak[k] = V_maxloc
endfor
Variable Average_AMPA
Average_AMPA= (mean(peak,0,Pi))*0.05
print average_AMPA

display /N=res, peak
end
//=====================================================================================================================//
Function AHP_area()
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) Area_ahp, Xpos, Ypos
 SetAxis left -0.004,0.004
 SetAxis bottom 0.45,1
showinfo /W = $agraph
ShowTools/A

NVAR slope = root:k1
usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
wave twave = waverefindexed(agraph,k, 1)
Wavestats /Z/R=(xcsr (A, agraph), xcsr (B, agraph)) twave

Xpos [k]= (xcsr (A, agraph))
Ypos [k] =(xcsr (B, agraph))
Area_ahp [k]=area (twave ,Xpos[k], Ypos[k])
endfor
Variable Average_area, final_average
Average_area= mean(Area_ahp,0,Pi)
Final_average =Average_area
Print "AHP_area is", final_average, "pA*ms"
Display /N=AMPA_area, Area_ampa
end
//=============================================================================================================//
Function measureCC()
//measures the series resistance jump specified by cursors for both baseline and the jump region,
// Then asks the user for the current -inj 
// Asks if any averages has to be done and then plots a graph
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) res,Vchange
showinfo /W = $agraph
variable choice=1
prompt choice, "Enter 1 if baseline need to be adjusted"
doprompt  "Your choice is", choice
if (choice==1)
adjustbaseline()
endif
NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
Print " Place the cursors in oder to measure the resistance jump"
usercursoradj(agraph);
variable k,currinj = 1E-12
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	Vchange[k] = V_avg;
	print " Maximum of the jump = ",V_max
	Print "Minimum of the jump = ", V_min
	print " The voltage change = ",Vchange[k]
endfor
prompt CurrInj,"Enter the Hyperpolarizing current injected"
doprompt  "Current(pA)", CurrInj 
res = (Vchange/CurrInj) *(-1)
display res


 Variable avgres=(mean(res,0,Pi))
print avgres
Variable avg
avg=100/avgres
print avg
make /O/D/N=(numtrace (waves)) newres

for(k=0;k<(numtrace(waves));k = k+1)


	newres=res*avg

	endfor

display res
 setscale /P   x, 0, 0.3, res

Appendtograph /W=ser_res, res
SetAxis left 80e6,120e6
Label left "\\Z14\\K(0,15872,65280) Series Resistance"
Label bottom "\\Z14\\K(0,15872,65280)Time (min)"
ModifyGraph nticks=5
ModifyGraph mode(res)=2,lsize(res)=6
ModifyGraph rgb(res)=(0,15872,65280)
ModifyGraph rgb=(0,15872,65280),lstyle(res)=2


SetAxis bottom 0,32 
// setscale /P   x, 0, nuslopes;
//Variable avgres1
//avgres1= mean(res,0,Pi)

//make  /O/D/N= (numtrace(waves)) 
//print avgres
end
//killwaves Vchange
//---------------------------------------------------------------------------------------------------
Function measVm()
//measures the series resistance jump specified by cursors for both baseline and the jump region,
// Then asks the user for the current -inj 
// Asks if any averages has to be done and then plots a graph
string agraph= winname(0,1)
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) res,Vchange
showinfo /W = $agraph
variable choice=1
prompt choice, "Enter 1 if baseline need to be adjusted"
doprompt  "Your choice is", choice
if (choice==1)
adjustbaseline()
endif
NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
Print " Place the cursors in oder to measure the resistance jump"
usercursoradj(agraph);
variable k,currinj = 1E-12
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	Vchange[k] = V_avg;
	print " Maximum of the jump = ",V_max
	Print "Minimum of the jump = ", V_min
	print " The voltage change = ",Vchange[k]
endfor
prompt CurrInj,"Enter the Hyperpolarizing current injected"
doprompt  "Current(pA)", CurrInj 
res = Vchange
display res


 
end
//killwaves Vchange

//------------------------------------------------------------------------------------------------------------------------------------------
Function Threshold()
string agraph= winname(0,1)
String trace= wavename(agraph,0,1)

Differentiate $trace/D=W_DIF
Display W_DIF vs $trace
//ModifyGraph axisOnTop(left)=1
//ModifyGraph grid(left)=1
ModifyGraph grid=1,nticks=10
ModifyGraph nticks(bottom)=20
setdatafolder root:
end


//--------------------------------------------------------------------------------------------------------------------------
//Function ZAP()
//SVAR c1folder = root:c1folder
//SVAR c2folder=root:c2folder
//variable folder
//prompt folder, "channel 1 or2"
//doprompt  "Your choice is", folder
//if (folder==1)
//setdatafolder c1folder
//else 
//setdatafolder c2folder
//endif
//string agraph= winname(0,1);
//string waves = tracenamelist(agraph,";",1)
//make/o/D/N=227 avg_phase
//make /O/D/N= (numtrace(waves)) zapfreq
//	wave atrace = waverefindexed(agraph,0,1)
//	variable m
//	for(m=0;m<(numtrace(waves));m = m+1)
//	wave twave = waverefindexed(agraph,m, 1)
//	
//	zapfreq[m]=Zapfrequency(twave)
//	if (m==0)
//	duplicate /o ZAPPhase, zphase1	
//	elseif (m==1)
//	duplicate /o ZAPPhase, zphase2	
//	elseif (m==2)
//	duplicate /o ZAPPhase, zphase3	
//	endif
//	endfor
//	
//	
//	variable trace
//	prompt trace, "number of traces"
//doprompt  "Your choice is", trace
//if (trace==3)
//	avg_phase= (zphase1+zphase2+zphase3)/3
//	elseif (trace==2)
//	avg_phase= (zphase1+zphase2)/2
//	elseif (trace==1)
//	avg_phase= (zphase1)
//	endif
//	SetScale/I x 0, 15,"",  avg_phase
//	Duplicate /O avg_phase, PhaseSmooth
//	smooth 4, PhaseSmooth
//	
//	setdatafolder root:
//	end
//	
//--------------------------------------------------------------------------------------------------------------------------	
	

 Function AP_props()


 string agraph= winname(0,1)
 String trace= wavename(agraph,0,1)
wave atrace = waverefindexed(agraph,0,1)
SVAR c1folder = root:c1folder
setdatafolder c1folder
//Zero_Sub()
showinfo /W = $agraph;

//automating to find the first AP
make /O /D /N= (2)  peaktime
variable  xpos1 =0.01;

	findpeak /Q /B=0.0001/M=5/R=(xpos1, pnt2x(atrace, numpnts(atrace))) atrace
	//findpeak /Q /B=1 /M=0.02 /R=(xpos1, 0.7) trace
	
	//setting the timewindow to measure the AP props
peaktime = V_peakLoc 
	variable x = peaktime[0]-2
	variable y= peaktime[0] +3
	setaxis bottom x-10, y+10
	
	

//usercursoradj(agraph);
//variable x= xcsr(A,agraph)
//variable y= xcsr(B,agraph)
variable z= x+0.5
 //Display W_DIF vs $atrace
 Differentiate $trace/D=W_DIF
  Differentiate W_DIF/D=W_DIF2
  Differentiate W_DIF2/D=W_DIF3
 //FindLevel /Q/R=(x,y) W_DIF, 20						/// finds level where first derivative = 20 mV/ms
	wavestats /Q/R=(x,z) W_DIF3
	variable q= V_max
	FindLevel /Q/R=(x,y) W_DIF3, q
					print "V_levelX", V_LevelX
					//variable xsample_interval= 100000* (root:,fADCSampleInterval)
					NVAR xsample_interval=root:fADCSampleInterval
					//variable a= 100000
					//NVAR xsample_interval =xsample_interval1
	Variable AP_thresh_wave=atrace[V_LevelX*1000/xsample_interval]	
//WaveStats/Q/R=(100,300) W_DIF
WaveStats/Q/R=(x,y) W_DIF							/// finds maximum rate of rise of the first spike
	Variable peak_dVdt=V_max	
make /O/D/N=2 spikehalfwidth
WaveStats/Q/R=(0,30) atrace							/// finds maximum rate of rise of the first spike
	Variable resting= V_Avg
//variable peak_width, peak_amp,Half_amp,V_rest,peak_amp_rest
showinfo /W = $agraph;
//usercursoradj(agraph);
//wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) atrace
//variable V_rest=V_avg
//usercursoradj(agraph);
	//findpeak /Q /B=0.5 /M=0.005 /R=(xcsr (A, agraph), xcsr (B, agraph)) trace
	//findpeak /Q /B=0.5 /M=0.001 /R=(190,210) trace
	
	//wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) atrace
	wavestats /Z /R=(x,y) atrace
	variable top = V_maxloc
	variable latency2spike= top-118
	variable peak_amp= V_max
	variable peak_amp_rest=peak_amp-(AP_thresh_wave)
	//variable peak_amp_rest=peak_amp-(AP_thresh_wave)
	variable half_amp= (0.5*peak_amp_rest)+ AP_thresh_wave
	//variable half_amp= (0.5*peak_amp_rest)+ resting
	//usercursoradj(agraph);
//Zero_Sub()
	FindLevels/D=spikehalfwidth/M=0.05/Q/R=(x,y)/T=0.01atrace half_amp 
	
variable half_amp_dur = spikehalfwidth[1] - spikehalfwidth[0]	
//usercursoradj(agraph);
wavestats /Z /R=(x,y) atrace
variable fAHP1= V_min
variable fAHP=fAHP1-AP_thresh_wave
variable fAHP_half= fAHP/2 + AP_thresh_wave

FindLevel /Q/R=(x,y) atrace, fAHP1
variable a = V_levelX
print a
variable b = a + 10
variable w= Ap_thresh_wave + 0.1*fAHP
FindLevel /Q/R=(a,b) atrace, w 
variable c = V_levelX
print c, x
variable half_fAHP = c - a	
make /O/D/N=10 result 
result[0] = peak_amp_rest
 result[1] = peak_dVdt
 result [2] = Ap_thresh_wave
 result [3] = half_amp_dur
 result [4] = fAHP
 result [5] = half_fAHP
 result [6] = latency2spike
	Print "peak_amp_rest, peak_dVdt, AP_thresh_wave, half_amp_dur, fAHP, half_fAHP"
	Print peak_amp_rest, peak_dVdt, AP_thresh_wave, half_amp_dur, fAHP, half_fAHP
	
	//Edit/K=0 root:result
	
	setdatafolder root:
	print peaktime[0], x, y
	
	end
	
	//----------------------------------------------------------------------
	

Function EPSPprops()
Setaxis bottom 48,60;
SVAR c1folder = root:c1folder
SVAR c2folder=root:c2folder
variable folder
prompt folder, "channel 1 or2"
doprompt  "Your choice is", folder
if (folder==1)
setdatafolder c1folder
else 
setdatafolder c2folder
endif
string agraph= winname(0,1),slopename,atrace
Setaxis bottom 48,60;
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves))  RMP, SSV, Inputres, peak, TOP1, TOP2, Timeofpeak,Peak1, amplitude2, peak2, timeofpeak2, TOP3, TOP4, slopes2, amplitude1, slopes
showinfo /W = $agraph
variable m
for(m=0;m<(numtrace(waves));m = m+1)
	wave twave = waverefindexed(agraph,m, 1)
	atrace = stringfromlist(m,waves,";")
wavestats /Z /R=(5,40) twave
	RMP[m] = V_avg;
	wavestats /Z/R= (1690,1699)twave
	SSV[m]= V_avg
	Inputres[m]= (RMP[m]-SSV[m])*(10)
	
	
	wavestats /Q/Z /R=(51,90) twave
	Peak1[m]=V_max
Timeofpeak[m]=V_maxLoc
TOP1 [m]=Timeofpeak [m]-0.0005
TOP2 [m]=Timeofpeak [m]+0.0005

Wavestats /Z/R=(TOP1[m],TOP2 [m]) twave
Peak1[m] = V_avg
Amplitude1[m]= Peak1[m]-RMP[m]

wavestats /Q/Z /R=(799,850) twave
	Peak2[m]=V_max
Timeofpeak2[m]=V_maxLoc
TOP3 [m]=Timeofpeak2 [m]-0.0005
TOP4 [m]=Timeofpeak2 [m]+0.0005

Wavestats /Z/R=(TOP3[m],TOP4 [m]) twave
Peak2[m] = V_avg
Amplitude2[m]= Peak2[m]-RMP[m]

endfor

Setaxis bottom 45,100;
NVAR slope = root:k1

variable x
prompt x, "slope1 x axis 1"
doprompt  "Your choice is", x
variable y=x+1
//prompt y, "slope x axis 2"
//doprompt  "Your choice is", y
variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	atrace = stringfromlist(k,waves,";")
	
	curvefit  /N/Q   line  twave(x,y) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	//removefromgraph $ ("fit_" + atrace)
	slopes[k] = slope;
endfor
//killwaves $ ("fit_" + atrace)
Setaxis bottom 795,820;
NVAR slope = root:k1


variable x1
prompt x1, "slope2 x axis 1"
doprompt  "Your choice is", x1
variable y1=x1+1
//prompt y1, "slope x axis 2"
//doprompt  "Your choice is", y1
for(k=0;k<(numtrace(waves));k = k+1)

	wave twave = waverefindexed(agraph,k, 1)
	atrace = stringfromlist(k,waves,";")
	
	curvefit  /N/Q   line  twave(x1,y1) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	slopes2[k] = slope;
endfor

killwaves SSV, peak,TOP1,TOP2, Timeofpeak, Peak1, Peak2, TOP3, TOP4

setdatafolder root:

end

//---------------------------------------------------------------
Function TBP_area()
SVAR c1folder = root:c1folder
setdatafolder c1folder

string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) Area_TBP, Xpos, Ypos
 //SetAxis left -0.004,0.004
 //etAxis bottom 0.45,1
showinfo /W = $agraph
//ShowTools/A

NVAR slope = root:k1
//usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)
wave twave = waverefindexed(agraph,k, 1)
Wavestats /Z/R=(0,20000) twave

//Xpos [k]= (xcsr (A, agraph))
//Ypos [k] =(xcsr (B, agraph))
Area_tbp [k]=area (twave ,0,20000)
endfor
Variable TBP_area
TBP_area= mean (Area_tbp,0,Pi)
Print "TBP_area for this cell is ", TBP_area

setdatafolder root:
end



//---------------------------
Function backAP_props()


 string agraph= winname(0,1)
 String trace= wavename(agraph,0,1)
wave atrace = waverefindexed(agraph,0,1)
SVAR c1folder = root:c1folder
SVAR c2folder=root:c2folder
variable folder
prompt folder, "channel 1 or2"
doprompt  "Your choice is", folder
if (folder==1)
setdatafolder c1folder
else 
setdatafolder c2folder
endif

//Zero_Sub()
showinfo /W = $agraph;

usercursoradj(agraph);
variable x= xcsr(A,agraph)
variable y= xcsr(b,agraph)
 //Display W_DIF vs $atrace
 Differentiate $trace/D=W_DIF
 FindLevel /Q/R=(x,y) W_DIF, 0.5						/// finds level where first derivative = 20 mV/ms
					print V_LevelX
					NVAR xsample_interval=root:xsample_interval
	Variable AP_thresh_wave=atrace[V_LevelX*1/xsample_interval]	
WaveStats/Q/R=(xcsr (A, agraph), xcsr (B, agraph)) W_DIF
//WaveStats/Q/R=(500,530) W_DIF							/// finds maximum rate of rise of the first spike
	Variable peak_dVdt=V_max	
make /O/D/N=2 spikehalfwidth
WaveStats/Q/R=(0,30) atrace							/// finds maximum rate of rise of the first spike
	Variable resting= V_Avg
//variable peak_width, peak_amp,Half_amp,V_rest,peak_amp_rest
showinfo /W = $agraph;
//usercursoradj(agraph);
//wavestats /Z /R=(750,790) atrace
//variable V_rest=V_avg
//usercursoradj(agraph);
	//findpeak /Q /B=0.5 /M=0.005 /R=(xcsr (A, agraph), xcsr (B, agraph)) trace
	//findpeak /Q /B=0.5 /M=0.001 /R=(190,210) trace
	wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) atrace
	
	variable peak_amp= V_max
	variable peak_amp_rest=peak_amp-resting
	//variable peak_amp_rest=peak_amp-(resting)
	variable half_amp= (0.5*peak_amp_rest)+ resting
	//variable half_amp= (0.5*peak_amp_rest)+ resting
	usercursoradj(agraph);
//Zero_Sub()
	FindLevels/D=spikehalfwidth/M=0.05/Q/R=(xcsr (A, agraph), xcsr (B, agraph))/T=0.01atrace half_amp 
	
variable half_amp_dur = spikehalfwidth[1] - spikehalfwidth[0]	
usercursoradj(agraph);
wavestats /Z /R=(xcsr (A, agraph), xcsr (B, agraph)) atrace
variable fAHP1= V_min
variable fAHP=fAHP1-AP_thresh_wave

	Print "peak_amp_rest, peak_dVdt, AP_thresh_wave, half_amp_dur, fAHP"
	Print peak_amp_rest, peak_dVdt, AP_thresh_wave, half_amp_dur, fAHP
	
	setdatafolder root:
	
	end
	
	//----------------------------------------------------------------------
	Function Paired_Pulse_Ratio_slope()
	
	SVAR c1folder = root:c1folder
setdatafolder c1folder
NVAR freq = s_freq
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) slopes, slopes2, PPR
variable interval
prompt interval, "50 100 150 200 250..."
doprompt  "Your choice is", interval
if (interval==50)
variable in=50
elseif (interval==100)
in=100
endif

variable x
prompt x, "slope x axis 1"
doprompt  "Your choice is", x
variable y
prompt y, "slope x axis 2"
doprompt  "Your choice is", y

showinfo /W = $agraph
NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
Print " Indicate the window by placing the cursors to calculate average Vm"
//usercursoradj(agraph);
variable k
for(k=0;k<(numtrace(waves));k = k+1)

	wave twave = waverefindexed(agraph,k, 1)
	 String atrace = stringfromlist(k,waves,";")
	
	curvefit  /N/Q   line  twave(x,y) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	//removefromgraph $ ("fit_" + atrace)
	slopes[k] = slope;
endfor




variable x1= x+in

variable y1= y+in

for(k=0;k<(numtrace(waves));k = k+1)

	wave twave = waverefindexed(agraph,k, 1)
	//atrace = stringfromlist(k,waves,";")
	
	curvefit  /N/Q   line  twave(x1,y1) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	slopes2[k] = slope;
endfor

for(k=0;k<(numtrace(waves));k = k+1)
PPR[K]= slopes2[k]/slopes[k]
endfor

Variable Paired_Pulse_Ratio, time_diff
Paired_Pulse_Ratio= mean(PPR,0,Pi)

Print " The average PPR for this cell is", Paired_Pulse_Ratio
setdatafolder root:

end
//-----------
Function Gradient()
//Setaxis bottom 48,60;

string agraph= winname(0,1),slopename,atrace
//Setaxis bottom 48,60;
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves))  Slopes, regression
showinfo /W = $agraph



NVAR slope = root:k1

variable x =-3.2

variable y = 2.75

variable k
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	atrace = stringfromlist(k,waves,";")
	
	curvefit  /N/Q   line  twave(0,111) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	//removefromgraph $ ("fit_" + atrace)
	slopes[k] = slope;
endfor


//killwaves SSV, peak,TOP1,TOP2, Timeofpeak, Peak1, Peak2, timeofpeak2, TOP3, TOP4

//setdatafolder root:

end


//=-------------------
Function PlateauHW()

SVAR c1folder = root:c1folder
setdatafolder c1folder

string agraph= winname(0,1);
Zero_Sub()
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) peak,half_peak, hw10,hw9,hw8,hw7,hw6,hw5,hw4,hw3,hw2,hw1 
 //SetAxis left -0.004,0.004
 //etAxis bottom 0.45,1
showinfo /W = $agraph



variable k
for(k=0;k<(numtrace(waves));k = k+1)
wave twave = waverefindexed(agraph,k, 1)

variable x=91
variable y=240

Variable z= x+10
Wavestats/Z/R= (z,y)twave	
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw1[k] = spikehalfwidth[0]-x

 x=291
 y=440

z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw2[k] = spikehalfwidth[0]-x
 
 x=491
y=640

z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw3[k] = spikehalfwidth[0]-x
 
 x=691
 y=840
z=x+10
Wavestats/Z/R= (z,y)twave	
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw4[k] = spikehalfwidth[0]-x
 
x=891
 y=1040

	
z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw5[k] = spikehalfwidth[0]-x
 
 
 x=1091
 y=1240

	
z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw6[k] = spikehalfwidth[0]-x

 x=1291
 y=1440

	
z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw7[k] = spikehalfwidth[0]-x
 
x=1490
 y=1640
	
z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
 hw8[k] = spikehalfwidth[0]-x
 
x=1691
 y=1840

	
z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
 hw9[k] = spikehalfwidth[0]-x
 
 x=1891
y=2000
	
z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw10[k] = spikehalfwidth[0]-x
 
 endfor
 
 setdatafolder root:
 end	

//------------------------
Function PlateauHW_filtered()

SVAR c1folder = root:c1folder
setdatafolder c1folder

string agraph= winname(0,1);
Zero_Sub()
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) peak,half_peak, hw10,hw9,hw8,hw7,hw6,hw5,hw4,hw3,hw2,hw1 
 //SetAxis left -0.004,0.004
 //etAxis bottom 0.45,1
showinfo /W = $agraph



variable k
for(k=0;k<(numtrace(waves));k = k+1)
wave swave = waverefindexed(agraph,k, 1)
variable smth


prompt smth, "do you want to smooth choose 1 the waves"
doprompt  "Your choice is 1 for smth", smth
if (smth==1)
Duplicate /O swave, twave
Smooth /M=0 500, twave

else
Duplicate /O swave, twave
endif


variable int
prompt int, "interval between 2 inputs"
doprompt  "Your choice is", int

variable x=91 +int
variable y=249

Variable z= x+10
Wavestats/Z/R= (z,y)twave	
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw1[k] = spikehalfwidth[0]-x

 x=291+int
 y=449

z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw2[k] = spikehalfwidth[0]-x
 
 x=491+int
y=649

z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw3[k] = spikehalfwidth[0]-x
 
 x=691+int
 y=849
z=x+10
Wavestats/Z/R= (z,y)twave	
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw4[k] = spikehalfwidth[0]-x
 
x=891+int
 y=1049

	
z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw5[k] = spikehalfwidth[0]-x
 
 
 x=1091+int
 y=1249

	
z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw6[k] = spikehalfwidth[0]-x

 x=1291+int
 y=1449

	
z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw7[k] = spikehalfwidth[0]-x
 
x=1490+int
 y=1649
	
z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
 hw8[k] = spikehalfwidth[0]-x
 
x=1691+int
 y=1849

	
z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
 hw9[k] = spikehalfwidth[0]-x
 
 x=1891+int
y=2009
	
z= x+10	
Wavestats/Z/R= (z,y)twave
peak[K]= V_max
half_peak[K]= peak[K]/2
make /O/D/N=2 spikehalfwidth
FindLevels/D=spikehalfwidth/M=0.01/Q/R=(z,y)/T=0.01twave half_peak[K]
	
 hw10[k] = spikehalfwidth[0]-x
 
 endfor
 
 setdatafolder root:
 end	

 //------------------------
 Function PPR_slope()
	
	SVAR c1folder = root:c1folder
setdatafolder c1folder
NVAR freq = s_freq
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) slope1, slope2, PPR

variable x
prompt x, "slope x axis 1"
doprompt  "Your choice is", x


showinfo /W = $agraph
NVAR slope = root:k1
//NVAR timeofpeak = root:V_maxloc
Print " Indicate the window by placing the cursors to calculate average Vm"
//usercursoradj(agraph);


	wave twave = waverefindexed(agraph,0, 1)
	 String atrace = stringfromlist(0,waves,";")
	
	curvefit  /N/Q   line  twave(x,x+1) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	//removefromgraph $ ("fit_" + atrace)
	slope1[0] = slope;
variable x1= x+25
wave twave = waverefindexed(agraph,0, 1)
	 atrace = stringfromlist(0,waves,";")
	
	curvefit  /N/Q   line  twave(x1,x1+1) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	//removefromgraph $ ("fit_" + atrace)
	slope2[0] = slope;

wave twave = waverefindexed(agraph,1, 1)
	atrace = stringfromlist(1,waves,";")
	
	curvefit  /N/Q   line  twave(x,x+1) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	//removefromgraph $ ("fit_" + atrace)
	slope1[1] = slope;
 x1= x+50
wave twave = waverefindexed(agraph,1, 1)
	  atrace = stringfromlist(1,waves,";")
	
	curvefit  /N/Q   line  twave(x1,x1+1) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	//removefromgraph $ ("fit_" + atrace)
	slope2[1] = slope;
	
	
	
wave twave = waverefindexed(agraph,2, 1)
	 atrace = stringfromlist(2,waves,";")
	
	curvefit  /N/Q   line  twave(x,x+1) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	//removefromgraph $ ("fit_" + atrace)
	slope1[2] = slope;
 x1= x+75
wave twave = waverefindexed(agraph,2, 1)
	 atrace = stringfromlist(2,waves,";")
	
	curvefit  /N/Q   line  twave(x1,x1+1) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	//removefromgraph $ ("fit_" + atrace)
	slope2[2] = slope;
	
	
	wave twave = waverefindexed(agraph,3, 1)
	 atrace = stringfromlist(3,waves,";")
	
	curvefit  /N/Q   line  twave(x,x+1) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	//removefromgraph $ ("fit_" + atrace)
	slope1[3] = slope;
x1= x+100
wave twave = waverefindexed(agraph,3, 1)
	atrace = stringfromlist(3,waves,";")
	
	curvefit  /N/Q   line  twave(x1,x1+1) /D
	modifygraph rgb($ ("fit_" + atrace))=(0,0,0);
	//removefromgraph $ ("fit_" + atrace)
	slope2[3] = slope;

PPR= slope2/slope1
Variable Paired_Pulse_Ratio, time_diff
Paired_Pulse_Ratio= mean(PPR,0,Pi)

Print " The average PPR for this cell is", Paired_Pulse_Ratio
setdatafolder root:

end

//----------
Function single_EPSC()


 string agraph= winname(0,1)
 String trace= wavename(agraph,0,1)
wave atrace = waverefindexed(agraph,0,1)
make /O/N =3 Output, Output2

//SVAR c1folder = root:c1folder

//setdatafolder c1folder


Zero_Sub()
//showinfo /W = $agraph;

//usercursoradj(agraph);

//variable x =xcsr(A,a

//variable y = 4300
//prompt x, "window"
//doprompt  "Your choice is", x

 //Display W_DIF vs $atrace
 
WaveStats/Q/R=(183,300) atrace
//WaveStats/Q/R=(500,530) W_DIF							/// finds maximum rate of rise of the first spike
	Variable peak=V_min
	variable half_amp1= (0.5*peak)
	make /O/D/N=2 spikehalfwidth1
	FindLevels/D=spikehalfwidth1/M=0.05/Q/R=(183,300)/T=0.01atrace half_amp1
	
variable half_amp_dur1 = spikehalfwidth1[1] - spikehalfwidth1[0]	

//atrace= atrace/peak
//WaveStats/Q/R=(x,y) atrace						/// finds maximum rate of rise of the first spike
//	Variable npeak= V_max
//variable peak_width, peak_amp,Half_amp,V_rest,peak_amp_rest
//showinfo /W = $agraph;
variable EPSCarea = area (atrace,183,300)
	
//	variable half_amp= (0.5*npeak)
//	make /O/D/N=2 spikehalfwidth
//	FindLevels/D=spikehalfwidth/M=0.05/Q/R=(x,x+50)/T=0.01atrace half_amp 
	
//variable half_amp_dur = spikehalfwidth[1] - spikehalfwidth[0]	

	Print peak,  half_amp_dur1, EPSCarea
	output [0]= peak
	output [1]= half_amp_dur1
	output [2] = EPSCarea
	
	
	// EPSC 2
	WaveStats/Q/R=(4183,4300) atrace
//WaveStats/Q/R=(500,530) W_DIF							/// finds maximum rate of rise of the first spike
	 peak=V_min
	 half_amp1= (0.5*peak)
	make /O/D/N=2 spikehalfwidth1
	FindLevels/D=spikehalfwidth1/M=0.05/Q/R=(4183,4300)/T=0.01atrace half_amp1
	
 half_amp_dur1 = spikehalfwidth1[1] - spikehalfwidth1[0]	

//atrace= atrace/peak
//WaveStats/Q/R=(x,y) atrace						/// finds maximum rate of rise of the first spike
//	Variable npeak= V_max
//variable peak_width, peak_amp,Half_amp,V_rest,peak_amp_rest
//showinfo /W = $agraph;
 EPSCarea = area (atrace,4183,4300)
	
//	variable half_amp= (0.5*npeak)
//	make /O/D/N=2 spikehalfwidth
//	FindLevels/D=spikehalfwidth/M=0.05/Q/R=(x,x+50)/T=0.01atrace half_amp 
	
//variable half_amp_dur = spikehalfwidth[1] - spikehalfwidth[0]	

	Print peak,  half_amp_dur1, EPSCarea
	output2 [0]= peak
	output2 [1]= half_amp_dur1
	output2 [2] = EPSCarea
	
	
	//setdatafolder root:
	
	end
	
	
//---------------------------------------------------------------------------//


//----------
Function single_IPSC()


 string agraph= winname(0,1)
 String trace= wavename(agraph,0,1)
wave atrace = waverefindexed(agraph,0,1)
make /O/N =3 Output, Output2

//SVAR c1folder = root:c1folder

//setdatafolder c1folder


Zero_Sub()
//showinfo /W = $agraph;

//usercursoradj(agraph);

//variable x =xcsr(A,a

//variable y = 4300
//prompt x, "window"
//doprompt  "Your choice is", x

 //Display W_DIF vs $atrace
 
WaveStats/Q/R=(183,600) atrace
//WaveStats/Q/R=(500,530) W_DIF							/// finds maximum rate of rise of the first spike
	Variable peak=V_max
	variable half_amp1= (0.5*peak)
	make /O/D/N=2 spikehalfwidth1
	FindLevels/D=spikehalfwidth1/M=0.05/Q/R=(183,600)/T=0.01atrace half_amp1
	
variable half_amp_dur1 = spikehalfwidth1[1] - spikehalfwidth1[0]	

//atrace= atrace/peak
//WaveStats/Q/R=(x,y) atrace						/// finds maximum rate of rise of the first spike
//	Variable npeak= V_max
//variable peak_width, peak_amp,Half_amp,V_rest,peak_amp_rest
//showinfo /W = $agraph;
variable EPSCarea = area (atrace,183,600)
	
//	variable half_amp= (0.5*npeak)
//	make /O/D/N=2 spikehalfwidth
//	FindLevels/D=spikehalfwidth/M=0.05/Q/R=(x,x+50)/T=0.01atrace half_amp 
	
//variable half_amp_dur = spikehalfwidth[1] - spikehalfwidth[0]	

	Print peak,  half_amp_dur1, EPSCarea
	output [0]= peak
	output [1]= half_amp_dur1
	output [2] = EPSCarea
	
	
	// EPSC 2
	WaveStats/Q/R=(4183,4600) atrace
//WaveStats/Q/R=(500,530) W_DIF							/// finds maximum rate of rise of the first spike
	 peak=V_max
	 half_amp1= (0.5*peak)
	make /O/D/N=2 spikehalfwidth1
	FindLevels/D=spikehalfwidth1/M=0.05/Q/R=(4183,4600)/T=0.01atrace half_amp1
	
 half_amp_dur1 = spikehalfwidth1[1] - spikehalfwidth1[0]	

//atrace= atrace/peak
//WaveStats/Q/R=(x,y) atrace						/// finds maximum rate of rise of the first spike
//	Variable npeak= V_max
//variable peak_width, peak_amp,Half_amp,V_rest,peak_amp_rest
//showinfo /W = $agraph;
 EPSCarea = area (atrace,4183,4600)
	
//	variable half_amp= (0.5*npeak)
//	make /O/D/N=2 spikehalfwidth
//	FindLevels/D=spikehalfwidth/M=0.05/Q/R=(x,x+50)/T=0.01atrace half_amp 
	
//variable half_amp_dur = spikehalfwidth[1] - spikehalfwidth[0]	

	Print peak,  half_amp_dur1, EPSCarea
	output2 [0]= peak
	output2 [1]= half_amp_dur1
	output2 [2] = EPSCarea
	
	
	//setdatafolder root:
	
	end
	
	
//---------------------------------------------------------------------------//


Function IPSCtrains()

//use this program to measure the amplitude and area under the curve of the Chr2 induced IPSC trains

// written for Pfizer drug experiments

	string agraph= winname(0,1)
 String trace= wavename(agraph,0,1)
wave atrace = waverefindexed(agraph,0,1)

Zero_Sub()
//remove300(0.1)
//Averagetraces()

showinfo /W = $agraph

variable int
prompt int, "interval between the stim pules, 100 ms, 50, 25 ms"
doprompt  "Your choice is", int

variable num =4
//prompt num, "number of stim pulses, usually 10"
//doprompt  "Your choice is",num
if (int==100)
string wavename = "curr_peak100"
elseif (int==50)
wavename = "curr_peak50"
else
wavename = "curr_peak25"
endif

make /O/D/N= (num) $wavename, curr_peak, curr_area

variable k
for (k=0; k<num; k=k+1)
 
variable x=102 +int*k
variable y= 102 + int*(k+1)



Wavestats/Z/R= (x,y)atrace
curr_peak[K]= V_max

//curr_area[k] = cur_ar
curr_area[k]= area (atrace,x,y)
endfor

Duplicate/O curr_peak, $wavename; //KillWaves data
end


//----------------------------------------
Function IPSCtimecourse()


 string agraph= winname(0,1)
// string waves = tracenamelist(agraph,";",1);
 String twave= wavename(agraph,0,1)
make /O/D/N= 1 peak, EPSCarea, risetime, decaytime

Zero_Sub()
//remove300()
//Averagetraces()

SetAxis bottom 115, 160
 
WaveStats/Q/R=(100,600) $twave

	 peak=V_max
	 variable x = V_maxloc
	 EPSCarea = area ($twave,121,200)
	 	
//half_amp1[k]= (0.5*peak[k])
//	make /O/D/N=2 spikehalfwidth1
	//FindLevels/D=spikehalfwidth1/M=0.05/Q/R=(100,300)/T=0.01twave half_amp1[k]
	
//halfwidth[k] = spikehalfwidth1[1] - spikehalfwidth1[0]	

//WaveStats/Q/R=(500,530) W_DIF							/// finds maximum rate of rise of the first spike
	usercursoradj(agraph);

variable a =x2pnt($twave,xcsr(A))
variable b =x2pnt($twave,xcsr(B))

	
	CurveFit/M=2/W=0 line, graph_avg[a,b]/D
	 risetime=K1
	 

	variable start = x2pnt($twave, x+2)
	variable stop = x2pnt($twave, x+200)
	

	CurveFit/NTHR=0 /M=2/W=0 exp_XOffset  $twave[start,stop] /D 
	decaytime= k2

End



	
//rserVC()
	

	//end
	
	
//---------------------------------------------------------------------------//
Function rserVC()
string agraph= winname(0,1)
 string waves = tracenamelist(agraph,";",1);
make /O/D/N= (numtrace(waves)) peakr, rse

Zero_Sub()

variable k
for(k=0;k<(numtrace(waves));k = k+1)
wave twave = waverefindexed(agraph,k, 1)
 
WaveStats/Q/R=(730,740) twave

	 peakr[k]=V_min	
//half_amp1[k]= (0.5*peak)
	//make /O/D/N=2 spikehalfwidth1
	//FindLevels/D=spikehalfwidth1/M=0.05/Q/R=(100,300)/T=0.01twave half_amp1
	
//halfwidth[k] = spikehalfwidth1[1] - spikehalfwidth1[0]	


//EPSCarea[k] = area (twave,100,400)
	rse[k]=10e3/(peakr[k]*-1)

	
endfor
	end
	
	//------------------------------------------------------------------
	
Function single_EPSPIPSP()


 string agraph= winname(0,1)
 String trace= wavename(agraph,0,1)
wave atrace = waverefindexed(agraph,0,1)
make /O/N =6 Output

//SVAR c1folder = root:c1folder

//setdatafolder c1folder

//Duplicate atrace, basewave
//basewave = 0

ModifyGraph grid=2

//adjustbaseline2(50,700)
Zero_Sub()

//showinfo /W = $agraph;
ModifyGraph grid=1
usercursoradj(agraph);

variable x =xcsr(A)
variable y =xcsr(B)


//variable y = 4300
//prompt x, "window"
//doprompt  "Your choice is", x

 //Display W_DIF vs $atrace
 
WaveStats/Q/R=(x,y) atrace
//WaveStats/Q/R=(500,530) W_DIF							/// finds maximum rate of rise of the first spike
	Variable EPSPpeak=V_max	
	variable half_amp1= (0.5*EPSPpeak)
	make /O/D/N=2 spikehalfwidth1
	FindLevels/D=spikehalfwidth1/M=10/Q/R=(x,y)/T=0.01atrace half_amp1
	
variable EPSPhw = spikehalfwidth1[1] - spikehalfwidth1[0]	

//atrace= atrace/peak
//WaveStats/Q/R=(x,y) atrace						/// finds maximum rate of rise of the first spike
//	Variable npeak= V_max
//variable peak_width, peak_amp,Half_amp,V_rest,peak_amp_rest
//showinfo /W = $agraph;
variable EPSParea = area (atrace,x,y)
	
//	variable half_amp= (0.5*npeak)
//	make /O/D/N=2 spikehalfwidth
//	FindLevels/D=spikehalfwidth/M=0.05/Q/R=(x,x+50)/T=0.01atrace half_amp 
	
//variable half_amp_dur = spikehalfwidth[1] - spikehalfwidth[0]	

	Print EPSPpeak,  EPSPhw, EPSParea
	output [0]= EPSPpeak
	output [1]= EPSPhw
	output [2]= EPSParea
	
	usercursoradj(agraph);

variable z =xcsr(A)
 y =xcsr(B)
	
	WaveStats/Q/R=(y,z) atrace
//WaveStats/Q/R=(500,530) W_DIF							/// finds maximum rate of rise of the first spike
	Variable IPSPpeak=V_min
	variable half_amp2= (0.5*IPSPpeak)
	make /O/D/N=2 spikehalfwidth2
	FindLevels/D=spikehalfwidth2/M=100/Q/R=(y,z)/T=0.01 atrace half_amp2
	
variable IPSPhw = spikehalfwidth2[1] - spikehalfwidth2[0]	

//atrace= atrace/peak
//WaveStats/Q/R=(x,y) atrace						/// finds maximum rate of rise of the first spike
//	Variable npeak= V_max
//variable peak_width, peak_amp,Half_amp,V_rest,peak_amp_rest
//showinfo /W = $agraph;
variable IPSParea = area (atrace,y, z)
	
	Print IPSPpeak,  IPSPhw, IPSParea
	output [3]= IPSPpeak
	output [4]= IPSPhw
	output [5]= IPSParea
	
	//setdatafolder root:
	
	end
	
	
//---------------------------------------------------------------------------//

Function normalizePSC()
string agraph= winname(0,1);
string waves = tracenamelist(agraph,";",1); 
make /O/D/N= (numtrace(waves)) slopes
showinfo /W = $agraph

variable x
		prompt x, "window"
doprompt  "Your choice is -1 for EPSC, +1 for IPSC", x
Print " Place the cursors inoder to adjust baseline"
usercursoradj(agraph);
variable k,firstavg,nextavg,nextavgch
for(k=0;k<(numtrace(waves));k = k+1)
	wave twave = waverefindexed(agraph,k, 1)
	wavestats /Q/Z /R=(xcsr (A, agraph), xcsr (B, agraph)) twave
	//wavestats /Q/Z /R=(0,80) twave
	
		firstavg = V_avg;
	
		twave = twave/firstavg
		

twave = x*twave
		

endfor
end

//--------------------------------
Function LRG_ThetaGamma(int)

variable int

// written for LRG experiments 032320. RM

string agraph= winname(0,1)

string waves = tracenamelist(agraph,";",1); 

Zero_Sub()


//prompt int, "interval between the stim pules, 25 ms or 12.5 ms"
//doprompt  "Your choice is", int

//make /O/D/N = (5,5) avg_peaks


variable m
for(m=0;m<(numtrace(waves));m = m+1)
	String trace= wavename(agraph,m,1)
	wave atrace = waverefindexed(agraph,m, 1)

//String trace= wavename(agraph,0,1)
//wave atrace = waverefindexed(agraph,0,1)

//Wavestats/Z/R= (0,100)atrace

//Variable baseline = V_avg
//print baseline 

//remove300(0.1)
//Averagetraces()

	duplicate/O atrace, postrace
	postrace = atrace >= 0 ? atrace : 0

	showinfo /W = $agraph
	
	variable bursts = 5


	variable num =5
	variable k


//prompt num, "number of stim pulses, usually 10"
////doprompt  "Your choice is",num


//if (int==25)
//string wavename = "curr_peak40Hz"
//elseif (int==12.5)
//wavename = "curr_peak80Hz"
//endif

	make /O/D/N= (num) curr_peak, curr_area, totalarea, decayarea

	

	variable l

	for (l=0; l<bursts; l=l+1) 

		variable startx = 202
		startx = startx + 200*l
			for (k=0; k<num; k=k+1)
 
			variable x=startx +int*k
			variable y= startx + int*(k+1)

			Wavestats/Z/R= (x,y)atrace
			curr_peak[K]= V_max
			
			

			//curr_area[k] = cur_ar
			Wavestats/Z/R= (x,y)postrace
			curr_area[k]= area (postrace,x,y) 
			
			if (k==4)
			variable peak_loc = V_maxloc
			endif
			
			
			if (l==4)
			variable end_loc = startx + 600
			else
			end_loc = startx + 190
			endif
			
			
			decayarea[l] = area (postrace,peak_loc,end_loc) 
			
			totalarea[l] = area(postrace,startx,end_loc)
			endfor


	string newwavename1 = "curr_peak" + num2str(l)
	string newwavename2 = "curr_area" + num2str(l)
	Duplicate/O curr_peak, $newwavename1 //KillWaves data
	Duplicate/O curr_area, $newwavename2 //KillWaves data
	
	
			
//			CurveFit/NTHR=0 exp_XOffset posttrace[peak_loc,end_loc] /D 
//		
//			decaytime[l]= k2
			
				
	endfor
	
	
	
string concatwavename = "EPSPpeak_out" + num2str(m)
Concatenate/O {curr_peak0,curr_peak1,curr_peak2,curr_peak3,curr_peak4},$concatwavename
killwaves curr_peak0 curr_peak1 curr_peak2 curr_peak3 curr_peak4 curr_peak

//Duplicate/O $concatwavename, old_peak


//avg_peaks = avg_peaks + old_peak



string newwavename3 = "plat_area" + num2str(m)
	Duplicate/O decayarea, $newwavename3 //KillWaves data
	string concatwavename3 = "platarea_out"
	Concatenate/O {plat_area0,plat_area1,plat_area2,plat_area3, plat_area4},$concatwavename3
	killwaves plat_area0 plat_area1 plat_area2 plat_area3 plat_area4 decayarea


	string newwavename4 = "burst_area" + num2str(m)
	Duplicate/O totalarea, $newwavename4 //KillWaves data
	string concatwavename4 = "burstarea_out"
	Concatenate/O {burst_area0,burst_area1,burst_area2,burst_area3, burst_area4},$concatwavename4
	killwaves burst_area0 burst_area1 burst_area2 burst_area3 burst_area4 totalarea

string concatwavename2 = "EPSParea_out" + num2str(m)

Concatenate/O {curr_area0,curr_area1,curr_area2,curr_area3, curr_area4},$concatwavename2
killwaves curr_area0 curr_area1 curr_area2 curr_area3 curr_area4 curr_area

endfor

make /O/D/N = (5,5) avg_peaks
Duplicate/O EPSPpeak_out0, wq0
Duplicate/O EPSPpeak_out1, wq1
Duplicate/O EPSPpeak_out2, wq2
Duplicate/O EPSPpeak_out3, wq3
Duplicate/O EPSPpeak_out4, wq4

avg_peaks = (wq0 + wq1 + wq2 + wq3 + wq4 )/5

killwaves wq0 wq1 wq2 wq3 wq4

//avg_peaks = EPSPpeak_out + EPSPpeak_out1 + EPSPpeak_out2 + EPSPpeak_out3 + EPSPpeak_out4
//avg_peaks = avg_peaks/5

end

//--------------------
Function TG_spiking(int)
variable int

string agraph= winname(0,1)
string waves = tracenamelist(agraph,";",1); 

//Zero_Sub()


//prompt int, "interval between the stim pules, 25 ms or 12.5 ms"
//doprompt  "Your choice is", int


variable m
for(m=0;m<(numtrace(waves));m = m+1)
	String trace= wavename(agraph,m,1)
	wave atrace = waverefindexed(agraph,m, 1)


//duplicate/O atrace, postrace
//postrace = atrace >= 0 ? atrace : 0

showinfo /W = $agraph

make /O /D /N= (20)  peaktime

variable bursts = 5


variable num =5
variable k


make /O/D/N= (num) curr_peak, curr_area

variable l

for (l=0; l<bursts; l=l+1) 

variable startx = 202
startx = startx + 200*l
//for (k=0; k<num; k=k+1)
 
 variable starty = 202 + 5*int
 starty = starty + 200*l
//variable x=startx +int*k
//variable y= startx + int*(k+1)

variable i = 0

do
	findpeak /Q /B=0.0001/M=10/R=(startx,starty) atrace
	if (V_flag != 0)  
		break
	endif
	
	i = i+1;
	peaktime[i] = V_peakLoc
	
	startx = V_peakLoc+1;
	
while(1)
variable AP=(i)
variable freq1 = AP*1000/(int*5)

make /O/D/N= (bursts) Freq, APnum
Freq [l] = freq1
APnum [l] = AP


string newwavename1 = "freq_out" + num2str(m)
string newwavename2 = "AP_out"+ num2str(m)
Duplicate/O Freq, $newwavename1 //KillWaves data
Duplicate/O APnum, $newwavename2 //KillWaves data
endfor

string concatwavename = "outFreq" 
Concatenate/O {freq_out0,freq_out1,freq_out2,freq_out3, freq_out4},$concatwavename
killwaves freq_out0 freq_out1 freq_out2 freq_out3 freq_out4

string concatwavename2 = "outAPnum" 
Concatenate/O {AP_out0,AP_out1,AP_out2,AP_out3, AP_out4},$concatwavename2
killwaves AP_out0 AP_out1 AP_out2 AP_out3 AP_out4

endfor
end


//----------------------------------------------------------------------//
Function smooth_waves()
string agraph= winname(0,1)

string waves = tracenamelist(agraph,";",1); 

Zero_Sub()
variable m
for(m=0;m<(numtrace(waves));m = m+1)
	String trace= wavename(agraph,m,1)
	wave atrace = waverefindexed(agraph,m, 1)
	string new_wavename = trace + "_smth"
	
Duplicate/O atrace, $new_wavename 

Smooth/M=0 100,$new_wavename
if (m==0)
Display $new_wavename
else
AppendToGraph $new_wavename
endif

endfor

end

//----------------------------------------------------------------------------//

Function runall()

variable int
prompt int, "interval between the stim pules, 25 ms or 12.5 ms"
doprompt  "Your choice is", int

TG_spiking(int)

sleep/s 1

smooth_waves()

LRG_ThetaGamma(int)



sleep/s 2

make /O/D/N = (3,5) fill1
make /O/D/N = (2,5) fill2
make /O/D/N = (4,5) fill3
Concatenate/O/NP=0 {outFreq,fill1, burstarea_out, fill2, platarea_out, fill3, avg_peaks},fin_result
killwaves fill1 fill2 fill3

end