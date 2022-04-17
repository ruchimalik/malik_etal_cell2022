#pragma rtGlobals=1		// Use modern global access method.
// written by Rick Gray, Center for Learning and Memory, UT Austin
//modified by Ruchi Malik, UCSF

//loads the older and integer format of abf files.
// description for head info taken from
// http://www.mathworks.com/matlabcentral/fileexchange/6190-abfload

// loads data from all acquisition channels and does not correct the gain, 
//use Duals_selectingchannels.ipf for selecting channels in 2 channel recordings
//corrects gain  

macro go()
	define_vars()
	main_panel()
	
end

function define_vars()
NewPath/O Desktop "RuchiHD480:Users:Ruchi:" // set the path to the data folder

execute "SetDataFolder root: "
	NewPath/O/M="Select Data Directory" desktop  // set the path variable to the data
	
	string/G MyHome=""  // initialize a global string to hold the pathname of your home directory
	pathinfo  Desktop  // Igor makes a path called "Desktop" when started; it expands to /YourHDname/Users/yourname/Desktop
	MyHome = ParseFilePath(1, S_path, ":", 0, 3) // get the first 3 parts of that path
	string/G dataroot= MyHome+"IgorData:" 	// used to make datapath for saving traces, logs, parameter files

	variable/g refNum
	variable/g save_traces
	string/g cell_num = "A"
	string now=Secs2Date(DateTime,-1)
	string/g dirname= now[3,4]+now[0,1]+now[8,9]+cell_num
	NewPath/C/O/Z datapath dataroot+dirname
	variable/g run
	string/g runstring
	variable/g trace
	string/g tracestring
	variable/g plot_waves
	string/g file_prefix
	string/g ABF_filename
	
	NewNotebook/N=filelogs /F=1/V=1/K=0/W=(549,196,1126,849) 
//SaveNotebook/P=desktop filelogs
				  										  // name, position in header in bytes, type, value
	variable/g fFileVersionNumber      //  'fFileVersionNumber',4,'float',-1;
	variable/g nOperationMode          //  'nOperationMode',8,'int16',-1;
	variable/g lActualAcqLength        //  'lActualAcqLength',10,'int32',-1;
	variable/g nNumPointsIgnored       //  'nNumPointsIgnored',14,'int16',-1;
	variable/g lActualEpisodes         //  'lActualEpisodes',16,'int32',-1;
	variable/g lFileStartTime          //  'lFileStartTime',24,'int32',-1;
	variable/g lDataSectionPtr         //  'lDataSectionPtr',40,'int32',-1;
	variable/g lSynchArrayPtr          //  'lSynchArrayPtr',92,'int32',-1;
	variable/g lSynchArraySize         //  'lSynchArraySize',96,'int32',-1;
	variable/g nDataFormat             //  'nDataFormat',100,'int16',-1;
	variable/g nADCNumChannels         //  'nADCNumChannels', 120, 'int16', -1;
	variable/g fADCSampleInterval      //  'fADCSampleInterval',122,'float', -1;
	variable/g fSynchTimeUnit          //  'fSynchTimeUnit',130,'float',-1;
	variable/g lNumSamplesPerEpisode   //  'lNumSamplesPerEpisode',138,'int32',-1;
	variable/g lPreTriggerSamples      //  'lPreTriggerSamples',142,'int32',-1;
	variable/g lEpisodesPerRun         //  'lEpisodesPerRun',146,'int32',-1;
	variable/g fADCRange               //  'fADCRange', 244, 'float', -1;
	variable/g lADCResolution          //  'lADCResolution', 252, 'int32', -1;
	variable/g nFileStartMillisecs     //  'nFileStartMillisecs', 366, 'int16', -1;
	variable/g nADCPtoLChannelMap      //  'nADCPtoLChannelMap', 378, 'int16', tmp;
	variable/g nADCSamplingSeq         //  'nADCSamplingSeq', 410, 'int16',  tmp;
	variable/g sADCChannelName         //  'sADCChannelName',442, 'uchar', repmat(tmp,1,10);
	variable/g fADCProgrammableGain    //  'fADCProgrammableGain', 730, 'float', tmp;
	variable/g fInstrumentScaleFactor  //  'fInstrumentScaleFactor', 922, 'float', tmp;
	variable/g fInstrumentOffset       //  'fInstrumentOffset', 986, 'float', tmp;
	variable/g fSignalGain             //  'fSignalGain', 1050, 'float', tmp;
	variable/g fSignalOffset           //  'fSignalOffset', 1114, 'float', tmp;
	variable/g nTelegraphEnable        //  'nTelegraphEnable',4512,'int16',tmp;
	variable/g fTelegraphAdditGain     //  'fTelegraphAdditGain',4576,'float',tmp

end

	Function AppendToLog(nb, str)
String nb // name of the notebook to log to
String str // the string to log
//Variable stampDateTime // non-zero if we want to include stamp
Silent 1
//Variable/D now
//String stamp
Notebook $nb  selection={endOfFile, endOfFile}
//if (stampDateTime)
//now = datetime
//stamp = Secs2Date(now,0) + ", " + Secs2Time(now,0) + "\r"
//Notebook $nb text=stamp
//endif
Notebook $nb text= str //+"\r"
End
	
	
Function ButtonProc(ctrlName) : ButtonControl
	String ctrlName
	NVAR refNum
	SVAR ABF_filename
	String fileFilters = "Data Files (*.abf):.abf;"
	Open/R/F=fileFilters/P=Desktop refNum
	ABF_filename = S_filename
	print ABF_filename, "opened"
	//Notebook filelogs selection={endOfFile, endOfFile}
	Appendtolog("filelogs",ABF_filename)
	SaveNotebook/P=desktop filelogs
End

Function CheckProc(ctrlName,checked) : CheckBoxControl
	String ctrlName
	Variable checked
	NVAR save_traces
	save_traces = checked
End

Function SetVarProc(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	SVAR file_prefix
	NVAR run
	SVAR cell_num
	SVAR dirname
	SVAR dataroot
	file_prefix = varStr
	cell_num = file_prefix[6]
	run=0  // reset run number when a new folder is created
	NewPath/C/O/Z datapath dataroot+dirname

End

Function ButtonProc_readfile(ctrlName) : ButtonControl
	String ctrlName
	NVAR refNum
	NVAR lNumSamplesPerEpisode
	NVAR lActualEpisodes
	NVAR fADCSampleInterval
	NVAR fADCRange
	NVAR lADCResolution
	SVAR ABF_filename
	NVAR run,trace,save_traces, plot_waves
	SVAR dirname,cell_num,runstring,tracestring
//	Open/R/P=Desktop refNum as "testfile.abf"
//	Open/R/P=Desktop refNum as "2016_11_07_0016a.abf"
///F=f	Controls the number of bytes read and how the bytes are interpreted.
//	f =0:	Native binary format of the object (default).
//	f =1:	Signed one-byte integer.
//	f =2:	Signed 16-bit word; two bytes.
//	f =3:	Signed 32-bit word; four bytes.
//	f =4:	32-bit IEEE floating point; four bytes.
//	f =5:	64-bit IEEE floating point; eight bytes.
//	f =6:	64-bit integer; eight bytes. Requires Igor Pro 7.00 or later.
///U	Integer formats (/F=1, 2, or 3) are unsigned. If /U is omitted, integers are signed.

//go to file position and read in the proper format
//	FSetPos refNum, 4
//	FBinRead/b=3/F=4 refNum, fFileVersionNumber
//	FSetPos refNum, 8
//	FBinRead/F=2 refNum, nOperationMode
	FSetPos refNum, 10
 	FBinRead/F=3 refNum, lActualAcqLength
	FSetPos refNum, 16
	FBinRead/F=3 refNum, lActualEpisodes
	FSetPos refNum, 122
	FBinRead/F=4 refNum, fADCSampleInterval
	FSetPos refNum, 138
	FBinRead/F=3 refNum, lNumSamplesPerEpisode
	FSetPos refNum, 244
	FBinRead/F=4 refNum, fADCRange
	FSetPos refNum, 252
	FBinRead/F=3 refNum, lADCResolution
	FSetPos refNum, 4576
	FBinRead/F=3 refNum, fTelegraphAdditGain
	FSetPos refNum, 120
	FBinRead/F=3 refNum, nADCNumChannels
	

	FSetPos refNum, 730
	FBinRead/F=3 refNum, fADCProgrammableGain
	FSetPos refNum, 922
	FBinRead/F=3 refNum, fInstrumentScaleFactor
	
	
	


	variable size_to_read = lNumSamplesPerEpisode * lActualEpisodes
	print "file contains ", lActualEpisodes, "traces of",lNumSamplesPerEpisode, "points sampled at ",fADCSampleInterval, "usec"
	trace=0
	string/g wavelist=""
	string cellnumber = "A"
		//prompt cellnumber, "enter cell number"
		//doprompt "A", cellnumber
	do
		sprintf runstring, "%03d", run
		sprintf tracestring, "%03d", trace
		
		//string wavename = cellnumber+"r"+runstring+"t"+tracestring+"c1"
		string wavename = cellnumber+"r"+runstring+"t"+tracestring+"c1"
		//string wavenamenew =cellnumber+"r"+runstring+"t"+tracestring+"c2"
		make/O/n=(lNumSamplesPerEpisode) $wavename
		FSetPos refNum, 6144+(lNumSamplesPerEpisode * 2 *trace)
		FBinRead/F=2 refNum, $wavename
		WAVE tmpwave = $wavename
		tmpwave *= (fADCRange/lADCResolution)*100
		SetScale/P x 0,(fADCSampleInterval/1000),"", tmpwave
		//Duplicate $wavenamenew, $wavename
		
		wavelist += wavename + ";"
		trace += 1
		if (save_traces)
			string notestr = ABF_filename + " trace " +tracestring
			Note  $wavename,  notestr

			save/O/P=datapath $wavename as wavename
		endif
	while (trace < lActualEpisodes)
 		if (plot_waves)
 			display $(StringFromList(0, wavelist))
 				variable twave = 1
 				do 
 					appendtograph $(StringFromList(twave, wavelist))
 					twave +=1
 				while (twave < lActualEpisodes)
 		endif
 		if (save_traces)
 			print "traces saved = ",wavelist
 			
 			
 		endif
	run += 1  // increment run number for next abf file read
	remove300()
	
End

Window main_panel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(508,45,808,245)
	Button button_pick_file,pos={79,6},size={132,25},proc=ButtonProc,title="Choose .abf file"
	Button button_pick_file,fSize=18,fColor=(1,16019,65535)
	CheckBox check_save,pos={96,61},size={132,19},proc=CheckProc,title="Save files to disk"
	CheckBox check_save,fSize=14,value= 0
	SetVariable setvar_file_prefix,pos={21,109},size={266,19},proc=SetVarProc,title="File prefix"
	SetVariable setvar_file_prefix,fSize=12,value= dirname
	Button button_doit,pos={57,137},size={190,56},proc=ButtonProc_readfile,title="Do it"
	Button button_doit,fSize=20,fColor=(0,26214,13293)
	CheckBox check_plot,pos={96,85},size={112,19},proc=CheckProc_1,title="Display waves"
	CheckBox check_plot,fSize=14,value= 0
	
EndMacro

Function CheckProc_1(ctrlName,checked) : CheckBoxControl
	String ctrlName
	Variable checked
	NVAR plot_waves
	plot_waves=checked
	
End

