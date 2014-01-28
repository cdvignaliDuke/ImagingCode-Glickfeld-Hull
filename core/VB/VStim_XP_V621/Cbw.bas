'*********************************************************************
'
'  File: CBW.BAS
'
'  (c) Copyright 1996 - 2004 by Measurement Computing Corp.
'      All Rights Reserved
'
' This file contains the Visual BASIC declarations for all Measurement 
' Computing library commands.   This file should be included in the
' project as a Global Module
'
'***********************************************************************

' Current Revision Number
Global Const CURRENTREVNUM = 5.53
' Error Codes
Global Const NOERRORS = 0
Global Const BADBOARD = 1
Global Const DEADDIGITALDEV = 2
Global Const DEADCOUNTERDEV = 3
Global Const DEADDADEV = 4
Global Const DEADADDEV = 5
Global Const NOTDIGITALCONF = 6
Global Const NOTCOUNTERCONF = 7
Global Const NOTDACONF = 8
Global Const NOTADCONF = 9
Global Const NOTMUXCONF = 10
Global Const BADPORTNUM = 11
Global Const BADCOUNTERDEVNUM = 12
Global Const BADDADEVNUM = 13
Global Const BADSAMPLEMODE = 14
Global Const BADINT = 15
Global Const BADADCHAN = 16
Global Const BADCOUNT = 17
Global Const BADCNTRCONFIG = 18
Global Const BADDAVAL = 19
Global Const BADDACHAN = 20
Global Const ALREADYACTIVE = 22
Global Const PAGEOVERRUN = 23
Global Const BADRATE = 24
Global Const COMPATMODE = 25
Global Const TRIGSTATE = 26
Global Const ADSTATUSHUNG = 27
Global Const TOOFEW = 28
Global Const OVERRUN = 29
Global Const BADRANGE = 30
Global Const NOPROGGAIN = 31
Global Const BADFILENAME = 32
Global Const DISKISFULL = 33
Global Const COMPATWARN = 34
Global Const BADPOINTER = 35
Global Const TOOMANYGAINS = 36
Global Const RATEWARNING = 37
Global Const CONVERTDMA = 38
Global Const DTCONNECTERR = 39
Global Const FORECONTINUOUS = 40
Global Const BADBOARDTYPE = 41
Global Const WRONGDIGCONFIG = 42
Global Const NOTCONFIGURABLE = 43
Global Const BADPORTCONFIG = 44
Global Const BADFIRSTPOINT = 45
Global Const ENDOFFILE = 46
Global Const NOT8254CTR = 47
Global Const NOT9513CTR = 48
Global Const BADTRIGTYPE = 49
Global Const BADTRIGVALUE = 50
Global Const BADOPTION = 52
Global Const BADPRETRIGCOUNT = 53
Global Const BADDIVIDER = 55
Global Const BADSOURCE = 56
Global Const BADCOMPARE = 57
Global Const BADTIMEOFDAY = 58
Global Const BADGATEINTERVAL = 59
Global Const BADGATECNTRL = 60
Global Const BADCOUNTEREDGE = 61
Global Const BADSPCLGATE = 62
Global Const BADRELOAD = 63
Global Const BADRECYCLEFLAG = 64
Global Const BADBCDFLAG = 65
Global Const BADDIRECTION = 66
Global Const BADOUTCONTROL = 67
Global Const BADBITNUMBER = 68
Global Const NONEENABLED = 69
Global Const BADCTRCONTROL = 70
Global Const BADEXPCHAN = 71
Global Const WRONGADRANGE = 72
Global Const OUTOFRANGE = 73
Global Const BADTEMPSCALE = 74
Global Const BADERRCODE = 75
Global Const NOQUEUE = 76
Global Const CONTINUOUSCOUNT = 77
Global Const UNDERRUN = 78
Global Const BADMEMMODE = 79
Global Const FREQOVERRUN = 80
Global Const NOCJCCHAN = 81
Global Const BADCHIPNUM = 82
Global Const DIGNOTENABLED = 83
Global Const CONVERT16BITS = 84
Global Const NOMEMBOARD = 85
Global Const DTACTIVE = 86
Global Const NOTMEMCONF = 87
Global Const ODDCHAN = 88
Global Const CTRNOINIT = 89
Global Const NOT8536CTR = 90
Global Const FREERUNNING = 91
Global Const INTERRUPTED = 92
Global Const NOSELECTORS = 93
Global Const NOBURSTMODE = 94
Global Const NOTWINDOWSFUNC = 95
Global Const NOTSIMULCONF = 96
Global Const EVENODDMISMATCH = 97
Global Const M1RATEWARNING = 98
Global Const NOTRS485 = 99
Global Const NOTDOSFUNC = 100
Global Const RANGEMISMATCH = 101
Global Const CLOCKTOOSLOW = 102
Global Const BADCALFACTORS = 103
Global Const BADCONFIGTYPE = 104
Global Const BADCONFIGITEM = 105
Global Const NOPCMCIABOARD = 106
Global Const NOBACKGROUND = 107
Global Const STRINGTOOSHORT = 108
Global Const CONVERTEXTMEM = 109
Global Const BADEUADD = 110
Global Const DAS16JRRATEWARNING = 111
Global Const DAS08TOOLOWRATE = 112
Global Const AMBIGSENSORTYPE = 114      ' more than one sensor type defined for EXP-GP (obsolete)
Global Const AMBIGSENSORONGP = 114      ' more than one sensor type defined for EXP-GP
Global Const NOSENSORTYPEONGP = 115             ' no sensor type defined for EXP-GP
Global Const NOCONVERSIONNEEDED = 116   ' 12 bit board without chan tags - converted in ISR
Global Const NOEXTCONTINUOUS = 117
Global Const INVALIDPRETRIGCONVERT = 118  ' cbConvertPretirg called after cbPretrigScan failed
Global Const BADCTRREG = 119            ' Bad arg to CLoad for 9513 }
Global Const BADTRIGTHRESHOLD = 120     ' Invalid trigger threshold specified in cbSetTrigger }
Global Const BADPCMSLOTREF = 121     ' Bad PCM Card slot reference
Global Const AMBIGPCMSLOTREF = 122  ' Ambiguous PCM Card slot reference
Global Const BADSENSORTYPE = 123 ' Bad sensor type selected in Instacal
Global Const DELBOARDNOTEXIST = 124     ' tried to delete board number which doesn't exist
Global Const NOBOARDNAMEFILE = 125      ' board name file not found
Global Const CFGFILENOTFOUND = 126      ' configuration file not found
Global Const NOVDDINSTALLED = 127     ' CBUL.386 device driver not installed
Global Const NOWINDOWSMEMORY = 128    ' No Windows memory available
Global Const OUTOFDOSMEMORY = 129    ' No DOS memory available
Global Const OBSOLETEOPTION = 130      ' Option on longer supporeted in cbGetConfig/cbSetConfig
Global Const NOPCMREGKEY = 131      ' No registry entry for this PCMCIA board
Global Const NOCBUL32SYS = 132        ' CBUL32.SYS device driver not installed
Global Const NODMAMEMEMORY = 133                ' No memory for device driver's DMA buffer
Global Const IRQNOTAVAILABLE = 134    ' IRQ in use by another device
Global Const NOT7266CTR = 135         ' This board does not have an LS7266 counter /
Global Const BADQUADRATURE = 136      ' Invalid quadrature specified
Global Const BADCOUNTMODE = 137       ' Invalid counting mode specified
Global Const BADENCODING = 138        ' Invalid data encoding specified
Global Const BADINDEXMODE = 139       ' Invalid index mode specified
Global Const BADINVERTINDEX = 140     ' Invalid invert index specified
Global Const BADFLAGPINS = 141        ' Invalid flag pins specified
Global Const NOCTRSTATUS = 142        ' This board does not support cbCStatus()
Global Const NOGATEALLOWED = 143      ' Gating and indexing not allowed simultaneously
Global Const NOINDEXALLOWED = 144     ' Indexing not allowed in non-quadratue mode
Global Const OPENCONNECTION = 145     ' Temperature input has open connection
Global Const BMCONTINUOUSCOUNT = 146  ' Count must be integer multiple of packetsize for recycle mode.
Global Const BADCALLBACKFUNC = 147    ' Invalid pointer to callback function passed as arg
Global Const MBUSINUSE = 148          ' Metrabus in use
Global Const MBUSNOCTLR = 149         ' MetraBus I/O card has no configured controller card
Global Const BADEVENTTYPE = 150       ' Invalid event type specified for this board.
Global Const ALREADYENABLED = 151     ' An event handler has already been enabled for this event type
Global Const BADEVENTSIZE = 152       ' Invalid event count specified.
Global Const CANTINSTALLEVENT = 153   ' Unable to install event handler
Global Const BADBUFFERSIZE = 154	     ' Buffer is too small for operation
Global Const BADAIMODE = 155          ' Invalid Analog Input Mode (RSE, NRSE, DIFF)
Global Const BADSIGNAL = 156          ' Invalid signal type specified. 
Global Const BADCONNECTION = 157      ' Invalid connection specified. 
Global Const BADINDEX = 158           ' Invalid index specified, or reached end of internal connection list. 
Global Const NOCONNECTION = 159       ' No connection is assigned to specified signal. 
Global Const BADBURSTIOCOUNT = 160    ' Count cannot be greater than the FIFO size for BURSTIO mode
Global Const DEADDEV = 161            ' Device has stopped responding. Please check connections.

Global Const INTERNALERR = 200         ' 200-299 = 16 bit library internal errors
Global Const CANT_LOCK_DMA_BUF = 201   ' DMA buffer could not be locked 
Global Const DMA_IN_USE = 202          ' DMA already controlled by another device
Global Const BAD_MEM_HANDLE = 203      ' Invalid Windows memory handle

Global Const INTERNALERR32 = 300          ' 300-399 = 32 bit library internal errors
Global Const CFG_FILE_READ_FAILURE = 304  ' Error reading from configuration file
Global Const CFG_FILE_WRITE_FAILURE = 305 ' Error writing to configuration file
Global Const CFGFILE_CANT_OPEN = 308      ' Cannot open configuration file
Global Const BAD_RTD_CONVERSION = 325     ' Overflow of RTD conversion
Global Const NO_PCI_BIOS = 326            ' PCI BIOS not present on the PC
Global Const BAD_PCI_INDEX = 327          ' Specified PCI board not detected
Global Const NO_PCI_BOARD = 328           ' Specified PCI board not detected
Global Const CANT_INSTALL_INT = 334       ' Cannot install interrupt handler. IRQ already in use

Global Const PCMCIAERRS = 400             ' 400-499 = PCMCIA errors

' These are the most commonly occurring remapped DOS error codes
Global Const DOSBADFUNC = 501
Global Const DOSFILENOTFOUND = 502
Global Const DOSPATHNOTFOUND = 503
Global Const DOSNOHANDLES = 504
Global Const DOSACCESSDENIED = 505
Global Const DOSINVALIDHANDLE = 506
Global Const DOSNOMEMORY = 507
Global Const DOSBADDRIVE = 515
Global Const DOSTOOMANYFILES = 518
Global Const DOSWRITEPROTECT = 519
Global Const DOSDRIVENOTREADY = 521
Global Const DOSSEEKERROR = 525
Global Const DOSWRITEFAULT = 529
Global Const DOSREADFAULT = 530
Global Const DOSGENERALFAULT = 531

Global Const WIN_CANNOT_ENABLE_INT = 603       ' Cannot enable interrupt. IRQ already in use
Global Const WIN_CANNOT_DISABLE_INT = 605      ' Cannot disable interrupts
Global Const WIN_CANT_PAGE_LOCK_BUFFER = 606   ' Insufficient memory to page lock data buffer
Global Const NO_PCM_CARD = 630                 ' PCM card not detected


' Types of operations or functions
Global Const AIFUNCTION = 1           ' Analog Input Function
Global Const AOFUNCTION = 2           ' Analog Output Function
Global Const DIFUNCTION = 3           ' Digital Input Function
Global Const DOFUNCTION = 4           ' Digital Output Function
Global Const CTRFUNCTION = 5          ' Counter Function


Global Const NotUsed = -1

' Maximum length of error string
Global Const ERRSTRLEN = 256

' Maximum length of board name string
Global Const BOARDNAMELEN = 25


' Status values
Global Const IDLE = 0
Global Const RUNNING = 1


Global Const CBENABLED = 1
Global Const CBDISABLED = 0

Global Const UPDATEIMMEDIATE = 0
Global Const UPDATEONCOMMAND = 1

' Types of error reporting
Global Const DONTPRINT = 0
Global Const PRINTWARNINGS = 1
Global Const PRINTFATAL = 2
Global Const PRINTALL = 3

' Types of error handling
Global Const DONTSTOP = 0
Global Const STOPFATAL = 1
Global Const STOPALL = 2

' Types of digital input ports
Global Const DIGITALOUT = 1
Global Const DIGITALIN = 2

' DT Modes for cbSetDTMode ()
Global Const DTIN = 0
Global Const DTOUT = 2

Global Const FROMHERE = -1
Global Const GETFIRST = -2
Global Const GETNEXT = -3

'  Temperature scales
Global Const CELSIUS = 0
Global Const FAHRENHEIT = 1
Global Const KELVIN = 2
Global Const VOLTS = 4

' Types of digital I/O Ports
Global Const AUXPORT = 1
Global Const FIRSTPORTA = 10
Global Const FIRSTPORTB = 11
Global Const FIRSTPORTCL = 12
Global Const FIRSTPORTCH = 13
Global Const SECONDPORTA = 14
Global Const SECONDPORTB = 15
Global Const SECONDPORTCL = 16
Global Const SECONDPORTCH = 17
Global Const THIRDPORTA = 18
Global Const THIRDPORTB = 19
Global Const THIRDPORTCL = 20
Global Const THIRDPORTCH = 21
Global Const FOURTHPORTA = 22
Global Const FOURTHPORTB = 23
Global Const FOURTHPORTCL = 24
Global Const FOURTHPORTCH = 25
Global Const FIFTHPORTA = 26
Global Const FIFTHPORTB = 27
Global Const FIFTHPORTCL = 28
Global Const FIFTHPORTCH = 29
Global Const SIXTHPORTA = 30
Global Const SIXTHPORTB = 31
Global Const SIXTHPORTCL = 32
Global Const SIXTHPORTCH = 33
Global Const SEVENTHPORTA = 34
Global Const SEVENTHPORTB = 35
Global Const SEVENTHPORTCL = 36
Global Const SEVENTHPORTCH = 37
Global Const EIGHTHPORTA = 38
Global Const EIGHTHPORTB = 39
Global Const EIGHTHPORTCL = 40
Global Const EIGHTHPORTCH = 41


' Selectable A/D Ranges codes
Global Const BIP20VOLTS = 15                    ' Bipolar Ranges (-20 to +20 Volts)
Global Const BIP10VOLTS = 1                     ' -10 to +10 Volts
Global Const BIP5VOLTS = 0                      ' -5 to +5 Volts
Global Const BIP4VOLTS = 16                     ' -4 to +4 Volts
Global Const BIP2PT5VOLTS = 2                   ' -2.5 to +2.5 Volts
Global Const BIP2VOLTS  =14                     ' -2 to +2 Volts
Global Const BIP1PT25VOLTS = 3                  ' -1.25 to +1.25 Volts
Global Const BIP1VOLTS = 4                      ' -1 to +1 Volt
Global Const BIPPT625VOLTS = 5                  ' -0.625 to + 0.625 Volt
Global Const BIPPT5VOLTS = 6                    ' -0.5 to +0.5 Volt
Global Const BIPPT25VOLTS = 12                  ' -0.25 to +0.25 Volt
Global Const BIPPT2VOLTS = 13                   ' -0.2 to +0.2 Volt
Global Const BIPPT1VOLTS = 7                    ' -0.1 to +0.1 Volt
Global Const BIPPT05VOLTS = 8                   ' -0.05 to +0.05 Volt
Global Const BIPPT01VOLTS = 9                   ' -0.01 to +0.01 Volt
Global Const BIPPT005VOLTS = 10                 ' -0.005 to +0.005 Volt
Global Const BIP1PT67VOLTS = 11                 ' -1.67 to +1.67 Volts

Global Const UNI10VOLTS = 100                  ' Unipolar Ranges (0 to 10 Volts)
Global Const UNI5VOLTS = 101                   ' 0 to 5 Volts
Global Const UNI4VOLTS = 114                   ' 0 to 4 Volts
Global Const UNI2PT5VOLTS = 102                ' 0 to 2.5 Volts
Global Const UNI2VOLTS = 103                   ' 0 to 2 Volts
Global Const UNI1PT67VOLTS = 109               ' 0 to 1.67 Volts
Global Const UNI1PT25VOLTS = 104               ' 0 to 1.25 Volts
Global Const UNI1VOLTS = 105                   ' 0 to 1 Volt
Global Const UNIPT5VOLTS = 110                 ' 0 to 0.5 Volt
Global Const UNIPT25VOLTS = 111                ' 0 to 0.25 Volt
Global Const UNIPT2VOLTS = 112                 ' 0 to 0.2 Volt
Global Const UNIPT1VOLTS = 106                 ' 0 to 0.1 Volt
Global Const UNIPT05VOLTS = 113                ' 0 to 0.05 Volt
Global Const UNIPT02VOLTS = 108                ' 0 to 0.02 Volt
Global Const UNIPT01VOLTS = 107                ' 0 to 0.01 Volt



Global Const MA4TO20 = 200                     ' Current Ranges (4 to 20 mA )
Global Const MA2to10 = 201                     ' 2 to 10 mA 
Global Const MA1TO5 = 202                      ' 1 to 5 mA 
Global Const MAPT5TO2PT5 = 203                 ' 0.5 to 2.5 mA 
Global Const MA0TO20 = 204                     ' 0 to 20 mA 

Global Const UNIPOLAR = 300                    ' Unipolar range
Global Const BIPOLAR = 301                     ' Bipolar range

' Types of D/A
Global Const ADDA1 = 0
Global Const ADDA2 = 1

' 8536 counter output 1 control
Global Const NOTLINKED = 0
Global Const GATECTR2 = 1
Global Const TRIGCTR2 = 2
Global Const INCTR2 = 3

' Types of 8254 Counter configurations
Global Const HIGHONLASTCOUNT = 0
Global Const ONESHOT = 1
Global Const RATEGENERATOR = 2
Global Const SQUAREWAVE = 3
Global Const SOFTWARESTROBE = 4
Global Const HARDWARESTROBE = 5

' Where to reload from for 9513 counters
Global Const LOADREG = 0
Global Const LOADANDHOLDREG = 1

' Counter recycle modes
Global Const ONETIME = 0
Global Const RECYCLE = 1

' Direction of counting for 9513 counters
Global Const COUNTDOWN = 0
Global Const COUNTUP = 1

' Types of count detection for 9513 counters
Global Const POSITIVEEDGE = 0
Global Const NEGATIVEEDGE = 1

' Counter output control
Global Const ALWAYSLOW = 0
Global Const HIGHPULSEONTC = 1
Global Const TOGGLEONTC = 2
Global Const DISCONNECTED = 4
Global Const LOWPULSEONTC = 5
Global Const HIGHUNTILTC = 6

' Counter input sources
Global Const TCPREVCTR = 0
Global Const CTRINPUT1 = 1
Global Const CTRINPUT2 = 2
Global Const CTRINPUT3 = 3
Global Const CTRINPUT4 = 4
Global Const CTRINPUT5 = 5
Global Const GATE1 = 6
Global Const GATE2 = 7
Global Const GATE3 = 8
Global Const GATE4 = 9
Global Const GATE5 = 10
Global Const FREQ1 = 11
Global Const FREQ2 = 12
Global Const FREQ3 = 13
Global Const FREQ4 = 14
Global Const FREQ5 = 15
Global Const CTRINPUT6 = 101
Global Const CTRINPUT7 = 102
Global Const CTRINPUT8 = 103
Global Const CTRINPUT9 = 104
Global Const CTRINPUT10 = 105
Global Const GATE6 = 106
Global Const GATE7 = 107
Global Const GATE8 = 108
Global Const GATE9 = 109
Global Const GATE10 = 110
Global Const FREQ6 = 111
Global Const FREQ7 = 112
Global Const FREQ8 = 113
Global Const FREQ9 = 114
Global Const FREQ10 = 115

Global Const CTRINPUT11 = 201
Global Const CTRINPUT12 = 202
Global Const CTRINPUT13 = 203
Global Const CTRINPUT14 = 204
Global Const CTRINPUT15 = 205
Global Const GATE11 = 206
Global Const GATE12 = 207
Global Const GATE13 = 208
Global Const GATE14 = 209
Global Const GATE15 = 210
Global Const FREQ11 = 211
Global Const FREQ12 = 212
Global Const FREQ13 = 213
Global Const FREQ14 = 214
Global Const FREQ15 = 215
Global Const CTRINPUT16 = 301
Global Const CTRINPUT17 = 302
Global Const CTRINPUT18 = 303
Global Const CTRINPUT19 = 304
Global Const CTRINPUT20 = 305
Global Const GATE16 = 306
Global Const GATE17 = 307
Global Const GATE18 = 308
Global Const GATE19 = 309
Global Const GATE20 = 310
Global Const FREQ16 = 311
Global Const FREQ17 = 312
Global Const FREQ18 = 313
Global Const FREQ19 = 314
Global Const FREQ20 = 315

' 9513 Counter registers
Global Const LOADREG1 = 1
Global Const LOADREG2 = 2
Global Const LOADREG3 = 3
Global Const LOADREG4 = 4
Global Const LOADREG5 = 5
Global Const LOADREG6 = 6
Global Const LOADREG7 = 7
Global Const LOADREG8 = 8
Global Const LOADREG9 = 9
Global Const LOADREG10 = 10

Global Const LOADREG11 = 11
Global Const LOADREG12 = 12
Global Const LOADREG13 = 13
Global Const LOADREG14 = 14
Global Const LOADREG15 = 15
Global Const LOADREG16 = 16
Global Const LOADREG17 = 17
Global Const LOADREG18 = 18
Global Const LOADREG19 = 19
Global Const LOADREG20 = 20

Global Const HOLDREG1 = 101
Global Const HOLDREG2 = 102
Global Const HOLDREG3 = 103
Global Const HOLDREG4 = 104
Global Const HOLDREG5 = 105
Global Const HOLDREG6 = 106
Global Const HOLDREG7 = 107
Global Const HOLDREG8 = 108
Global Const HOLDREG9 = 109
Global Const HOLDREG10 = 110

Global Const HOLDREG11 = 111
Global Const HOLDREG12 = 112
Global Const HOLDREG13 = 113
Global Const HOLDREG14 = 114
Global Const HOLDREG15 = 115
Global Const HOLDREG16 = 116
Global Const HOLDREG17 = 117
Global Const HOLDREG18 = 118
Global Const HOLDREG19 = 119
Global Const HOLDREG20 = 120

Global Const ALARM1CHIP1 = 201
Global Const ALARM2CHIP1 = 202
Global Const ALARM1CHIP2 = 301
Global Const ALARM2CHIP2 = 302
Global Const ALARM1CHIP3 = 401
Global Const ALARM2CHIP3 = 402
Global Const ALARM1CHIP4 = 501
Global Const ALARM2CHIP4 = 502

' LS7266 Counter registers
Global Const COUNT1 = 601
Global Const COUNT2 = 602
Global Const COUNT3 = 603
Global Const COUNT4 = 604

Global Const PRESET1 = 701
Global Const PRESET2 = 702
Global Const PRESET3 = 703
Global Const PRESET4 = 704

Global Const PRESCALER1 = 801
Global Const PRESCALER2 = 802
Global Const PRESCALER3 = 803
Global Const PRESCALER4 = 804


'  Counter Gate Control
Global Const NOGATE = 0
Global Const AHLTCPREVCTR = 1
Global Const AHLNEXTGATE = 2
Global Const AHLPREVGATE = 3
Global Const AHLGATE = 4
Global Const ALLGATE = 5
Global Const AHEGATE = 6
Global Const ALEGATE = 7

' 7266 Counter Quadrature values
Global Const NO_QUAD = 0
Global Const X1_QUAD = 1
Global Const X2_QUAD = 2
Global Const X4_QUAD = 4

' 7266 Counter Counting Modes
Global Const NORMAL_MODE = 0
Global Const RANGE_LIMIT = 1
Global Const NO_RECYCLE = 2
Global Const MODULO_N = 3

' 7266 Counter encodings
Global Const BCD_ENCODING = 1
Global Const BINARY_ENCODING = 2

' 7266 Counter Index Modes
Global Const INDEX_DISABLED = 0
Global Const LOAD_CTR = 1
Global Const LOAD_OUT_LATCH = 2
Global Const RESET_CTR = 3

' 7266 Counter Flag Pins
Global Const CARRY_BORROW = 1
Global Const COMPARE_BORROW = 2
Global Const CARRYBORROW_UPDOWN = 3
Global Const INDEX_ERROR = 4

' Counter status bits
Global Const C_UNDERFLOW = &H1
Global Const C_OVERFLOW = &H2
Global Const C_COMPARE = &H4
Global Const C_SIGN = &H8
Global Const C_ERROR = &H10
Global Const C_UP_DOWN = &H20
Global Const C_INDEX = &H40


' Types of triggers
Global Const TRIGABOVE = 0
Global Const TRIGBELOW = 1
Global Const GATENEGHYS = (TRIGBELOW + 1)
Global Const GATEPOSHYS = (TRIGBELOW + 2)
Global Const GATEABOVE = (TRIGBELOW + 3)
Global Const GATEBELOW = (TRIGBELOW + 4)
Global Const GATEINWINDOW = (TRIGBELOW + 5)
Global Const GATEOUTWINDOW = (TRIGBELOW + 6)
Global Const GATEHIGH = (TRIGBELOW + 7)
Global Const GATELOW = (TRIGBELOW + 8)
Global Const TRIGHIGH = (TRIGBELOW + 9)
Global Const TRIGLOW = (TRIGBELOW + 10)
Global Const TRIGPOSEDGE = (TRIGBELOW + 11)
Global Const TRIGNEGEDGE = (TRIGBELOW + 12)


' Signal I/O Configuration Parameters 
' --Connections 
Global Const AUXIN0&  = &H01          
Global Const AUXIN1&  = &H02
Global Const AUXIN2&  = &H04
Global Const AUXIN3&  = &H08
Global Const AUXIN4&  = &H10
Global Const AUXIN5&  = &H20

Global Const AUXOUT0&  = &H0100
Global Const AUXOUT1&  = &H0200
Global Const AUXOUT2& = &H0400

Global Const DS_CONNECTOR& = &H01000

Global Const MAX_CONNECTIONS& = 4  'maximum  number connections per output signal

' --Signal Types 
Global Const ADC_CONVERT& = &H0001    
Global Const ADC_GATE& = &H0002   
Global Const ADC_START_TRIG& = &H0004   
Global Const ADC_STOP_TRIG& = &H0008   
Global Const ADC_TB_SRC& = &H0010
Global Const ADC_SCANCLK& = &H0020&
Global Const ADC_SSH& = &H0040&
Global Const ADC_STARTSCAN& = &H0080&
Global Const ADC_SCAN_STOP& = &H0100&

Global Const DAC_UPDATE& = &H0200&
Global Const DAC_TB_SRC& = &H0400&
Global Const DAC_START_TRIG& = &H0800&

Global Const SYNC_CLK& = &H1000&

Global Const CTR1_CLK& = &H2000&
Global Const CTR2_CLK& = &H4000&

Global Const DGND& = &H8000&

' -- Signal Direction 
Global Const SIGNAL_IN& = 2
Global Const SIGNAL_OUT& = 4

Global Const INVERTED& = 1
Global Const NONINVERTED& = 0




' Types of configuration information
Global Const GLOBALINFO = 1
Global Const BOARDINFO = 2
Global Const DIGITALINFO = 3
Global Const COUNTERINFO = 4
Global Const EXPANSIONINFO = 5
Global Const MISCINFO = 6
Global Const EXPINFOARRAY = 7
Global Const MEMINFO = 8

' Types of global configuration information
Global Const GIVERSION = 36        ' Config file format version number
Global Const GINUMBOARDS = 38      ' Maximum number of boards
Global Const GINUMEXPBOARDS = 40   ' Maximum number of expansion boards

' Types of board configuration information
Global Const BIBASEADR = 0         ' Base Address
Global Const BIBOARDTYPE = 1       ' Board Type (0x101 - 0x7FFF)
Global Const BIINTLEVEL = 2        ' Interrupt level
Global Const BIDMACHAN = 3         ' DMA channel
Global Const BIINITIALIZED = 4     ' TRUE or FALSE
Global Const BICLOCK = 5           ' Clock freq (1, 10 or bus)
Global Const BIRANGE = 6           ' Switch selectable range
Global Const BINUMADCHANS = 7      ' Number of A/D channels
Global Const BIUSESEXPS = 8        ' Supports expansion boards TRUE/FALSE
Global Const BIDINUMDEVS = 9       ' Number of digital devices
Global Const BIDIDEVNUM = 10       ' Index into digital information
Global Const BICINUMDEVS = 11      ' Number of counter devices
Global Const BICIDEVNUM = 12       ' Index into counter information
Global Const BINUMDACHANS = 13     ' Number of D/A channels
Global Const BIWAITSTATE = 14      ' Wait state enabled TRUE/FALSE
Global Const BINUMIOPORTS = 15     ' I/O address space used by board
Global Const BIPARENTBOARD = 16    ' Board number of parent board
Global Const BIDTBOARD = 17        ' Board number of connected DT board
Global Const BINUMEXPS = 18        ' Number of EXPs attached to board
Global Const BISERIALNUM = 214     ' Serial Number for USB boards 
Global Const BIDACUPDATEMODE = 215 ' Update immediately or upon AOUPDATE command
Global Const BIDACUPDATECMD = 216  ' Issue D/A UPDATE command 
Global Const BIDACSTARTUP = 217    ' Restore last value written for startup  
Global Const BIADTRIGCOUNT = 219	  ' Number of samples to acquire per trigger in retrigger mode
Global Const BIADFIFOSIZE = 220    ' Set FIFO override size for retrigger mode
Global Const BIADSOURCE = 221      ' Set A/D source to internal reference(>=0) or external connector(-1)
Global Const BICALOUTPUT  = 222    ' CAL output pin setting 
Global Const BISRCADPACER = 223     ' Source A/D Pacer output 

' Types of digital device information
Global Const DIBASEADR = 0         ' Base address
Global Const DIINITIALIZED = 1     ' TRUE or FALSE
Global Const DIDEVTYPE = 2         ' AUXPORT or xPORTA - CH
Global Const DIMASK = 3            ' Bit mask for this port
Global Const DIREADWRITE = 4       ' Read required before write
Global Const DICONFIG = 5          ' Current configuration
Global Const DINUMBITS = 6         ' Number of bits in port
Global Const DICURVAL = 7          ' Current value of outputs
Global Const DIINMASK = 8		   ' Input bit mask for port
Global Const DIOUTMASK = 9			' Output bit mask for port

' Types of counter device information
Global Const CIBASEADR = 0         ' Base address
Global Const CIINITIALIZED = 1     ' TRUE or FALSE
Global Const CICTRTYPE = 2         ' Counter type 8254, 9513 or 8536
Global Const CICTRNUM = 3          ' Which counter on chip
Global Const CICONFIGBYTE = 4      ' Configuration byte

' Types of expansion board information
Global Const XIBOARDTYPE = 0       ' Expansion board type
Global Const XIMUXADCHAN1 = 1      ' 0 - 15
Global Const XIMUXADCHAN2 = 2      ' 0 - 15 or NOTUSED
Global Const XIRANGE1 = 3          ' Range (gain) of low 16 chans
Global Const XIRANGE2 = 4          ' Range (gain) of high 16 chans
Global Const XICJCCHAN = 5         ' 0 - 15 or NOTUSED
Global Const XITHERMTYPE = 6       ' TYPEJ, TYPEK, TYPEB, TYPET, TYPEE, TYPER, or TYPES
Global Const XINUMEXPCHANS = 7     ' Number of expansion channels on board
Global Const XIPARENTBOARD = 8     ' Board number of parent A/D board
Global Const XISPARE0 = 9          ' 16 words of misc options

' Types of Events
Global Const ON_SCAN_ERROR = &H1&
Global Const ON_EXTERNAL_INTERRUPT = &H2&
Global Const ON_PRETRIGGER = &H4&
Global Const ON_DATA_AVAILABLE = &H8&
Global Const ON_END_OF_AI_SCAN = &H10&
Global Const ON_END_OF_AO_SCAN = &H20&
Global Const ALL_EVENT_TYPES = &HFF&


' If you are using an older version of Visual BASIC that does not recognize the #If
' statement below then remove all of the lines from #If to #Else and also remove
' the "End IF statement below
#If Win32 Then

' Option Flags
'
Global Const FOREGROUND = &H0&
Global Const BACKGROUND = &H1&

Global Const SINGLEEXEC = &H0&
Global Const CONTINUOUS = &H2&

Global Const TIMED = &H0&
Global Const EXTCLOCK = &H4&

Global Const NOCONVERTDATA = &H0&
Global Const CONVERTDATA = &H8&

Global Const NODTCONNECT = &H0&
Global Const DTCONNECT = &H10&

Global Const DEFAULTIO = &H0&
Global Const SINGLEIO = &H20&
Global Const DMAIO = &H40&
Global Const BLOCKIO = &H60&
Global Const BURSTIO = &H10000&    ' Transfer upon scan completion
Global Const RETRIGMODE = &H20000&  ' Re-arm trigger upon acquiring trigger count samples

Global Const BYTEXFER = &H0&
Global Const WORDXFER = &H100&

Global Const INDIVIDUAL = &H0&
Global Const SIMULTANEOUS = &H200&

Global Const FILTER = &H0&
Global Const NOFILTER = &H400&

Global Const NORMMEMORY = &H0&
Global Const EXTMEMORY = &H800&

Global Const BURSTMODE = &H1000&

Global Const NOTODINTS = &H2000&

Global Const EXTTRIGGER = &H4000&

Global Const NOCALIBRATEDATA = &H8000&
Global Const CALIBRATEDATA = &H0&
        '
        ' 32-bit function prototypes
        '
Declare Function cbLoadConfig Lib "cbw32.dll" (ByVal CfgFileName$) As Long
Declare Function cbSaveConfig Lib "cbw32.dll" (ByVal CfgFileName$) As Long
Declare Function cbAConvertData Lib "cbw32.dll" (ByVal BoardNum&, ByVal NumPoints&, ADData%, ChanTags%) As Long
Declare Function cbACalibrateData Lib "cbw32.dll" (ByVal BoardNum&, ByVal NumPoints&, ByVal Gain&, ADData%) As Long
Declare Function cbAConvertPretrigData Lib "cbw32.dll" (ByVal BoardNum&, ByVal PretrigCount&, ByVal TotalCount&, ADData%, ChanTags%) As Long
Declare Function cbAIn Lib "cbw32.dll" (ByVal BoardNum&, ByVal Chan&, ByVal Gain&, DataValue%) As Long
Declare Function cbAInScan Lib "cbw32.dll" (ByVal BoardNum&, ByVal LowChan&, ByVal HighChan&, ByVal CBCount&, CBRate&, ByVal Gain&, ByVal MemHandle&, ByVal Options&) As Long
Declare Function cbALoadQueue Lib "cbw32.dll" (ByVal BoardNum&, ChanArray%, GainArray%, ByVal NumChans&) As Long
Declare Function cbAOut Lib "cbw32.dll" (ByVal BoardNum&, ByVal Chan&, ByVal Gain&, ByVal DataValue%) As Long
Declare Function cbAOutScan Lib "cbw32.dll" (ByVal BoardNum&, ByVal LowChan&, ByVal HighChan&, ByVal CBCount&, CBRate&, ByVal Gain&, ByVal MemHandle&, ByVal Options&) As Long
Declare Function cbAPretrig Lib "cbw32.dll" (ByVal BoardNum&, ByVal LowChan&, ByVal HighChan&, PretrigCount&, CBCount&, CBRate&, ByVal Gain&, ByVal MemHandle&, ByVal Options&) As Long
Declare Function cbATrig Lib "cbw32.dll" (ByVal BoardNum&, ByVal Chan&, ByVal TrigType&, ByVal TrigValue%, ByVal Gain&, DataValue%) As Long
Declare Function cbC7266Config Lib "cbw32.dll" (ByVal BoardNum&, ByVal CounterNum&, ByVal Quadrature&, ByVal CountingMode&, ByVal DataEncoding&, ByVal IndexMode&, ByVal InvertIndex&, ByVal FlagPins&, ByVal GateEnable&) As Long
Declare Function cbC8254Config Lib "cbw32.dll" (ByVal BoardNum&, ByVal CounterNum&, ByVal Config&) As Long
Declare Function cbC8536Config Lib "cbw32.dll" (ByVal BoardNum&, ByVal CounterNum&, ByVal OutputControl&, ByVal RecycleMode&, ByVal Retrigger&) As Long
Declare Function cbC9513Config Lib "cbw32.dll" (ByVal BoardNum&, ByVal CounterNum&, ByVal GateControl&, ByVal CounterEdge&, ByVal CountSource&, ByVal SpecialGate&, ByVal Reload&, ByVal RecycleMode&, ByVal BCDMode&, ByVal CountDirec&, ByVal OutputCtrl&) As Long
Declare Function cbC8536Init Lib "cbw32.dll" (ByVal BoardNum&, ByVal ChipNum&, ByVal Ctr1Output&) As Long
Declare Function cbC9513Init Lib "cbw32.dll" (ByVal BoardNum&, ByVal ChipNum&, ByVal FOutDivider&, ByVal FOutSource&, ByVal Compare1&, ByVal Compare2&, ByVal TimeOfDay&) As Long
Declare Function cbCStoreOnInt Lib "cbw32.dll" (ByVal BoardNum&, ByVal IntCount&, CntrControl%, ByVal DataBuffer&) As Long
Declare Function cbCFreqIn Lib "cbw32.dll" (ByVal BoardNum&, ByVal SigSource&, ByVal GateInterval&, CBCount%, Freq&) As Long

'*****************************************************************************
'   Legacy Function Prototypes: to revert to legacy calls, un-comment the
'          prototypes immediately below.
'
'      Declare Function cbCIn Lib "cbw32.dll" (ByVal BoardNum&, ByVal CounterNum&, CBCount&) As Long
' 
'   Remove the following if using the above legacy function prototypes.
Declare Function cbCIn Lib "cbw32.dll" (ByVal BoardNum&, ByVal CounterNum&, CBCount%) As Long
'*****************************************************************************

Declare Function cbCIn32 Lib "cbw32.dll" (ByVal BoardNum&, ByVal CounterNum&, CBCount&) As Long
Declare Function cbCLoad Lib "cbw32.dll" (ByVal BoardNum&, ByVal RegNum&, ByVal LoadValue&) As Long
Declare Function cbCLoad32 Lib "cbw32.dll" (ByVal BoardNum&, ByVal RegNum&, ByVal LoadValue&) As Long
Declare Function cbCStatus Lib "cbw32.dll" (ByVal BoardNum&, ByVal CounterNum&, StatusBits&) As Long
Declare Function cbDBitIn Lib "cbw32.dll" (ByVal BoardNum&, ByVal PortType&, ByVal BitNum&, BitValue%) As Long
Declare Function cbDBitOut Lib "cbw32.dll" (ByVal BoardNum&, ByVal PortType&, ByVal BitNum&, ByVal BitValue&) As Long
Declare Function cbDConfigPort Lib "cbw32.dll" (ByVal BoardNum&, ByVal PortNum&, ByVal Direction&) As Long
Declare Function cbDConfigBit Lib "cbw32.dll" (ByVal BoardNum&, ByVal PortNum&, ByVal BitNum&, ByVal Direction&) As Long
Declare Function cbDeclareRevision Lib "cbw32.dll" (RevNum!) As Long
Declare Function cbDIn Lib "cbw32.dll" (ByVal BoardNum&, ByVal PortNum&, DataValue%) As Long
Declare Function cbDInScan Lib "cbw32.dll" (ByVal BoardNum&, ByVal PortNum&, ByVal CBCount&, CBRate&, ByVal MemHandle&, ByVal Options&) As Long
Declare Function cbDOut Lib "cbw32.dll" (ByVal BoardNum&, ByVal PortNum&, ByVal DataValue%) As Long
Declare Function cbDOutScan Lib "cbw32.dll" (ByVal BoardNum&, ByVal PortNum&, ByVal CBCount&, CBRate&, ByVal MemHandle&, ByVal Options&) As Long
Declare Function cbErrHandling Lib "cbw32.dll" (ByVal ErrReporting&, ByVal ErrHandling&) As Long
Declare Function cbFileAInScan Lib "cbw32.dll" (ByVal BoardNum&, ByVal LowChan&, ByVal HighChan&, ByVal CBCount&, CBRate&, ByVal Gain&, ByVal FileName$, ByVal Options&) As Long
Declare Function cbFileGetInfo Lib "cbw32.dll" (ByVal FileName$, LowChan%, HighChan%, PretrigCount&, TotalCount&, CBRate&, Gain&) As Long
Declare Function cbFilePretrig Lib "cbw32.dll" (ByVal BoardNum&, ByVal LowChan&, ByVal HighChan&, PretrigCount&, CBCount&, CBRate&, ByVal Gain&, ByVal FileName$, ByVal Options&) As Long
Declare Function cbFileRead Lib "cbw32.dll" (ByVal FileName$, ByVal FirstPoint&, NumPoints&, DataBuffer%) As Long
Declare Function cbFlashLED Lib "cbw32.dll" (ByVal BoardNum&) As Long
Declare Function cbGetErrMsg Lib "cbw32.dll" (ByVal ErrCode&, ByVal ErrMsg$) As Long
Declare Function cbGetRevision Lib "cbw32.dll" (DLLRevNum!, VXDRevNum!) As Long

'*****************************************************************************
'   Legacy Function Prototypes: to revert to legacy calls, un-comment the
'          prototypes immediately below.
'
'       Declare Function cbGetStatus Lib "cbw32.dll" (ByVal BoardNum&, Status%, CurCount&, CurIndex&) As Long
'       Declare Function cbStopBackground Lib "cbw32.dll" (ByVal BoardNum&) As Long
'
'   Remove the following if using the above legacy function prototypes.
Declare Function cbGetStatus Lib "cbw32.dll" Alias "cbGetIOStatus" (ByVal BoardNum&, Status%, CurCount&, CurIndex&, ByVal FunctionType&) As Long
Declare Function cbStopBackground Lib "cbw32.dll" Alias "cbStopIOBackground" (ByVal BoardNum&, ByVal FunctionType&) As Long
'*****************************************************************************

Declare Function cbMemSetDTMode Lib "cbw32.dll" (ByVal BoardNum&, ByVal Mode&) As Long
Declare Function cbMemReset Lib "cbw32.dll" (ByVal BoardNum&) As Long
Declare Function cbMemRead Lib "cbw32.dll" (ByVal BoardNum&, DataBuffer%, ByVal FirstPoint&, ByVal CBCount&) As Long
Declare Function cbMemWrite Lib "cbw32.dll" (ByVal BoardNum&, DataBuffer%, ByVal FirstPoint&, ByVal CBCount&) As Long
Declare Function cbMemReadPretrig Lib "cbw32.dll" (ByVal BoardNum&, DataBuffer%, ByVal FirstPoint&, ByVal CBCount&) As Long
Declare Function cbRS485 Lib "cbw32.dll" (ByVal BoardNum&, ByVal Transmit&, ByVal Receive&) As Long
Declare Function cbTIn Lib "cbw32.dll" (ByVal BoardNum&, ByVal Chan&, ByVal CBScale&, TempValue!, ByVal Options&) As Long
Declare Function cbTInScan Lib "cbw32.dll" (ByVal BoardNum&, ByVal LowChan&, ByVal HighChan&, ByVal CBScale&, DataBuffer!, ByVal Options&) As Long
Declare Function cbWinBufToArray Lib "cbw32.dll" (ByVal MemHandle&, DataBuffer%, ByVal FirstPoint&, ByVal CBCount&) As Long
Declare Function cbWinArrayToBuf Lib "cbw32.dll" (DataBuffer%, ByVal MemHandle&, ByVal FirstPoint&, ByVal CBCount&) As Long
Declare Function cbWinBufAlloc Lib "cbw32.dll" (ByVal NumPoints&) As Long
Declare Function cbWinBufFree Lib "cbw32.dll" (ByVal MemHandle&) As Long
Declare Function cbInByte Lib "cbw32.dll" (ByVal BoardNum&, ByVal PortNum&) As Long
Declare Function cbOutByte Lib "cbw32.dll" (ByVal BoardNum&, ByVal PortNum&, ByVal PortVal%) As Long
Declare Function cbInWord Lib "cbw32.dll" (ByVal BoardNum&, ByVal PortNum&) As Long
Declare Function cbOutWord Lib "cbw32.dll" (ByVal BoardNum&, ByVal PortNum&, ByVal PortVal%) As Long

'*****************************************************************************
'   Legacy Function Prototypes: to revert to legacy calls, un-comment the
'          prototypes immediately below.
'
'      Declare Function cbGetConfig Lib "cbw32.dll" (ByVal InfoType&, ByVal BoardNum&, ByVal DevNum&, ByVal ConfigItem&, ByRef ConfigVal%) As Long
'
'   Remove the following if using the above legacy function prototypes.
Declare Function cbGetConfig Lib "cbw32.dll" (ByVal InfoType&, ByVal BoardNum&, ByVal DevNum&, ByVal ConfigItem&, ByRef ConfigVal&) As Long
'*****************************************************************************

Declare Function cbSetConfig Lib "cbw32.dll" (ByVal InfoType&, ByVal BoardNum&, ByVal DevNum&, ByVal ConfigItem&, ByVal ConfigVal&) As Long
Declare Function cbToEngUnits Lib "cbw32.dll" (ByVal BoardNum&, ByVal Range&, ByVal DataVal%, EngUnits!) As Long
Declare Function cbFromEngUnits Lib "cbw32.dll" (ByVal BoardNum&, ByVal Range&, ByVal EngUnits!, DataVal%) As Long
Declare Function cbGetBoardName Lib "cbw32.dll" (ByVal BoardNum&, ByVal BoardName$) As Long
Declare Function cbSetTrigger Lib "cbw32.dll" (ByVal BoardNum&, ByVal TrigType&, ByVal LowThreshold%, ByVal HighThreshold%) As Long
Declare Function cbEnableEvent Lib "cbw32.dll" (ByVal BoardNum&, ByVal EventType&, ByVal EventSize&, ByVal Callback As Long, ByRef UserData As Any) As Long
Declare Function cbDisableEvent Lib "cbw32.dll" (ByVal BoardNum&, ByVal EventType&) As Long
Declare Function cbSelectSignal Lib "cbw32.dll" (ByVal BoardNum&, ByVal Direction&, ByVal Signal&, ByVal Connection&, ByVal Polarity&) As Long
Declare Function cbGetSignal Lib "cbw32.dll" (ByVal BoardNum&, ByVal Direction&, ByVal Signal&, ByVal Index&, ByRef Connection&, ByRef Polarity&) As Long

#Else
'
'  16-bit
'

' Option Flags
'
Global Const FOREGROUND = &H0
Global Const BACKGROUND = &H1

Global Const SINGLEEXEC = &H0
Global Const CONTINUOUS = &H2

Global Const TIMED = &H0
Global Const EXTCLOCK = &H4

Global Const NOCONVERTDATA = &H0
Global Const CONVERTDATA = &H8

Global Const NODTCONNECT = &H0
Global Const DTCONNECT = &H10

Global Const DEFAULTIO = &H0
Global Const SINGLEIO = &H20
Global Const DMAIO = &H40
Global Const BLOCKIO = &H60

Global Const BYTEXFER = &H0
Global Const WORDXFER = &H100

Global Const INDIVIDUAL = &H0
Global Const SIMULTANEOUS = &H200

Global Const FILTER = &H0
Global Const NOFILTER = &H400

Global Const NORMMEMORY = &H0
Global Const EXTMEMORY = &H800

Global Const BURSTMODE = &H1000

Global Const NOTODINTS = &H2000

Global Const EXTTRIGGER = &H4000

Global Const NOCALIBRATEDATA = &H8000
Global Const CALIBRATEDATA = &H0
        '
        ' 16-bit function prototypes
        '
Declare Function cbAConvertData Lib "cbw.dll" (ByVal BoardNum%, ByVal NumPoints&, ADData%, ChanTags%) As Integer
Declare Function cbACalibrateData Lib "cbw.dll" (ByVal BoardNum%, ByVal NumPoints&, ByVal Gain%, ADData%) As Integer
Declare Function cbAConvertPretrigData Lib "cbw.dll" (ByVal BoardNum%, ByVal PretrigCount&, ByVal TotalCount&, ADData%, ChanTags%) As Integer
Declare Function cbAIn Lib "cbw.dll" (ByVal BoardNum%, ByVal Chan%, ByVal Gain%, DataValue%) As Integer
Declare Function cbAInScan Lib "cbw.dll" (ByVal BoardNum%, ByVal LowChan%, ByVal HighChan%, ByVal CBCount&, CBRate&, ByVal Gain%, ByVal MemHandle%, ByVal Options%) As Integer
Declare Function cbALoadQueue Lib "cbw.dll" (ByVal BoardNum%, ChanArray%, GainArray%, ByVal NumChans%) As Integer
Declare Function cbAOut Lib "cbw.dll" (ByVal BoardNum%, ByVal Chan%, ByVal Gain%, ByVal DataValue%) As Integer
Declare Function cbAOutScan Lib "cbw.dll" (ByVal BoardNum%, ByVal LowChan%, ByVal HighChan%, ByVal CBCount&, CBRate&, ByVal Gain%, ByVal MemHandle%, ByVal Options%) As Integer
Declare Function cbAPretrig Lib "cbw.dll" (ByVal BoardNum%, ByVal LowChan%, ByVal HighChan%, PretrigCount&, CBCount&, CBRate&, ByVal Gain%, ByVal MemHandle%, ByVal Options%) As Integer
Declare Function cbATrig Lib "cbw.dll" (ByVal BoardNum%, ByVal Chan%, ByVal TrigType%, ByVal TrigValue%, ByVal Gain%, DataValue%) As Integer
Declare Function cbC8254Config Lib "cbw.dll" (ByVal BoardNum%, ByVal CounterNum%, ByVal Config%) As Integer
Declare Function cbC8536Config Lib "cbw.dll" (ByVal BoardNum%, ByVal CounterNum%, ByVal OutputControl%, ByVal RecycleMode%, ByVal Retrigger%) As Integer
Declare Function cbC9513Config Lib "cbw.dll" (ByVal BoardNum%, ByVal CounterNum%, ByVal GateControl%, ByVal CounterEdge%, ByVal CountSource%, ByVal SpecialGate%, ByVal Reload%, ByVal RecycleMode%, ByVal BCDMode%, ByVal CountDirec%, ByVal OutputCtrl%) As Integer
Declare Function cbC8536Init Lib "cbw.dll" (ByVal BoardNum%, ByVal ChipNum%, ByVal Ctr1Output%) As Integer
Declare Function cbC9513Init Lib "cbw.dll" (ByVal BoardNum%, ByVal ChipNum%, ByVal FOutDivider%, ByVal FOutSource%, ByVal Compare1%, ByVal Compare2%, ByVal TimeOfDay%) As Integer
Declare Function cbCStoreOnInt Lib "cbw.dll" (ByVal BoardNum%, ByVal IntCount%, CntrControl%, ByVal DataBuffer%) As Integer
Declare Function cbCFreqIn Lib "cbw.dll" (ByVal BoardNum%, ByVal SigSource%, ByVal GateInterval%, CBCount%, Freq&) As Integer
Declare Function cbCIn Lib "cbw.dll" (ByVal BoardNum%, ByVal CounterNum%, CBCount As Any) As Integer
Declare Function cbCLoad Lib "cbw.dll" (ByVal BoardNum%, ByVal RegNum%, ByVal LoadValue%) As Integer
Declare Function cbDBitIn Lib "cbw.dll" (ByVal BoardNum%, ByVal PortType%, ByVal BitNum%, BitValue%) As Integer
Declare Function cbDBitOut Lib "cbw.dll" (ByVal BoardNum%, ByVal PortType%, ByVal BitNum%, ByVal BitValue%) As Integer
Declare Function cbDConfigPort Lib "cbw.dll" (ByVal BoardNum%, ByVal PortNum%, ByVal Direction%) As Integer
Declare Function cbDeclareRevision Lib "cbw.dll" (RevNum!) As Integer
Declare Function cbDIn Lib "cbw.dll" (ByVal BoardNum%, ByVal PortNum%, DataValue%) As Integer
Declare Function cbDInScan Lib "cbw.dll" (ByVal BoardNum%, ByVal PortNum%, ByVal CBCount&, CBRate&, ByVal MemHandle%, ByVal Options%) As Integer
Declare Function cbDOut Lib "cbw.dll" (ByVal BoardNum%, ByVal PortNum%, ByVal DataValue%) As Integer
Declare Function cbDOutScan Lib "cbw.dll" (ByVal BoardNum%, ByVal PortNum%, ByVal CBCount&, CBRate&, ByVal MemHandle%, ByVal Options%) As Integer
Declare Function cbErrHandling Lib "cbw.dll" (ByVal ErrReporting%, ByVal ErrHandling%) As Integer
Declare Function cbFileAInScan Lib "cbw.dll" (ByVal BoardNum%, ByVal LowChan%, ByVal HighChan%, ByVal CBCount&, CBRate&, ByVal Gain%, ByVal FileName$, ByVal Options%) As Integer
Declare Function cbFileGetInfo Lib "cbw.dll" (ByVal FileName$, LowChan%, HighChan%, PretrigCount&, TotalCount&, CBRate&, Gain As Any) As Integer
Declare Function cbFilePretrig Lib "cbw.dll" (ByVal BoardNum%, ByVal LowChan%, ByVal HighChan%, PretrigCount&, CBCount&, CBRate&, ByVal Gain%, ByVal FileName$, ByVal Options%) As Integer
Declare Function cbFileRead Lib "cbw.dll" (ByVal FileName$, ByVal FirstPoint&, NumPoints&, DataBuffer%) As Integer
Declare Function cbGetErrMsg Lib "cbw.dll" (ByVal ErrCode%, ByVal ErrMsg$) As Integer
Declare Function cbGetRevision Lib "cbw.dll" (DLLRevNum!, VXDRevNum!) As Integer
Declare Function cbGetStatus Lib "cbw.dll" (ByVal BoardNum%, Status%, CurCount&, CurIndex&) As Integer
Declare Function cbStopBackground Lib "cbw.dll" (ByVal BoardNum%) As Integer
Declare Function cbMemSetDTMode Lib "cbw.dll" (ByVal BoardNum%, ByVal Mode%) As Integer
Declare Function cbMemReset Lib "cbw.dll" (ByVal BoardNum%) As Integer
Declare Function cbMemRead Lib "cbw.dll" (ByVal BoardNum%, DataBuffer%, ByVal FirstPoint&, ByVal CBCount&) As Integer
Declare Function cbMemWrite Lib "cbw.dll" (ByVal BoardNum%, DataBuffer%, ByVal FirstPoint&, ByVal CBCount&) As Integer
Declare Function cbMemReadPretrig Lib "cbw.dll" (ByVal BoardNum%, DataBuffer%, ByVal FirstPoint&, ByVal CBCount&) As Integer
Declare Function cbRS485 Lib "cbw.dll" (ByVal BoardNum%, ByVal Transmit%, ByVal Receive%) As Integer
Declare Function cbTIn Lib "cbw.dll" (ByVal BoardNum%, ByVal Chan%, ByVal CBScale%, TempValue!, ByVal Options%) As Integer
Declare Function cbTInScan Lib "cbw.dll" (ByVal BoardNum%, ByVal LowChan%, ByVal HighChan%, ByVal CBScale%, DataBuffer!, ByVal Options%) As Integer
Declare Function cbWinBufToArray Lib "cbw.dll" (ByVal MemHandle%, DataBuffer%, ByVal FirstPoint&, ByVal CBCount&) As Integer
Declare Function cbWinArrayToBuf Lib "cbw.dll" (DataBuffer%, ByVal MemHandle%, ByVal FirstPoint&, ByVal CBCount&) As Integer
Declare Function cbWinBufAlloc Lib "cbw.dll" (ByVal NumPoints&) As Integer
Declare Function cbWinBufFree Lib "cbw.dll" (ByVal MemHandle%) As Integer
Declare Function cbInByte Lib "cbw.dll" (ByVal BoardNum%, ByVal PortNum%) As Integer
Declare Function cbOutByte Lib "cbw.dll" (ByVal BoardNum%, ByVal PortNum%, ByVal PortVal%) As Integer
Declare Function cbInWord Lib "cbw.dll" (ByVal BoardNum%, ByVal PortNum%) As Integer
Declare Function cbOutWord Lib "cbw.dll" (ByVal BoardNum%, ByVal PortNum%, ByVal PortVal%) As Integer
Declare Function cbGetConfig Lib "cbw.dll" (ByVal InfoType%, ByVal BoardNum%, ByVal DevNum%, ByVal ConfigItem%, ConfigVal%) As Integer
Declare Function cbSetConfig Lib "cbw.dll" (ByVal InfoType%, ByVal BoardNum%, ByVal DevNum%, ByVal ConfigItem%, ByVal ConfigVal%) As Integer
Declare Function cbToEngUnits Lib "cbw.dll" (ByVal BoardNum%, ByVal Range%, ByVal DataVal%, EngUnits!) As Integer
Declare Function cbFromEngUnits Lib "cbw.dll" (ByVal BoardNum%, ByVal Range%, ByVal EngUnits!, DataVal%) As Integer
Declare Function cbGetBoardName Lib "cbw.dll" (ByVal BoardNum%, ByVal BoardName$) As Integer
Declare Function cbSetTrigger Lib "cbw.dll" (ByVal BoardNum%, ByVal TrigType%, ByVal LowThreshold%, ByVal HighThreshold%) As Integer
#End If

