VERSION 5.00
Begin VB.Form frmMain 
   AutoRedraw      =   -1  'True
   Caption         =   "Bar & Wave V6.1.2"
   ClientHeight    =   8670
   ClientLeft      =   1320
   ClientTop       =   1380
   ClientWidth     =   15240
   Icon            =   "Form1.frx":0000
   KeyPreview      =   -1  'True
   LinkTopic       =   "Form1"
   PaletteMode     =   1  'UseZOrder
   ScaleHeight     =   578
   ScaleMode       =   3  'Pixel
   ScaleWidth      =   1016
   Begin VB.Frame Frame9 
      Caption         =   "Stim Position"
      Height          =   1095
      Left            =   8640
      TabIndex        =   88
      Top             =   4920
      Width           =   6495
      Begin VB.Label lblStimPosY 
         Caption         =   "N/A"
         Height          =   255
         Left            =   2880
         TabIndex        =   90
         Top             =   360
         Width           =   3255
      End
      Begin VB.Label lblStimPosX 
         Caption         =   "N/A"
         Height          =   255
         Left            =   120
         TabIndex        =   89
         Top             =   360
         Width           =   2535
      End
   End
   Begin VB.Frame Frame7 
      Caption         =   "Avg Counts"
      Height          =   4335
      Left            =   8640
      TabIndex        =   84
      Top             =   480
      Width           =   6495
      Begin VB.CommandButton cmdSaveCounts 
         Caption         =   "Save As..."
         Height          =   495
         Left            =   240
         TabIndex        =   86
         Top             =   3720
         Width           =   2415
      End
      Begin VB.PictureBox pctCounts 
         AutoRedraw      =   -1  'True
         BackColor       =   &H80000007&
         Height          =   3255
         Left            =   120
         ScaleHeight     =   3195
         ScaleWidth      =   6195
         TabIndex        =   85
         Top             =   240
         Width           =   6255
      End
   End
   Begin VB.CheckBox chkTrigger 
      Caption         =   "Start On Trigger"
      Height          =   255
      Left            =   240
      TabIndex        =   83
      Top             =   6360
      Width           =   1935
   End
   Begin VB.CheckBox chkUseLUT 
      Caption         =   "Use color LUT"
      Height          =   255
      Left            =   360
      TabIndex        =   53
      Top             =   4560
      Width           =   1575
   End
   Begin VB.Timer Timer1 
      Enabled         =   0   'False
      Interval        =   1000
      Left            =   0
      Top             =   840
   End
   Begin VB.TextBox txtFramesOfStim 
      Height          =   315
      Left            =   1680
      TabIndex        =   43
      Text            =   "10"
      Top             =   7320
      Width           =   555
   End
   Begin VB.TextBox txtFramesBeforeStim 
      Height          =   315
      Left            =   1680
      TabIndex        =   42
      Text            =   "10"
      Top             =   6840
      Width           =   555
   End
   Begin VB.CommandButton cmdLoadLUT 
      Caption         =   "Load Color LUT"
      Height          =   495
      Left            =   480
      TabIndex        =   38
      TabStop         =   0   'False
      Top             =   4920
      Width           =   1215
   End
   Begin VB.TextBox txtDistance 
      Height          =   315
      Left            =   600
      TabIndex        =   34
      TabStop         =   0   'False
      Text            =   "114"
      Top             =   4080
      Width           =   975
   End
   Begin VB.TextBox txtHSize 
      Height          =   315
      Left            =   600
      TabIndex        =   32
      TabStop         =   0   'False
      Text            =   "6"
      Top             =   3360
      Width           =   975
   End
   Begin VB.CommandButton cmdGetFrameRate 
      Caption         =   "Get Rate"
      Height          =   495
      Left            =   600
      TabIndex        =   30
      TabStop         =   0   'False
      Top             =   2400
      Width           =   1215
   End
   Begin VB.TextBox txtFrameRate 
      Height          =   315
      Left            =   600
      TabIndex        =   29
      TabStop         =   0   'False
      Text            =   "60"
      Top             =   1560
      Width           =   975
   End
   Begin VB.CommandButton cmdQuit 
      Caption         =   "Quit"
      Height          =   495
      Left            =   240
      TabIndex        =   15
      TabStop         =   0   'False
      Top             =   7920
      Width           =   1575
   End
   Begin VB.Frame Frame1 
      Caption         =   "Square Wave"
      Height          =   8115
      Left            =   2400
      TabIndex        =   5
      Top             =   480
      Width           =   6135
      Begin VB.Frame Frame8 
         Caption         =   "Mask"
         Height          =   2295
         Left            =   4200
         TabIndex        =   87
         Top             =   960
         Width           =   1815
         Begin VB.CheckBox chkLockPhaseToMask 
            Caption         =   "Lock Phase to Mask"
            Height          =   495
            Left            =   120
            TabIndex        =   2
            TabStop         =   0   'False
            Top             =   1680
            Width           =   1455
         End
         Begin VB.CheckBox chkUseStimOri 
            Caption         =   "Use Stim Orientation"
            Height          =   375
            Left            =   120
            TabIndex        =   1
            TabStop         =   0   'False
            Top             =   960
            Width           =   1575
         End
         Begin VB.CheckBox chkSquareMask 
            Caption         =   "Use square mask"
            Height          =   255
            Left            =   120
            TabIndex        =   0
            TabStop         =   0   'False
            Top             =   360
            Width           =   1575
         End
      End
      Begin VB.Frame Frame6 
         Caption         =   "Wave Type"
         Height          =   1215
         Left            =   4200
         TabIndex        =   80
         Top             =   3360
         Width           =   1815
         Begin VB.OptionButton optWaveType 
            Caption         =   "Sin"
            Height          =   255
            Index           =   1
            Left            =   120
            TabIndex        =   82
            Top             =   720
            Width           =   1455
         End
         Begin VB.OptionButton optWaveType 
            Caption         =   "Square"
            Height          =   255
            Index           =   0
            Left            =   120
            TabIndex        =   81
            Top             =   360
            Value           =   -1  'True
            Width           =   1455
         End
      End
      Begin VB.Frame Frame5 
         Caption         =   "Info Counters"
         Height          =   1695
         Left            =   4800
         TabIndex        =   76
         Top             =   4800
         Width           =   1215
         Begin VB.Label Label33 
            Caption         =   "N/A"
            Height          =   255
            Left            =   120
            TabIndex        =   79
            Top             =   1320
            Width           =   855
         End
         Begin VB.Label Label32 
            Caption         =   "N/A"
            Height          =   255
            Left            =   120
            TabIndex        =   78
            Top             =   840
            Width           =   855
         End
         Begin VB.Label Label31 
            Caption         =   "N/A"
            Height          =   255
            Left            =   120
            TabIndex        =   77
            Top             =   360
            Width           =   855
         End
      End
      Begin VB.TextBox txtFlashPeriodInFrames 
         Enabled         =   0   'False
         Height          =   285
         Left            =   120
         TabIndex        =   74
         Text            =   "0"
         Top             =   6360
         Width           =   1095
      End
      Begin VB.TextBox txtRotationPeriodInFrames 
         Enabled         =   0   'False
         Height          =   285
         Left            =   120
         TabIndex        =   72
         Text            =   "60"
         Top             =   5160
         Width           =   1095
      End
      Begin VB.CheckBox chkDoubleStim 
         Caption         =   "Double Stim"
         Height          =   255
         Left            =   2640
         TabIndex        =   68
         Top             =   4440
         Width           =   1695
      End
      Begin VB.Frame Frame4 
         Caption         =   "Flash"
         Height          =   1335
         Left            =   2640
         TabIndex        =   64
         Top             =   960
         Width           =   1455
         Begin VB.OptionButton optFlashType 
            Caption         =   "Inverse"
            Height          =   255
            Index           =   2
            Left            =   120
            TabIndex        =   67
            Top             =   960
            Value           =   -1  'True
            Width           =   975
         End
         Begin VB.OptionButton optFlashType 
            Caption         =   "Max - Mean"
            Height          =   255
            Index           =   1
            Left            =   120
            TabIndex        =   66
            Top             =   600
            Width           =   1215
         End
         Begin VB.OptionButton optFlashType 
            Caption         =   "Max - Min"
            Height          =   255
            Index           =   0
            Left            =   120
            TabIndex        =   65
            Top             =   240
            Width           =   1095
         End
      End
      Begin VB.TextBox txtFlashPeriodInHz 
         Height          =   285
         Left            =   120
         TabIndex        =   62
         Text            =   "0"
         Top             =   6000
         Width           =   1095
      End
      Begin VB.Frame Frame3 
         Caption         =   "Frame Timing"
         Height          =   1695
         Left            =   2640
         TabIndex        =   56
         Top             =   4800
         Width           =   2055
         Begin VB.Label lblMissedFrames 
            Caption         =   "0"
            Height          =   255
            Left            =   1320
            TabIndex        =   71
            Top             =   1080
            Width           =   615
         End
         Begin VB.Label Label7 
            Caption         =   "Missed Frames:"
            Height          =   255
            Left            =   120
            TabIndex        =   70
            Top             =   1080
            Width           =   1215
         End
         Begin VB.Label lblMaxdT 
            Caption         =   "0"
            Height          =   255
            Left            =   840
            TabIndex        =   60
            Top             =   720
            Width           =   1095
         End
         Begin VB.Label lblTotaldT 
            Caption         =   "0"
            Height          =   255
            Left            =   840
            TabIndex        =   59
            Top             =   360
            Width           =   1095
         End
         Begin VB.Label Label28 
            Caption         =   "Max dT:"
            Height          =   255
            Left            =   120
            TabIndex        =   58
            Top             =   720
            Width           =   735
         End
         Begin VB.Label Label27 
            Caption         =   "Total dT:"
            Height          =   255
            Left            =   120
            TabIndex        =   57
            Top             =   360
            Width           =   735
         End
      End
      Begin VB.Frame Frame2 
         Caption         =   "Drawing time"
         Height          =   615
         Left            =   2640
         TabIndex        =   54
         Top             =   6600
         Width           =   2055
         Begin VB.Label lblDrT 
            Caption         =   "0"
            Height          =   255
            Left            =   120
            TabIndex        =   55
            Top             =   240
            Width           =   1455
         End
      End
      Begin VB.TextBox txtBKGColor 
         Height          =   285
         Left            =   2640
         TabIndex        =   51
         Text            =   "0.5"
         Top             =   3960
         Width           =   615
      End
      Begin VB.TextBox txtBarW 
         Height          =   285
         Left            =   120
         TabIndex        =   47
         Text            =   "1"
         Top             =   1560
         Width           =   1095
      End
      Begin VB.TextBox txtRotationPeriodInHz 
         Height          =   285
         Left            =   120
         TabIndex        =   44
         Text            =   "1"
         Top             =   4800
         Width           =   1095
      End
      Begin VB.CommandButton cmdLoadParamW 
         Caption         =   "Load Parameters"
         Height          =   555
         Left            =   2640
         TabIndex        =   39
         TabStop         =   0   'False
         Top             =   240
         Width           =   3255
      End
      Begin VB.CommandButton cmdStop 
         Caption         =   "Stop"
         Height          =   495
         Left            =   3240
         TabIndex        =   4
         TabStop         =   0   'False
         Top             =   7440
         Width           =   2655
      End
      Begin VB.CommandButton cmdMotion 
         Caption         =   "Start"
         Height          =   495
         Left            =   120
         TabIndex        =   3
         TabStop         =   0   'False
         Top             =   7440
         Width           =   2655
      End
      Begin VB.TextBox txtTempPeriodInFrames 
         Enabled         =   0   'False
         Height          =   285
         Left            =   120
         TabIndex        =   21
         TabStop         =   0   'False
         Text            =   "170"
         Top             =   3360
         Width           =   1095
      End
      Begin VB.TextBox txtSpPeriod 
         Height          =   285
         Left            =   120
         TabIndex        =   18
         TabStop         =   0   'False
         Text            =   "5"
         Top             =   960
         Width           =   1095
      End
      Begin VB.TextBox txtPhase 
         Height          =   285
         Left            =   120
         TabIndex        =   17
         TabStop         =   0   'False
         Text            =   "0"
         Top             =   2400
         Width           =   1095
      End
      Begin VB.TextBox txtMean 
         Height          =   285
         Left            =   2640
         TabIndex        =   14
         TabStop         =   0   'False
         Text            =   "0.5"
         Top             =   3360
         Width           =   495
      End
      Begin VB.TextBox txtContrast 
         Height          =   285
         Left            =   2640
         TabIndex        =   12
         TabStop         =   0   'False
         Text            =   "1"
         Top             =   2760
         Width           =   495
      End
      Begin VB.TextBox txtTempPeriodInHz 
         Height          =   285
         Left            =   120
         TabIndex        =   10
         TabStop         =   0   'False
         Text            =   "1"
         Top             =   3000
         Width           =   1095
      End
      Begin VB.TextBox txtOrientation 
         Height          =   285
         Left            =   120
         TabIndex        =   8
         TabStop         =   0   'False
         Text            =   "0"
         Top             =   4200
         Width           =   1095
      End
      Begin VB.Label Label35 
         Caption         =   "frames"
         Height          =   255
         Left            =   1320
         TabIndex        =   75
         Top             =   6360
         Width           =   495
      End
      Begin VB.Label Label34 
         Caption         =   "frames"
         Height          =   255
         Left            =   1320
         TabIndex        =   73
         Top             =   5160
         Width           =   495
      End
      Begin VB.Label Label30 
         Caption         =   "Hz"
         Height          =   255
         Left            =   1320
         TabIndex        =   63
         Top             =   6000
         Width           =   375
      End
      Begin VB.Label Label29 
         Caption         =   "Flash Period"
         Height          =   255
         Left            =   120
         TabIndex        =   61
         Top             =   5760
         Width           =   975
      End
      Begin VB.Label Label26 
         Caption         =   "0-1"
         Height          =   255
         Left            =   3360
         TabIndex        =   52
         Top             =   3960
         Width           =   375
      End
      Begin VB.Label Label25 
         Caption         =   "BKG Color"
         Height          =   255
         Left            =   2640
         TabIndex        =   50
         Top             =   3720
         Width           =   735
      End
      Begin VB.Label Label24 
         Caption         =   "Bar Width"
         Height          =   255
         Left            =   120
         TabIndex        =   49
         Top             =   1320
         Width           =   975
      End
      Begin VB.Label Label21 
         Caption         =   "degree"
         Height          =   255
         Left            =   1320
         TabIndex        =   48
         Top             =   1560
         Width           =   615
      End
      Begin VB.Label Label20 
         Caption         =   "Hz"
         Height          =   255
         Left            =   1320
         TabIndex        =   46
         Top             =   4800
         Width           =   375
      End
      Begin VB.Label Label19 
         Caption         =   "Rotational Period"
         Height          =   255
         Left            =   120
         TabIndex        =   45
         Top             =   4560
         Width           =   1215
      End
      Begin VB.Label Label14 
         Caption         =   "degree"
         Height          =   195
         Left            =   1320
         TabIndex        =   26
         Top             =   2400
         Width           =   495
      End
      Begin VB.Label Label13 
         Caption         =   "degree"
         Height          =   195
         Left            =   1320
         TabIndex        =   25
         Top             =   4200
         Width           =   495
      End
      Begin VB.Label Label12 
         Caption         =   "0-1"
         Height          =   195
         Left            =   3360
         TabIndex        =   24
         Top             =   3360
         Width           =   315
      End
      Begin VB.Label Label11 
         Caption         =   "0-1"
         Height          =   255
         Left            =   3360
         TabIndex        =   23
         Top             =   2760
         Width           =   315
      End
      Begin VB.Label Label10 
         Caption         =   "frames"
         Height          =   255
         Left            =   1320
         TabIndex        =   22
         Top             =   3360
         Width           =   555
      End
      Begin VB.Label Label9 
         Caption         =   "Hz"
         Height          =   195
         Left            =   1320
         TabIndex        =   20
         Top             =   3000
         Width           =   435
      End
      Begin VB.Label Label8 
         Caption         =   "degree"
         Height          =   195
         Left            =   1320
         TabIndex        =   19
         Top             =   960
         Width           =   555
      End
      Begin VB.Label Label5 
         Caption         =   "Phase:"
         Height          =   195
         Left            =   120
         TabIndex        =   16
         Top             =   2160
         Width           =   735
      End
      Begin VB.Label lblMean 
         Caption         =   "Mean:"
         Height          =   255
         Left            =   2640
         TabIndex        =   13
         Top             =   3120
         Width           =   495
      End
      Begin VB.Label Label4 
         Caption         =   "Contrast:"
         Height          =   195
         Left            =   2640
         TabIndex        =   11
         Top             =   2520
         Width           =   615
      End
      Begin VB.Label Label3 
         Caption         =   "Temporal period:"
         Height          =   255
         Left            =   120
         TabIndex        =   9
         Top             =   2760
         Width           =   1215
      End
      Begin VB.Label Label2 
         Caption         =   "Spacial period:"
         Height          =   255
         Left            =   120
         TabIndex        =   7
         Top             =   720
         Width           =   1155
      End
      Begin VB.Label Label1 
         Caption         =   "Orientation:"
         Height          =   255
         Left            =   120
         TabIndex        =   6
         Top             =   3960
         Width           =   915
      End
   End
   Begin VB.Label lbl2 
      Height          =   735
      Left            =   8760
      TabIndex        =   92
      Top             =   7440
      Width           =   5535
   End
   Begin VB.Label lbl1 
      Height          =   495
      Left            =   8880
      TabIndex        =   91
      Top             =   6480
      Width           =   5295
   End
   Begin VB.Label lblFrames 
      Caption         =   "0"
      Height          =   255
      Left            =   600
      TabIndex        =   69
      Top             =   2040
      Width           =   855
   End
   Begin VB.Label Label23 
      Caption         =   "Frames of Stim:"
      Height          =   255
      Left            =   120
      TabIndex        =   41
      Top             =   7320
      Width           =   1395
   End
   Begin VB.Label Label22 
      Caption         =   "Frames before Stim:"
      Height          =   255
      Left            =   120
      TabIndex        =   40
      Top             =   6840
      Width           =   1395
   End
   Begin VB.Label lblParamFileName 
      Caption         =   "N/A"
      Height          =   255
      Left            =   4200
      TabIndex        =   37
      Top             =   120
      Width           =   10095
   End
   Begin VB.Label Label18 
      Caption         =   "Parameter File:"
      Height          =   255
      Left            =   2400
      TabIndex        =   36
      Top             =   120
      Width           =   1095
   End
   Begin VB.Label lblResolution 
      Caption         =   "NA"
      Height          =   255
      Left            =   600
      TabIndex        =   35
      Top             =   720
      Width           =   975
   End
   Begin VB.Label Label17 
      Caption         =   "Distance"
      Height          =   255
      Left            =   600
      TabIndex        =   33
      Top             =   3720
      Width           =   1035
   End
   Begin VB.Label Label16 
      Caption         =   "Horizontal Size"
      Height          =   315
      Left            =   600
      TabIndex        =   31
      Top             =   3000
      Width           =   1155
   End
   Begin VB.Label Label6 
      Caption         =   "Frame Rate:"
      Height          =   255
      Left            =   600
      TabIndex        =   28
      Top             =   1200
      Width           =   1035
   End
   Begin VB.Label Label15 
      Caption         =   "Resolution:"
      Height          =   255
      Left            =   600
      TabIndex        =   27
      Top             =   240
      Width           =   975
   End
End
Attribute VB_Name = "frmMain"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Option Explicit


Private Sub chkLockPhaseToMask_Click()
frmMain.SetFocus
End Sub

Private Sub chkSquareMask_Click()
frmMain.SetFocus
End Sub

Private Sub chkUseStimOri_Click()
frmMain.SetFocus
End Sub

Private Sub cmdGetFrameRate_Click()
Dim j As Integer
Dim Trash As Long
Dim StartT As Currency
Dim StopT As Currency
Dim Mfr As Currency
Dim N As Integer


Trash = QueryPerformanceFrequency(Mfr)

frmStimulus.ZOrder
bStop = False

FillTwoSurf RGB(0, 0, 0)
ddsBack.SetForeColor RGB(255, 0, 0)
ddsBack.setDrawWidth 50
ddsBack.DrawCircle ResolutionX / 2 + 20, ResolutionY / 2, 25

Flip
WaitFrameStart

ddsBack.SetForeColor RGB(0, 255, 0)
ddsBack.setDrawWidth 50
ddsBack.DrawCircle ResolutionX / 2 - 20, ResolutionY / 2, 25

WaitFrameStart

WaitFrameStart
Trash = QueryPerformanceCounter(StartT)
For j = 1 To 2048
    DoEvents
    Flip
    If bStop Then
        Exit For
    End If
    frmMain.lblFrames.Caption = CStr(j)
    WaitFrameStart
Next j
Trash = QueryPerformanceCounter(StopT)

FrameRate = ((Mfr * 2048) / (StopT - StartT))
txtFrameRate.Text = Format(FrameRate, "###0.000000")
RecalculateTempPeriod
ShowParametersW WaveLoopSequence(0)

End Sub

Private Sub cmdLoadLUT_Click()
Dim ssFile As String
On Error GoTo ErrorExit
    frmCDHolder.CommonDialog1.FILTER = "LUT file (*.lut)|*.lut|All files(*.*)|*.*"
    frmCDHolder.CommonDialog1.ShowOpen
    ssFile = frmCDHolder.CommonDialog1.FileName
    LUTFname = ssFile
    ReadLUT (ssFile)

ErrorExit:

End Sub


Private Sub cmdLoadParamW_Click()
bStop = True
On Error GoTo Err:
    frmCDHolder.CommonDialog1.FILTER = "Parameter file (*.x)|*.x|All files(*.*)|*.*"
    frmCDHolder.CommonDialog1.ShowOpen
    ParamFileName = frmCDHolder.CommonDialog1.FileName
    frmMain.cmdMotion.Enabled = False
    ReadParamW (ParamFileName)
    frmMain.lblParamFileName = ParamFileName
 '   MsgBox "Parameter File Loaded"
Err:
    frmMain.cmdMotion.Enabled = True
End Sub



Private Sub cmdQuit_Click()
    DestroySurface
    RestoreDispMode
Unload frmStimulus
Set objDX = Nothing
Set objDD = Nothing
End
End Sub

Private Sub cmdMotion_Click()
Dim ii As Long
Dim Trash As Long
Dim TrialIndex As Long
Dim Angle As Double
Dim dL As Double

Dim i As Integer
Dim j As Long
Dim ExitLoop As Boolean

Dim BarIsOn As Boolean 'to create flashing bar

Dim CurrentCounterReading As Long
Dim NewCounterReading As Long
Dim BaseCounterReading As Long

Dim CurrentCounterReadingAvg As Long
Dim NewCounterReadingAvg As Long
Dim avgIndex As Integer
Dim avgstartTime As Currency
Dim avgstopTime As Currency

Dim Frame As Long
Dim scFrames As Long

Dim loopStr As String
Dim cumLoopSize As Long

Dim getPuls As Boolean
Dim FlashPhase As Double


Dim startTime5 As Currency
Dim stopTime5 As Currency


frmMain.SetFocus

'Exit Sub

scFrames = 0
frmStimulus.WindowState = 2
RecalculateTempPeriod
WriteIni (App.Path & "\" & IniFileName)
WriteLog ParamFileName
Trash = QueryPerformanceFrequency(tikPerSec)
bStop = False
CalculateCompPar

If loopsizeToAVG > 0 Then
    For i = 0 To loopsizeToAVG * 2 - 1
        TotalCounts(i) = 0
        AvgCounts(i) = 0
        TotalTime(i) = 0
    Next i
    ShowAvg
    SetCounter CounterNforAvg, 0
End If
SetCounter CounterN, 0
SetCounter CounterNTrigger, 0



'set all TTL outs to initial state
For i = 0 To MarkNumber - 1
    Call SetTTLMark(MarkArray(i).InitialState, MarkArray(i).MarkCh)
Next i
WaitMs 20

CurrentLUTStart = WaveLoopSequence(0).LUTStart
CurrentLUTEnd = WaveLoopSequence(0).LUTEnd

GCon.SetGammaRamp DDSGR_DEFAULT, InitialRamp

FillTwoSurf (MakeColor(WaveLoopSequence(0).BKGColor))


'wait for trigger
If frmMain.chkTrigger.Value = 1 Then
    CurrentCounterReading = ReadCounter(CounterNTrigger)
    If CurrentCounterReading = 100 Then
        SetCounter CounterNTrigger, 101
    Else
        SetCounter CounterNTrigger, 100
    End If
        CurrentCounterReading = ReadCounter(CounterNTrigger)
        ExitLoop = False
        Do While ExitLoop = False
            DoEvents
            If bStop Then
                Exit Do
            End If
            NewCounterReading = ReadCounter(CounterNTrigger)
            frmMain.Caption = CStr(NewCounterReading)
            If NewCounterReading <> CurrentCounterReading Then
                ExitLoop = True
            End If
        Loop
End If
ResetCounter CounterN
BaseCounterReading = 0

'set start TTL out
For i = 0 To MarkNumber - 1
    If InStr(MarkArray(i).MarkTiming, "Start") > 0 Then
        Call SetTTLMark(MarkArray(i).InitialState + 1, MarkArray(i).MarkCh)
    End If
Next i

Dim TrIndStep As Integer
TrIndStep = 1 + frmMain.chkDoubleStim.Value

For TrialIndex = 0 To LoopSequenceSize - 1 Step TrIndStep
'    CurrentCounterReading = ReadCounter(CounterN)
 '   ResetCounter CounterN
  '  getPuls = False



If loopsizeToAVG > 0 Then
    avgIndex = TrialIndex Mod loopsizeToAVG
    CurrentCounterReadingAvg = ReadCounter(CounterNforAvg)
    If CurrentCounterReadingAvg > 10000 Then
        Call SetCounter(CounterNforAvg, 0)
    End If
End If
    Trash = QueryPerformanceCounter(avgstartTime)
    
    GCon.SetGammaRamp DDSGR_DEFAULT, InitialRamp

    CurrentLUTStart = WaveLoopSequence(TrialIndex).LUTStart
    CurrentLUTEnd = WaveLoopSequence(TrialIndex).LUTEnd

    'set loop marks on start of blank
    For i = 0 To MarkNumber - 1
        If MarkArray(i).Position = 0 Then
            'Loop0 mark - each new parameter set
            'Loop1... - mark at start of corresponding loop
            cumLoopSize = 1
            For ii = 0 To LoopsNumber - 1
                loopStr = "Loop" & CStr(ii)
                If InStr(MarkArray(i).MarkTiming, loopStr) > 0 Then
                    If TrialIndex Mod (cumLoopSize * MarkArray(i).MarkPeriod) = 0 Then
                        If MarkArray(i).MarkType = 0 Then
                            Call SetTTLMark(MarkArray(i).InitialState + 1, MarkArray(i).MarkCh)
                        End If
                        If MarkArray(i).MarkType = 1 Then
                            Call SetTTLMark(((TrialIndex / (cumLoopSize * MarkArray(i).MarkPeriod)) Mod 2) + MarkArray(i).InitialState + 1, MarkArray(i).MarkCh)
                        End If
                    End If
                End If
                cumLoopSize = cumLoopSize * LoopSizes(ii)
            Next ii
        End If
    Next i

    Call ShowParametersW(WaveLoopSequence(TrialIndex))
    
    
  '''a  CurrentCounterReading = ReadCounter(CounterN)
  '''a  ResetCounter CounterN
    getPuls = False
    If FramesBeforeStim > 0 Then 'rev 621
        FillTwoSurf (MakeColor(WaveLoopSequence(TrialIndex).BKGColor))
    End If
    

    If frmMain.optWaveType(1).Value = True Then
        'prepare sin wave on BKGImage surface
        Call PreparePictureWaveForBlt(TrialIndex, dL, Angle, 1)
    End If
        
    'waiting N pulses to start
    If FramesBeforeStim > 0 Then
 'moved up  rev621     FillTwoSurf (MakeColor(WaveLoopSequence(TrialIndex).BKGColor))
        'reset all TTL outs to initial state if puls type=0 (short puls)
        For i = 0 To MarkNumber - 1
            If MarkArray(i).MarkType = 0 Then
                Call SetTTLMark(MarkArray(i).InitialState, MarkArray(i).MarkCh)
            End If
        Next i

        ExitLoop = False
        Trash = QueryPerformanceCounter(startTime5)
        Do While ExitLoop = False
            DoEvents
            If bStop Then
                Exit Do
            End If
            NewCounterReading = ReadCounter(CounterN)
            frmMain.Caption = CStr(NewCounterReading - BaseCounterReading)
'''a            If NewCounterReading <> CurrentCounterReading Then
            
'''a                 Trash = QueryPerformanceCounter(stopTime5)
'''a                 getPuls = True
'''a                 frmMain.lbl1.Caption = Format((stopTime5 - startTime5) * 1000# / tikPerSec, "###0.00000")
'''a                 startTime5 = stopTime5
                
'''a             End If
'''a             CurrentCounterReading = NewCounterReading 'SY
            
'''a             If getPuls = True Then
'''a                 If NewCounterReading > FramesBeforeStim - 1 Then
'''a                     ExitLoop = True
'''a                 End If
'''a             End If
            If NewCounterReading > BaseCounterReading + FramesBeforeStim - 1 + (FramesBeforeStim + FramesOfStim) * TrialIndex / TrIndStep Then
                ExitLoop = True
            End If
        Loop
    End If

    TotaldT = 0
    MaxdT = 0
    frmMain.lblMaxdT.Caption = Format(MaxdT, "###0.00000")
    frmMain.lblTotaldT.Caption = Format(TotaldT, "###0.00000")
    frmMain.lblMissedFrames.Caption = CStr(scFrames)
    
    Trash = QueryPerformanceCounter(startTime)

    Frame = 0
    ExitLoop = False
    
    If loopsizeToAVG > 0 Then
        NewCounterReadingAvg = ReadCounter(CounterNforAvg)
        If CurrentCounterReadingAvg > NewCounterReadingAvg Then
            TotalCounts(avgIndex * 2) = TotalCounts(avgIndex * 2) + NewCounterReadingAvg
        Else
            TotalCounts(avgIndex * 2) = TotalCounts(avgIndex * 2) + NewCounterReadingAvg - CurrentCounterReadingAvg
        End If
        Trash = QueryPerformanceCounter(avgstopTime)
        TotalTime(avgIndex * 2) = TotalTime(avgIndex * 2) + CDbl((avgstopTime - avgstartTime) / tikPerSec)
        AvgCounts(avgIndex * 2) = TotalCounts(avgIndex * 2) / TotalTime(avgIndex * 2)
        CurrentCounterReading = NewCounterReading
        avgstartTime = avgstopTime
        ShowAvg
    End If
    
    'set loop marks on start of stimuli
    For i = 0 To MarkNumber - 1
        If MarkArray(i).Position = 1 Then
            'Loop0 mark - each new parameter set
            'Loop1... - mark at start of corresponding loop
            cumLoopSize = 1
            For ii = 0 To LoopsNumber - 1
                loopStr = "Loop" & CStr(ii)
                If InStr(MarkArray(i).MarkTiming, loopStr) > 0 Then
                    If TrialIndex Mod (cumLoopSize * MarkArray(i).MarkPeriod) = 0 Then
                        If MarkArray(i).MarkType = 0 Then
                            Call SetTTLMark(MarkArray(i).InitialState + 1, MarkArray(i).MarkCh)
                        End If
                        If MarkArray(i).MarkType = 1 Then
                            Call SetTTLMark(((TrialIndex / (cumLoopSize * MarkArray(i).MarkPeriod)) Mod 2) + MarkArray(i).InitialState + 1, MarkArray(i).MarkCh)
                        End If
                    End If
                End If
                cumLoopSize = cumLoopSize * LoopSizes(ii)
            Next ii
        End If
    Next i

    
    
    Do While ExitLoop = False
        DoEvents

        
        If FramesOfStim > 0 Then
            NewCounterReading = ReadCounter(CounterN)
            frmMain.Caption = CStr(NewCounterReading - BaseCounterReading)
'''a            If NewCounterReading <> CurrentCounterReading Then
'''a                getPuls = True
'''a            End If

'''a            If getPuls = True Then
'''a                If NewCounterReading > FramesBeforeStim + FramesOfStim - 1 Then
'''a                    ExitLoop = True
'                    frmMain.lbl1.Caption = CStr(NewCounterReading)
'''a                End If
'''a            End If
            If NewCounterReading > BaseCounterReading - 1 + (FramesBeforeStim + FramesOfStim) * (1 + TrialIndex / TrIndStep) Then
                ExitLoop = True
            Else
                If NewCounterReading > 20000 Then
                      '''a  CurrentCounterReading = ReadCounter(CounterN)
                    ResetCounter CounterN
                    BaseCounterReading = BaseCounterReading - NewCounterReading
''''                MsgBox (CStr(ReadCounter(CounterN)))
                End If
                
            End If

            
            
            If ExitLoop = True Then
                Exit Do
            End If
        End If
        
        If bStop = True Then
            Exit Do
        End If
        
        frmMain.Label31.Caption = CStr(Frame)
        frmMain.Label32.Caption = CStr(TrialIndex)
        
        Call CalculateParamByFrame(TrialIndex, Frame, dL, Angle, BarIsOn, FlashPhase)
        '************* start period marks
        For i = 0 To MarkNumber - 1
            If InStr(MarkArray(i).MarkTiming, "DriftPeriod") > 0 Then
                If WaveLoopSequence(TrialIndex).WaveTempPeriodInFrames <> 0 Then
                    If Frame Mod (WaveLoopSequence(TrialIndex).WaveTempPeriodInFrames * MarkArray(i).MarkPeriod) = 0 Then
                        If MarkArray(i).MarkType = 0 Then
                            Call SetTTLMark(MarkArray(i).InitialState + 1, MarkArray(i).MarkCh)
                        End If
                        If MarkArray(i).MarkType = 1 Then
                            Call SetTTLMark(((TrialIndex / (cumLoopSize * MarkArray(i).MarkPeriod)) Mod 2) + MarkArray(i).InitialState + 1, MarkArray(i).MarkCh)
                        End If
                    End If
                End If
            End If
            If InStr(MarkArray(i).MarkTiming, "RotationPeriod") > 0 Then
                If WaveLoopSequence(TrialIndex).RotationPeriodInFrames <> 0 Then
                    If Frame Mod (WaveLoopSequence(TrialIndex).RotationPeriodInFrames * MarkArray(i).MarkPeriod) = 0 Then
                        If MarkArray(i).MarkType = 0 Then
                            Call SetTTLMark(MarkArray(i).InitialState + 1, MarkArray(i).MarkCh)
                        End If
                        If MarkArray(i).MarkType = 1 Then
                            Call SetTTLMark(((TrialIndex / (cumLoopSize * MarkArray(i).MarkPeriod)) Mod 2) + MarkArray(i).InitialState + 1, MarkArray(i).MarkCh)
                        End If
                    End If
                End If
            End If
            If InStr(MarkArray(i).MarkTiming, "FlashPeriod") > 0 Then
                If WaveLoopSequence(TrialIndex).FlashPeriodInFrames <> 0 Then
                    If Frame Mod (WaveLoopSequence(TrialIndex).FlashPeriodInFrames * MarkArray(i).MarkPeriod) = 0 Then
                        If MarkArray(i).MarkType = 0 Then
                            Call SetTTLMark(MarkArray(i).InitialState + 1, MarkArray(i).MarkCh)
                        End If
                        If MarkArray(i).MarkType = 1 Then
                            Call SetTTLMark(((TrialIndex / (cumLoopSize * MarkArray(i).MarkPeriod)) Mod 2) + MarkArray(i).InitialState + 1, MarkArray(i).MarkCh)
                        End If
                    End If
                End If
            End If
        Next i
    '***** end pariod marks

        If WaveLoopSequence(TrialIndex).WaveTempPeriodInFrames <> 0 Then
            If Frame > WaveLoopSequence(TrialIndex).WaveTempPeriodInFrames * WaveLoopSequence(TrialIndex).PeriodsToShow - 1 Then
                Exit Do
            End If
        End If

        If WaveLoopSequence(TrialIndex).RotationPeriodInFrames <> 0 Then
            If Frame > WaveLoopSequence(TrialIndex).RotationPeriodInFrames * WaveLoopSequence(TrialIndex).PeriodsToShow - 1 Then
                Exit Do
            End If
        End If

        If WaveLoopSequence(TrialIndex).FlashPeriodInFrames <> 0 Then
            If Frame > WaveLoopSequence(TrialIndex).FlashPeriodInFrames * WaveLoopSequence(TrialIndex).PeriodsToShow - 1 Then
                Exit Do
            End If
        End If

        If BarIsOn Then
            If frmMain.optWaveType(0).Value = True Then
                FillBkgSurf (MakeColor(WaveLoopSequence(TrialIndex).Mean * (1 - WaveLoopSequence(TrialIndex).Contrast)))
                Call CreatePictureWaveDirect(TrialIndex, dL, Angle, 1, False)
            End If
            If frmMain.optWaveType(1).Value = True Then
                Call CreatePictureWaveByBlt(TrialIndex, dL, Angle, 1)
                If WaveLoopSequence(TrialIndex).FlashPeriodInFrames > 0 Then
                    Call ChangeRamp(WaveLoopSequence(TrialIndex).Mean, Cos(cPi * FlashPhase / 180))
                End If
            
            End If
        Else
            If frmMain.optFlashType(0).Value = True Then
                FillBkgSurf (MakeColor(WaveLoopSequence(TrialIndex).Mean * (1 - WaveLoopSequence(TrialIndex).Contrast)))
'                ddsBack.DrawCircle ResolutionX / 2 + WaveLoopSequence(TrialIndex).StimXposition * AnglToPixelCoeff, ResolutionY / 2 + WaveLoopSequence(TrialIndex).StimYPosition * AnglToPixelCoeff, WaveLoopSequence(TrialIndex).StimDiameter * AnglToPixelCoeff / 2 + 400
                Call CreatePictureWaveDirect(TrialIndex, dL, Angle, 1, True)
            End If
            If frmMain.optFlashType(1).Value = True Then
                FillBkgSurf (MakeColor(WaveLoopSequence(TrialIndex).Mean))
'                ddsBack.DrawCircle ResolutionX / 2 + WaveLoopSequence(TrialIndex).StimXposition * AnglToPixelCoeff, ResolutionY / 2 + WaveLoopSequence(TrialIndex).StimYPosition * AnglToPixelCoeff, WaveLoopSequence(TrialIndex).StimDiameter * AnglToPixelCoeff / 2 + 400
                Call CreatePictureWaveDirect(TrialIndex, dL, Angle, 1, True)
            End If
            If frmMain.optFlashType(2).Value = True Then
                If frmMain.optWaveType(0).Value = True Then
                    FillBkgSurf (MakeColor(WaveLoopSequence(TrialIndex).Mean * (1 + WaveLoopSequence(TrialIndex).Contrast)))
                    Call CreatePictureWaveDirect(TrialIndex, dL, Angle, -1, False)
                End If
                If frmMain.optWaveType(1).Value = True Then
                    Call CreatePictureWaveByBlt(TrialIndex, dL, Angle, 1)
                    If WaveLoopSequence(TrialIndex).FlashPeriodInFrames > 0 Then
                        Call ChangeRamp(WaveLoopSequence(TrialIndex).Mean, Cos(cPi * FlashPhase / 180))
                    End If
                End If
            End If
        End If
        
        If TrIndStep = 2 Then
            CurrentLUTStart = WaveLoopSequence(TrialIndex + 1).LUTStart
            CurrentLUTEnd = WaveLoopSequence(TrialIndex + 1).LUTEnd

            Call CalculateParamByFrame(TrialIndex + 1, Frame, dL, Angle, BarIsOn, FlashPhase)
            Call ChangeRamp(WaveLoopSequence(TrialIndex + 1).Mean, Sin(cPi * FlashPhase / 180))
            If BarIsOn Then
                Call CreatePictureWaveDirect(TrialIndex + 1, dL, Angle, 1, False)
            Else
                If frmMain.optFlashType(2).Value = True Then
                    Call CreatePictureWaveDirect(TrialIndex + 1, dL, Angle, -1, False)
                End If
            End If
        End If
    
        Trash = QueryPerformanceCounter(stopTime)
        frmMain.lblDrT.Caption = Format((stopTime - startTime) * 1000# / tikPerSec, "###0.00000")
        
        Frame = Frame + 1

        If frmMain.optWaveType(1).Value = True And WaveLoopSequence(TrialIndex).FlashPeriodInFrames > 0 Then
            'seems to me setGammaRamp itself syncronized with vertical blank
            'FlipNoVSync
          '  ddsFront.Flip Nothing, DDFLIP_DONOTWAIT
 '''         WaitFrameStart
            If Frame = 1 Then
                Flip
            Else
                If WaveLoopSequence(TrialIndex).WaveTempPeriodInFrames > 0 Then
                    WaitFrameStart
                    Flip
                End If
            End If
 '           Trash = QueryPerformanceCounter(stopTime2)
 '           If 1 / FrameRate > CDbl((stopTime2 - startTime) / tikPerSec) Then
                WaitFrameStart
 '           End If
        Else
''            Flip
''            WaitFrameEnd
 ''           WaitFrameStart
            
            WaitFrameEnd
            FlipNoVSync
 '''           WaitFrameEnd

        End If

        
        'reset all TTL outs to initial state if puls type=0 (short puls)
        For i = 0 To MarkNumber - 1
            If MarkArray(i).MarkType = 0 Then
                Call SetTTLMark(MarkArray(i).InitialState, MarkArray(i).MarkCh)
            End If
        Next i
        Trash = QueryPerformanceCounter(stopTime2)
        If 1.5 / FrameRate < CDbl((stopTime2 - startTime) / tikPerSec) Then
            scFrames = scFrames + 1
            frmMain.lblMissedFrames.Caption = CStr(scFrames)
            Frame = Frame + 1 'compensation for missed frame
        End If
        
''        If Frame = 1 Then
''        frmMain.lbl1.Caption = Format((stopTime2 - startTime) * 1000# / tikPerSec, "###0.00000")
''        End If
       
        If MaxdT < (CDbl((stopTime2 - startTime) / tikPerSec) - 1 / FrameRate) * 1000# Then
            MaxdT = (CDbl((stopTime2 - startTime) / tikPerSec) - 1 / FrameRate) * 1000#
                frmMain.lblMaxdT.Caption = Format(MaxdT, "###0.00000")
        End If
        TotaldT = TotaldT + (CDbl((stopTime2 - startTime) / tikPerSec) - 1 / FrameRate) * 1000#
        TotaldT = (CDbl((stopTime2 - startTime) / tikPerSec)) * 1000#

        frmMain.lblTotaldT.Caption = Format(TotaldT, "###0.00000")
        
        startTime = stopTime2
        DoEvents
        If bStop = True Then
            Exit For
        End If
    Loop

    If bStop = True Then
        Exit For
    End If
    
    If loopsizeToAVG > 0 Then
    NewCounterReadingAvg = ReadCounter(CounterNforAvg)
    If CurrentCounterReadingAvg > NewCounterReadingAvg Then
        TotalCounts(avgIndex * 2 + 1) = TotalCounts(avgIndex * 2 + 1) + NewCounterReadingAvg
    Else
        TotalCounts(avgIndex * 2 + 1) = TotalCounts(avgIndex * 2 + 1) + NewCounterReadingAvg - CurrentCounterReadingAvg
    End If
    Trash = QueryPerformanceCounter(avgstopTime)
    TotalTime(avgIndex * 2 + 1) = TotalTime(avgIndex * 2 + 1) + CDbl((avgstopTime - avgstartTime) / tikPerSec)
    AvgCounts(avgIndex * 2 + 1) = TotalCounts(avgIndex * 2 + 1) / TotalTime(avgIndex * 2 + 1)
    CurrentCounterReading = NewCounterReading
    avgstartTime = avgstopTime
    ShowAvg
    End If
Next TrialIndex

'reset all TTL outs to initial state
For i = 0 To MarkNumber - 1
    Call SetTTLMark(MarkArray(i).InitialState, MarkArray(i).MarkCh)
Next i
GCon.SetGammaRamp DDSGR_DEFAULT, InitialRamp

FillTwoSurf (MakeColor(WaveLoopSequence(0).BKGColor))
End Sub

Private Sub cmdSaveCounts_Click()
Dim FileName As String
On Error GoTo Err:
    frmCDHolder.CommonDialog1.FILTER = "Counts file (*.cnt)|*.cnt|All files(*.*)|*.*"
    frmCDHolder.CommonDialog1.FileName = ""
    frmCDHolder.CommonDialog1.ShowSave
    FileName = frmCDHolder.CommonDialog1.FileName
    SaveAvg (FileName)
Err:
    
End Sub

Private Sub cmdStop_Click()
    bStop = True
End Sub

Private Sub Form_Activate()
    If needsRestore = True Then
        frmStimulus.WindowState = 2
    End If

End Sub

Private Sub Form_KeyDown(KeyCode As Integer, Shift As Integer)


Select Case KeyCode
Case vbKeyEscape
    bStop = True
Case vbKeyUp
    StimYPosition = StimYPosition - 0.1
Case vbKeyDown
    StimYPosition = StimYPosition + 0.1
Case vbKeyLeft
    StimXposition = StimXposition - 0.1
Case vbKeyRight
    StimXposition = StimXposition + 0.1
Case vbKeyN
    Exit Sub
End Select
    WaveDefault.StimXposition = StimXposition
    WaveDefault.StimYPosition = StimYPosition
    'WaveLoopSequence(0) = WaveDefault
    WaveLoopSequence(0).StimXposition = StimXposition
    WaveLoopSequence(0).StimYPosition = StimYPosition
    
frmMain.lblStimPosX.Caption = Format(WaveLoopSequence(0).StimXposition, "###0.00")
frmMain.lblStimPosY.Caption = Format(WaveLoopSequence(0).StimYPosition, "###0.00")
    
End Sub

Private Sub Form_Load()
Dim i As Integer

hProcess = GetCurrentProcess()
hThread = GetCurrentThread()
'SetPriorityClass hProcess, HIGH_PRIORITY_CLASS
SetPriorityClass hProcess, REALTIME_PRIORITY_CLASS
SetThreadPriority hThread, 2
Load frmStimulus


IniParam
ReadIni (App.Path & "\" & IniFileName)
AnglToPixelCoeff = (cPi / 180) * Distance * ResolutionX / dHMonSize

IniBoard



For i = 0 To MarkNumber - 1
    Call SetTTLMark(MarkArray(i).InitialState, MarkArray(i).MarkCh)
Next i


'Call PositionForms(frmMain, 2, frmStimulus, 3)
'Call PositionForm(frmMain, 2, 10, 10)
'Call PositionForm(frmStimulus, 3, 10, 10)


If StimMonitor = 2 Then
    Call PositionForm(frmMain, 2, 10, 10)
    Call PositionForm(frmStimulus, 3, 10, 10)
Else
    Call PositionForm(frmMain, 3, 10, 10)
    Call PositionForm(frmStimulus, 2, 10, 10)
End If

frmMain.WindowState = 2
frmStimulus.WindowState = 2
frmStimulus.Show


    ResolutionX = frmStimulus.ScaleWidth
    ResolutionY = frmStimulus.ScaleHeight

frmMain.lblResolution.Caption = CStr(ResolutionX) & " x " & CStr(ResolutionY)
Timer1.Enabled = True
End Sub


Private Sub Form_QueryUnload(Cancel As Integer, UnloadMode As Integer)
    DestroySurface
    RestoreDispMode
Unload frmStimulus
Set objDX = Nothing
Set objDD = Nothing
End

End Sub



Private Sub optWaveType_Click(Index As Integer)
If Index = 0 Then
    frmMain.txtRotationPeriodInHz.Enabled = True
    frmMain.chkDoubleStim.Enabled = True
End If
If Index = 1 Then
    frmMain.txtRotationPeriodInHz.Enabled = False
    frmMain.chkDoubleStim.Enabled = False
End If
End Sub



Private Sub Timer1_Timer()

 '   SetDispMode
 If StimMonitor = 2 Then
    SetFullScreenDispMode frmStimulus, 3
 Else
    SetFullScreenDispMode frmStimulus, 2
 End If
    CreateMainSurfacesWithClipper frmStimulus
    CreateSecondarySurfaces 3, 3

    Timer1.Enabled = False
    frmMain.WindowState = 2
End Sub


Private Sub txtOrientation_LostFocus()
    If IsNumeric(txtOrientation.Text) Then
        WaveLoopSequence(0).Orientation = CDbl(txtOrientation.Text) * cPi / 180
    End If
    ShowParametersW WaveLoopSequence(0)
End Sub


Private Sub txtBarW_LostFocus()
    If IsNumeric(txtBarW.Text) Then
        WaveLoopSequence(0).BarW = CDbl(txtBarW.Text)
    End If
    ShowParametersW WaveLoopSequence(0)
End Sub

Private Sub txtBKGColor_LostFocus()
    If IsNumeric(txtBKGColor.Text) Then
        WaveLoopSequence(0).BKGColor = CDbl(txtBKGColor.Text)
    End If
    ShowParametersW WaveLoopSequence(0)

End Sub

Private Sub txtContrast_LostFocus()

    If IsNumeric(txtContrast.Text) Then
        WaveLoopSequence(0).Contrast = CDbl(txtContrast.Text)
    End If

    If WaveLoopSequence(0).Contrast > 1 Then
        WaveLoopSequence(0).Contrast = 1
    End If
    ShowParametersW WaveLoopSequence(0)
End Sub

Private Sub txtDistance_Change()
    If IsNumeric(txtDistance.Text) Then
        Distance = CDbl(txtDistance.Text)
        AnglToPixelCoeff = (cPi / 180) * Distance * ResolutionX / dHMonSize
    End If
    txtDistance.Text = CStr(Distance)
End Sub


Private Sub txtFlashPeriodInHz_LostFocus()
    If IsNumeric(txtFlashPeriodInHz.Text) Then
        WaveLoopSequence(0).FlashPeriodInHz = CDbl(txtFlashPeriodInHz.Text)
    End If
RecalculateTempPeriod
ShowParametersW WaveLoopSequence(0)
End Sub

Private Sub txtFrameRate_LostFocus()
    If IsNumeric(txtFrameRate.Text) Then
        FrameRate = CStr(txtFrameRate.Text)
        RecalculateTempPeriod
    End If

ShowParametersW WaveLoopSequence(0)

End Sub

Private Sub txtFramesBeforeStim_LostFocus()
    If IsNumeric(txtFramesBeforeStim.Text) Then
        FramesBeforeStim = CInt(txtFramesBeforeStim.Text)
    End If
    txtFramesBeforeStim.Text = CStr(FramesBeforeStim)

End Sub

Private Sub txtFramesOfStim_LostFocus()
    If IsNumeric(txtFramesOfStim.Text) Then
        FramesOfStim = CInt(txtFramesOfStim.Text)
    End If
    txtFramesOfStim.Text = CStr(FramesOfStim)

End Sub

Private Sub txtHSize_Change()
    If IsNumeric(txtHSize.Text) Then
        dHMonSize = CDbl(txtHSize.Text)
        AnglToPixelCoeff = (cPi / 180) * Distance * ResolutionX / dHMonSize
    End If
    txtHSize.Text = CStr(dHMonSize)
End Sub



Private Sub txtMean_LostFocus()
    If IsNumeric(txtMean.Text) Then
        WaveLoopSequence(0).Mean = CDbl(txtMean.Text)
    End If
    ShowParametersW WaveLoopSequence(0)

End Sub




Private Sub txtPhase_LostFocus()
    If IsNumeric(txtPhase.Text) Then
        WaveLoopSequence(0).Phase = CDbl(txtPhase.Text) * cPi / 180
    End If
    ShowParametersW WaveLoopSequence(0)

End Sub



Private Sub txtRotationPeriodInHz_LostFocus()
    If IsNumeric(txtRotationPeriodInHz.Text) Then
        WaveLoopSequence(0).RotationPeriodInHz = CDbl(txtRotationPeriodInHz.Text)
    End If
RecalculateTempPeriod
ShowParametersW WaveLoopSequence(0)
End Sub

Private Sub txtSpPeriod_LostFocus()
    If IsNumeric(txtSpPeriod.Text) Then
        WaveLoopSequence(0).WaveSpPeriod = CDbl(txtSpPeriod.Text)
    End If

    ShowParametersW WaveLoopSequence(0)
End Sub

Private Sub txtTempPeriodInHz_LostFocus()
If IsNumeric(txtTempPeriodInHz.Text) Then
'    If CDbl(txtTempPeriodInHz.Text) > 0 Then
        WaveLoopSequence(0).WaveTempPeriodInHz = CDbl(txtTempPeriodInHz.Text)
'    End If
End If
RecalculateTempPeriod
ShowParametersW WaveLoopSequence(0)
End Sub

Private Sub vsScreen1_KeyDown(KeyCode As Integer, Shift As Integer)
Select Case KeyCode
Case vbKeyEscape
bStop = True
End Select

End Sub

