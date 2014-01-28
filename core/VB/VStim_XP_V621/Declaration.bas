Attribute VB_Name = "Declaration"
Option Explicit

Public ResolutionX As Integer
Public ResolutionY As Integer

Public Type WParameters
     StimDiameter As Double
     StimLength As Double
     WaveSpPeriod As Double
     Orientation As Double
     Mean As Double
     Phase As Double
     WaveTempPeriodInHz As Double
     WaveTempPeriodInFrames As Long
     Contrast As Double
     StimXposition As Double
     StimYPosition As Double
     PeriodsToShow As Long
     BarW As Double
     RotationPeriodInHz As Double
     RotationPeriodInFrames As Long
     BKGColor As Double
     FlashPeriodInHz As Double
     FlashPeriodInFrames As Long
     FlashType As Integer
     LUTStart As Integer
     LUTEnd As Integer
End Type

Public CounterNforAvg As Long
Public loopsizeToAVG As Integer
Public TotalCounts() As Long
Public AvgCounts() As Double
Public TotalTime() As Double

Public WaveLoopSequence() As WParameters
Public WaveDefault As WParameters
Public LoopsNumber As Integer
Public LoopSizes() As Integer   ' size of each loop
Public LoopSequenceSize As Long ' total size of WaveLoopSequence

Public FrameRate As Double
Public AnglToPixelCoeff As Double
Public dHMonSize As Double
Public Distance As Double

Public StimXposition As Double
Public StimYPosition As Double

Public Declare Function GetPrivateProfileString Lib "kernel32" Alias "GetPrivateProfileStringA" (ByVal lpApplicationName As String, ByVal lpKeyName As Any, ByVal lpDefault As String, ByVal lpReturnedString As String, ByVal nSize As Long, ByVal lpFileName As String) As Long
Public Declare Function WritePrivateProfileString Lib "kernel32" Alias "WritePrivateProfileStringA" (ByVal lpApplicationName As String, ByVal lpKeyName As Any, ByVal lpString As Any, ByVal lpFileName As String) As Long


Public bStop As Boolean

Public Declare Function QueryPerformanceCounter Lib "kernel32" (lpPerformanceCount As Currency) As Long
Public Declare Function QueryPerformanceFrequency Lib "kernel32" (lpFrequency As Currency) As Long

Public Declare Function BitBlt Lib "gdi32" (ByVal hDestDC As Long, ByVal x As Long, ByVal y As Long, ByVal nWidth As Long, ByVal nHeight As Long, ByVal hSrcDC As Long, ByVal xSrc As Long, ByVal ySrc As Long, ByVal dwRop As Long) As Long
Global Const cPi = 3.1415927


'Public Const COUNTER_NUM = 1
'Public Const BOARD_NUM = 0


'I/O parameters

Public BoardN As Integer
Public CounterN As Long
Public CounterNTrigger As Long
Public CounterNPol As Long
Public CounterNTriggerPol As Long
Public FramesBeforeStim As Integer
Public FramesOfStim As Integer



Public Type TTLMark
     MarkCh As Integer
     MarkType As Integer    '0-short puls, 1 -square wave
     MarkPeriod As Integer  ' if MarkType =1
     MarkTiming As String  ' Start, DriftPeriod, RotationPeriod, FlashPeriod,
                            ' Loop0 (each new parameter set), Loop1, Loop2, Loop3...
     InitialState As Integer
     Position As Integer '0 - before blank, 1 - before stim
End Type


Public MarkArray() As TTLMark
Public MarkNumber As Integer


Public ParamFileName As String

'LUT table
Public iColorLUT() As Integer
Public LUTLen As Long
Public LUTFname As String
Public CurrentLUTStart As Integer
Public CurrentLUTEnd As Integer


Public objDX As New DirectX7
Public objDD As DirectDraw7

Public Declare Function GetCurrentProcess Lib "kernel32" () As Long
Public Declare Function SetPriorityClass Lib "kernel32" (ByVal hProcess As Long, ByVal dwPriorityClass As Long) As Long
Public Declare Function GetCurrentThread Lib "kernel32" () As Long
Public Declare Function SetThreadPriority Lib "kernel32" (ByVal hThread As Long, ByVal nPriority As Long) As Long

Public Const HIGH_PRIORITY_CLASS = &H80
Public Const REALTIME_PRIORITY_CLASS = &H100


Public hProcess As Long
Public hThread As Long


Public Const IniFileName = "VSim.ini"

'*********************************************************************************************************************


'Declare Function Ellipse Lib "gdi32" (ByVal hdc As Long, ByVal X1 As Long, ByVal Y1 As Long, ByVal X2 As Long, ByVal Y2 As Long) As Long
'Declare Function CreatePen Lib "gdi32" (ByVal nPenStyle As Long, ByVal nWidth As Long, ByVal crColor As Long) As Long
'Declare Function CreateSolidBrush Lib "gdi32" (ByVal crColor As Long) As Long
'Declare Function DeleteObject Lib "gdi32" (ByVal hObject As Long) As Long
'Declare Function SelectObject Lib "gdi32" (ByVal hdc As Long, ByVal hObject As Long) As Long


'Public Declare Function BitBlt Lib "gdi32" (ByVal hDestDC As Long, ByVal X As Long, ByVal Y As Long, ByVal nWidth As Long, ByVal nHeight As Long, ByVal hSrcDC As Long, ByVal xSrc As Long, ByVal ySrc As Long, ByVal dwRop As Long) As Long
'Public Declare Function StretchBlt Lib "gdi32" (ByVal hdc As Long, ByVal X As Long, ByVal Y As Long, ByVal nWidth As Long, ByVal nHeight As Long, ByVal hSrcDC As Long, ByVal xSrc As Long, ByVal ySrc As Long, ByVal nSrcWidth As Long, ByVal nSrcHeight As Long, ByVal dwRop As Long) As Long

'Public Declare Sub Sleep Lib "kernel32" (ByVal dwMilliseconds As Long)
'Public Declare Function ShowCursor Lib "user32" (ByVal bShow As Long) As Long


'conversion coeffs
'Public DegreeToPix As Double


Public tikPerSec As Currency
Public startTime As Currency
Public stopTime As Currency
Public stopTime2 As Currency


Public TotaldT As Double
Public MaxdT As Double



Public Type POINTAPI
        x As Long
        y As Long
End Type

Public objPos() As POINTAPI
Public objVel() As POINTAPI

Public StimMonitor As Integer
