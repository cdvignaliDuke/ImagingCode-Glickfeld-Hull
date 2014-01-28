Attribute VB_Name = "Misc"
Option Explicit

Public Sub CalculateCompPar()
    AnglToPixelCoeff = (cPi / 180) * Distance * ResolutionX / dHMonSize
End Sub

Public Sub IniBoard()
Dim i As Integer
Call cbC9513Init(BoardN, 1, 1, FREQ4, CBDISABLED, CBDISABLED, CBDISABLED)
Call cbC9513Config(BoardN, 1, NOGATE, POSITIVEEDGE, CTRINPUT1, CBDISABLED, LOADREG, RECYCLE, CBDISABLED, COUNTUP, ALWAYSLOW)
Call cbC9513Config(BoardN, 2, NOGATE, POSITIVEEDGE, CTRINPUT2, CBDISABLED, LOADREG, RECYCLE, CBDISABLED, COUNTUP, ALWAYSLOW)
Call cbC9513Config(BoardN, 3, NOGATE, POSITIVEEDGE, CTRINPUT3, CBDISABLED, LOADREG, RECYCLE, CBDISABLED, COUNTUP, ALWAYSLOW)

If CounterNTriggerPol = -1 Then
    Call cbC9513Config(BoardN, CounterNTrigger, NOGATE, NEGATIVEEDGE, CounterNTrigger, CBDISABLED, LOADREG, RECYCLE, CBDISABLED, COUNTUP, ALWAYSLOW)
End If

If CounterNPol = -1 Then
    Call cbC9513Config(BoardN, CounterN, NOGATE, NEGATIVEEDGE, CounterN, CBDISABLED, LOADREG, RECYCLE, CBDISABLED, COUNTUP, ALWAYSLOW)
End If

End Sub


Public Sub ReadParamW(ParamFileName As String)
' reads parameters from .x file for wave
Dim RetString As String
Dim Size As Long
Dim RetSize As Integer
Dim i As Integer

Dim NameString As String

Dim CarVarName As String
Dim j As Long
Dim k As Integer
Dim m As Integer
Dim NameString2 As String
Dim LocalVarInCurrentLoop As Double
Dim LoopsSizeBefore As Long

Dim BarWProblem As Boolean
BarWProblem = False

RetString = Space(255)
Size = Len(RetString)
'1
'1- read all parameters, if parameter was not in .x file - keep old parameter
RetSize = GetPrivateProfileString("parameters", "StimDiameter", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
WaveDefault.StimDiameter = CDbl(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "StimLength", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
WaveDefault.StimLength = CDbl(RetString)
Else
WaveDefault.StimLength = WaveDefault.StimDiameter
End If

RetSize = GetPrivateProfileString("parameters", "WaveSpPeriod", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
WaveDefault.WaveSpPeriod = CDbl(RetString)
End If
RetSize = GetPrivateProfileString("parameters", "Orientation", "1000", RetString, Size, ParamFileName)
If CSng(RetString) < 900 Then
    WaveDefault.Orientation = CDbl(RetString) * cPi / 180
End If

RetSize = GetPrivateProfileString("parameters", "Mean", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
WaveDefault.Mean = CDbl(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "Phase", "1000", RetString, Size, ParamFileName)
If CSng(RetString) < 900 Then
WaveDefault.Phase = CDbl(RetString) * cPi / 180
End If

RetSize = GetPrivateProfileString("parameters", "WaveTempPeriodInHz", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
WaveDefault.WaveTempPeriodInHz = CDbl(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "Contrast", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
WaveDefault.Contrast = CDbl(RetString)

End If

RetSize = GetPrivateProfileString("parameters", "StimXPosition", "1000", RetString, Size, ParamFileName)
If CSng(RetString) < 900 Then
WaveDefault.StimXposition = CSng(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "StimYPosition", "1000", RetString, Size, ParamFileName)
If CSng(RetString) < 900 Then
WaveDefault.StimYPosition = CSng(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "PeriodsToShow", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
WaveDefault.PeriodsToShow = CInt(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "FramesBeforeStim", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
    FramesBeforeStim = CInt(RetString)
    frmMain.txtFramesBeforeStim = CStr(FramesBeforeStim)
End If


RetSize = GetPrivateProfileString("parameters", "FramesOfStim", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
    FramesOfStim = CInt(RetString)
    frmMain.txtFramesOfStim = CStr(FramesOfStim)
End If


RetSize = GetPrivateProfileString("parameters", "BarWidth", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
WaveDefault.BarW = CDbl(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "RotationPeriodInHz", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
WaveDefault.RotationPeriodInHz = CDbl(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "BKGColor", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
WaveDefault.BKGColor = CDbl(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "FlashPeriodInHz", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
WaveDefault.FlashPeriodInHz = CDbl(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "LUTStart", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
WaveDefault.LUTStart = CInt(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "LUTEnd", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
WaveDefault.LUTEnd = CInt(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "WaveType", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
    frmMain.optWaveType(CInt(RetString)).Value = True
End If

RetSize = GetPrivateProfileString("parameters", "ModulationType", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
    frmMain.optFlashType(CInt(RetString)).Value = True
End If


RetSize = GetPrivateProfileString("parameters", "MaskShape", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
    frmMain.chkSquareMask.Value = CInt(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "MaskOrientation", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
    frmMain.chkUseStimOri.Value = CInt(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "LockPhaseToMask", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
    frmMain.chkLockPhaseToMask.Value = CInt(RetString)
End If




If Abs(WaveDefault.WaveSpPeriod / 2 - WaveDefault.BarW) > 0.01 Then
    BarWProblem = True
End If


RetSize = GetPrivateProfileString("online", "CounterNforAvg", "-1", RetString, Size, ParamFileName)
If CLng(RetString) >= 0 Then
    CounterNforAvg = CLng(RetString)
End If

RetSize = GetPrivateProfileString("online", "loopsizeToAVG", "-1", RetString, Size, ParamFileName)
If CInt(RetString) >= 0 Then
    loopsizeToAVG = CInt(RetString)
End If
If loopsizeToAVG > 0 Then
    ReDim TotalCounts(loopsizeToAVG * 2 - 1) As Long
    ReDim AvgCounts(loopsizeToAVG * 2 - 1) As Double
    ReDim TotalTime(loopsizeToAVG * 2 - 1) As Double
End If


LoopSequenceSize = 1
RetSize = GetPrivateProfileString("loops", "TotalNumber", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
    LoopsNumber = CInt(RetString)
Else
    ReDim WaveLoopSequence(0) As WParameters
    WaveLoopSequence(0) = WaveDefault
    Exit Sub
End If

ReDim LoopSizes(LoopsNumber - 1) As Integer
For i = 1 To LoopsNumber
    NameString = "loop" & CStr(i) & "_0"
    RetSize = GetPrivateProfileString(NameString, "loopsize", "-1", RetString, Size, ParamFileName)
    If CSng(RetString) >= 0 Then
        LoopSizes(i - 1) = CInt(RetString)
        LoopSequenceSize = LoopSequenceSize * LoopSizes(i - 1)
    End If
Next i


ReDim WaveLoopSequence(LoopSequenceSize - 1) As WParameters

For i = 0 To LoopSequenceSize - 1
    WaveLoopSequence(i) = WaveDefault
Next i
LoopsSizeBefore = 1

    For i = 1 To LoopsNumber
        For m = 0 To 9
            NameString = "loop" & CStr(i) & "_" & CStr(m)
            RetSize = GetPrivateProfileString(NameString, "varname", "BlaBla", RetString, Size, ParamFileName)
            'If InStr(ExpType, "BlaBla") > 0 Then
            'Else
                CarVarName = Left(RetString, RetSize)
            'End If

            For k = 1 To LoopSizes(i - 1)
                NameString2 = "val" & CStr(k)
                RetSize = GetPrivateProfileString(NameString, NameString2, "-10000", RetString, Size, ParamFileName)
                'MsgBox NameString & NameString2 & RetString
                If CDbl(RetString) > -9000 Then
                    LocalVarInCurrentLoop = CDbl(RetString)
                End If

                Select Case CarVarName
                Case "StimDiameter"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).StimDiameter = LocalVarInCurrentLoop
                        'MsgBox CStr(LocalVarInCurrentLoop)
                        End If
                    Next j
                Case "StimLength"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).StimLength = LocalVarInCurrentLoop
                        'MsgBox CStr(LocalVarInCurrentLoop)
                        End If
                    Next j
                Case "WaveSpPeriod"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).WaveSpPeriod = LocalVarInCurrentLoop
                        End If
                    Next j
                Case "Orientation"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).Orientation = LocalVarInCurrentLoop * cPi / 180
                        End If
                    Next j
                Case "Mean"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).Mean = LocalVarInCurrentLoop
                        End If
                    Next j
                Case "Phase"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).Phase = LocalVarInCurrentLoop * cPi / 180
                        End If
                    Next j
                Case "WaveTempPeriodInHz"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).WaveTempPeriodInHz = LocalVarInCurrentLoop
                        End If
                    Next j
                Case "Contrast"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).Contrast = LocalVarInCurrentLoop
                        End If
                    Next j
                Case "StimXPosition"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).StimXposition = LocalVarInCurrentLoop
                        End If
                    Next j
                Case "StimYPosition"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).StimYPosition = LocalVarInCurrentLoop
                        End If
                    Next j
                Case "PeriodsToShow"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).PeriodsToShow = LocalVarInCurrentLoop
                        End If
                    Next j
                Case "BarWidth"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).BarW = LocalVarInCurrentLoop
                            If Abs(WaveLoopSequence(j).WaveSpPeriod / 2 - WaveLoopSequence(j).BarW) > 0.01 Then
                                BarWProblem = True
                            End If
                        End If
                    Next j
                    
                Case "RotationPeriodInHz"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).RotationPeriodInHz = LocalVarInCurrentLoop
                        End If
                    Next j
                Case "BKGColor"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).BKGColor = LocalVarInCurrentLoop
                        End If
                    Next j
                Case "FlashPeriodInHz"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).FlashPeriodInHz = LocalVarInCurrentLoop
                        End If
                    Next j
                Case "LUTStart"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).LUTStart = LocalVarInCurrentLoop
                        End If
                    Next j
                Case "LUTEnd"
                    For j = 0 To LoopSequenceSize - 1
                        If k = (j \ LoopsSizeBefore) Mod LoopSizes(i - 1) + 1 Then
                            WaveLoopSequence(j).LUTEnd = LocalVarInCurrentLoop
                        End If
                    Next j
                End Select
            Next k
        Next m
        LoopsSizeBefore = LoopsSizeBefore * LoopSizes(i - 1)
    Next i

' FramesBeforeStim, FramesOfStim are ini general .ini and here too
RetSize = GetPrivateProfileString("parameters", "FramesBeforeStim", "-1", RetString, Size, IniFileName)
If CSng(RetString) >= 0 Then
    FramesBeforeStim = CInt(RetString)
    frmMain.txtFramesBeforeStim = CStr(FramesBeforeStim)
End If


RetSize = GetPrivateProfileString("parameters", "FramesOfStim", "-1", RetString, Size, IniFileName)
If CSng(RetString) >= 0 Then
    FramesOfStim = CInt(RetString)
    frmMain.txtFramesOfStim = CStr(FramesOfStim)
End If

RetSize = GetPrivateProfileString("parameters", "LUTfile", "", RetString, Size, IniFileName)
If CSng(RetSize) > 0 Then
    LUTFname = CInt(RetString)
    ReadLUT LUTFname
End If


'****************

RetSize = GetPrivateProfileString("TTLout", "TotalNumber", "-1", RetString, Size, ParamFileName)
If CSng(RetString) >= 0 Then
    MarkNumber = CInt(RetString)
    ReDim MarkArray(MarkNumber - 1) As TTLMark
    For i = 1 To MarkNumber
        MarkArray(i - 1).MarkTiming = "Start"
        MarkArray(i - 1).MarkCh = 0
        MarkArray(i - 1).MarkType = 0
        MarkArray(i - 1).MarkPeriod = 0
        MarkArray(i - 1).InitialState = 0
        MarkArray(i - 1).Position = 0
    
        NameString = "Mark" & CStr(i)
        
        RetSize = GetPrivateProfileString(NameString, "MarkTiming", "Start", RetString, Size, ParamFileName)
        MarkArray(i - 1).MarkTiming = Trim(RetString)
        
        RetSize = GetPrivateProfileString(NameString, "MarkCh", "-1", RetString, Size, ParamFileName)
        If CSng(RetString) >= 0 Then
            MarkArray(i - 1).MarkCh = CInt(RetString)
        End If
        
        RetSize = GetPrivateProfileString(NameString, "MarkType", "-1", RetString, Size, ParamFileName)
        If CSng(RetString) >= 0 Then
            MarkArray(i - 1).MarkType = CInt(RetString)
        End If
        RetSize = GetPrivateProfileString(NameString, "MarkPeriod", "-1", RetString, Size, ParamFileName)
        If CSng(RetString) >= 0 Then
            MarkArray(i - 1).MarkPeriod = CInt(RetString)
        End If
        
        RetSize = GetPrivateProfileString(NameString, "InitialState", "-1", RetString, Size, ParamFileName)
        If CSng(RetString) >= 0 Then
            MarkArray(i - 1).InitialState = CInt(RetString)
        End If
        RetSize = GetPrivateProfileString(NameString, "Position", "-1", RetString, Size, ParamFileName)
        If CSng(RetString) >= 0 Then
            MarkArray(i - 1).Position = CInt(RetString)
        End If

    Next i
End If

For i = 0 To MarkNumber - 1
    Call SetTTLMark(MarkArray(i).InitialState, MarkArray(i).MarkCh)
Next i


'********************


RecalculateTempPeriod

ShowParametersW WaveLoopSequence(0)



If BarWProblem = True Then
    MsgBox "Check BarWidth"
End If

If CounterNforAvg = CounterN Then
    MsgBox "Check CounterNforAvg"
End If

If CounterNforAvg = CounterNTrigger Then
    MsgBox "Check CounterNforAvg"
End If

End Sub

Public Sub ReadIni(IniFileName As String)

' reads parameters from .x file for wave
Dim RetString As String
Dim Size As Long
Dim RetSize As Integer
Dim i As Integer
Dim NameString As String


RetString = Space(255)
Size = Len(RetString)

RetSize = GetPrivateProfileString("parameters", "FrameRate", "-1", RetString, Size, IniFileName)

If CSng(RetString) >= 0 Then
    FrameRate = CDbl(RetString)
    frmMain.txtFrameRate = CStr(FrameRate)
End If

RetSize = GetPrivateProfileString("parameters", "HMonSize", "-1", RetString, Size, IniFileName)
If CSng(RetString) >= 0 Then
    dHMonSize = CDbl(RetString)
    frmMain.txtHSize = CStr(dHMonSize)
End If


RetSize = GetPrivateProfileString("parameters", "Distance", "-1", RetString, Size, IniFileName)
If CSng(RetString) >= 0 Then
    Distance = CDbl(RetString)
    frmMain.txtDistance = CStr(Distance)
End If


RetSize = GetPrivateProfileString("parameters", "BoardN", "-1", RetString, Size, IniFileName)
If CSng(RetString) >= 0 Then
    BoardN = CInt(RetString)
End If


RetSize = GetPrivateProfileString("parameters", "CounterN", "-1", RetString, Size, IniFileName)
If CLng(RetString) >= 0 Then
    CounterN = CLng(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "CounterNTrigger", "-1", RetString, Size, IniFileName)
If CLng(RetString) >= 0 Then
    CounterNTrigger = CLng(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "CounterNPol", "0", RetString, Size, IniFileName)
If CLng(RetString) <> 0 Then
    CounterNPol = CLng(RetString)
End If

RetSize = GetPrivateProfileString("parameters", "CounterNTriggerPol", "0", RetString, Size, IniFileName)
If CLng(RetString) <> 0 Then
    CounterNTriggerPol = CLng(RetString)
End If



RetSize = GetPrivateProfileString("parameters", "FramesBeforeStim", "-1", RetString, Size, IniFileName)
If CSng(RetString) >= 0 Then
    FramesBeforeStim = CInt(RetString)
    frmMain.txtFramesBeforeStim = CStr(FramesBeforeStim)
End If


RetSize = GetPrivateProfileString("parameters", "FramesOfStim", "-1", RetString, Size, IniFileName)
If CSng(RetString) >= 0 Then
    FramesOfStim = CInt(RetString)
    frmMain.txtFramesOfStim = CStr(FramesOfStim)
End If

RetSize = GetPrivateProfileString("parameters", "StimMonitor", "-1", RetString, Size, IniFileName)
If CSng(RetString) >= 0 Then
    StimMonitor = CInt(RetString)
End If

'****************

RetSize = GetPrivateProfileString("TTLout", "TotalNumber", "-1", RetString, Size, IniFileName)

If CSng(RetString) >= 0 Then
    MarkNumber = CInt(RetString)
    ReDim MarkArray(MarkNumber - 1) As TTLMark
    For i = 1 To MarkNumber
        MarkArray(i - 1).MarkTiming = "Start"
        MarkArray(i - 1).MarkCh = 0
        MarkArray(i - 1).MarkType = 0
        MarkArray(i - 1).MarkPeriod = 0
        MarkArray(i - 1).InitialState = 0
    
        NameString = "Mark" & CStr(i)
        
        RetSize = GetPrivateProfileString(NameString, "MarkTiming", "Start", RetString, Size, IniFileName)
        MarkArray(i - 1).MarkTiming = Trim(RetString)
        
        RetSize = GetPrivateProfileString(NameString, "MarkCh", "-1", RetString, Size, IniFileName)
        If CSng(RetString) >= 0 Then
            MarkArray(i - 1).MarkCh = CInt(RetString)
        End If
        
        RetSize = GetPrivateProfileString(NameString, "MarkType", "-1", RetString, Size, IniFileName)
        If CSng(RetString) >= 0 Then
            MarkArray(i - 1).MarkType = CInt(RetString)
        End If
        RetSize = GetPrivateProfileString(NameString, "MarkPeriod", "-1", RetString, Size, IniFileName)
        If CSng(RetString) >= 0 Then
            MarkArray(i - 1).MarkPeriod = CInt(RetString)
        End If
        RetSize = GetPrivateProfileString(NameString, "InitialState", "-1", RetString, Size, IniFileName)
        If CSng(RetString) >= 0 Then
            MarkArray(i - 1).InitialState = CInt(RetString)
        End If
        RetSize = GetPrivateProfileString(NameString, "Position", "-1", RetString, Size, IniFileName)
        If CSng(RetString) >= 0 Then
            MarkArray(i - 1).Position = CInt(RetString)
        End If

    Next i
End If

For i = 0 To MarkNumber - 1
    Call SetTTLMark(MarkArray(i).InitialState, MarkArray(i).MarkCh)
Next i


'********************



RecalculateTempPeriod
ShowParametersW WaveLoopSequence(0)


End Sub


Public Sub WriteIni(IniFileName As String)

WritePrivateProfileString "parameters", "FrameRate", CStr(FrameRate), IniFileName
WritePrivateProfileString "parameters", "HMonSize", CStr(dHMonSize), IniFileName
WritePrivateProfileString "parameters", "Distance", CStr(Distance), IniFileName
WritePrivateProfileString "parameters", "FramesBeforeStim", CStr(FramesBeforeStim), IniFileName
WritePrivateProfileString "parameters", "FramesOfStim", CStr(FramesOfStim), IniFileName

End Sub

Public Sub SaveCurrentParamStim(TrialIndex As Long, FileName As String)
Dim i As Integer
Dim string1 As String
i = InStr(FileName, ".")
If i > 0 Then
FileName = Left(FileName, i - 1) & ".x"
Else
FileName = FileName & ".x"
End If

string1 = "Run" & CStr(TrialIndex)
With WaveLoopSequence(TrialIndex)
    WritePrivateProfileString string1, "StimDiameter", CStr(.StimDiameter), FileName
    WritePrivateProfileString string1, "WaveSpPeriod", CStr(.WaveSpPeriod), FileName
    WritePrivateProfileString string1, "Orientation", CStr(.Orientation / (cPi / 180)), FileName
    WritePrivateProfileString string1, "Mean", CStr(.Mean), FileName
    WritePrivateProfileString string1, "Phase", CStr(.Phase / (cPi / 180)), FileName
    WritePrivateProfileString string1, "WaveTempPeriodInHz", CStr(.WaveTempPeriodInHz), FileName
    WritePrivateProfileString string1, "WaveTempPeriodInFrames", CStr(.WaveTempPeriodInFrames), FileName
    WritePrivateProfileString string1, "Contrast", CStr(.Contrast), FileName
    WritePrivateProfileString string1, "StimXposition", CStr(.StimXposition), FileName
    WritePrivateProfileString string1, "StimYPosition", CStr(.StimYPosition), FileName
    WritePrivateProfileString string1, "PeriodsToShow", CStr(.PeriodsToShow), FileName
    WritePrivateProfileString string1, "BarWidth", CStr(.BarW), FileName
    WritePrivateProfileString string1, "RotationPeriodInHz", CStr(.RotationPeriodInHz), FileName
    WritePrivateProfileString string1, "BKGColor", CStr(.BKGColor), FileName
    WritePrivateProfileString string1, "FlashPeriodInHz", CStr(.FlashPeriodInHz), FileName
End With

End Sub

Public Sub SetTTLMark(Mark As Integer, Ch As Integer)
Dim Trash As Long
'commented dor demo 08/06/02
If Ch >= 0 Then
    Trash = cbDBitOut(BoardN, AUXPORT, Ch, Mark Mod 2)
End If
End Sub


Public Sub IniParam()
Dim i As Integer
FrameRate = 60
dHMonSize = 41
Distance = 114
AnglToPixelCoeff = (cPi / 180) * Distance * ResolutionX / dHMonSize

StimXposition = 0
StimYPosition = 0
LUTLen = 2
ReDim iColorLUT(2, LUTLen - 1) As Integer
iColorLUT(0, 0) = 0
iColorLUT(0, 1) = 255
iColorLUT(1, 0) = 0
iColorLUT(1, 1) = 255
iColorLUT(2, 0) = 0
iColorLUT(2, 1) = 255


WaveDefault.StimDiameter = 160
WaveDefault.StimLength = 120
WaveDefault.WaveSpPeriod = 5
WaveDefault.Orientation = 0
WaveDefault.Mean = 0.5
WaveDefault.Phase = 0
WaveDefault.WaveTempPeriodInHz = 1
WaveDefault.WaveTempPeriodInFrames = 170
WaveDefault.Contrast = 1
WaveDefault.PeriodsToShow = 30000
WaveDefault.StimXposition = 0
WaveDefault.StimYPosition = 0
WaveDefault.BKGColor = 0.5

WaveDefault.BarW = 1
WaveDefault.RotationPeriodInHz = 0.1
WaveDefault.FlashPeriodInHz = 0
WaveDefault.LUTStart = 1
WaveDefault.LUTEnd = 10000

CurrentLUTStart = 1
CurrentLUTEnd = 10000

LoopSequenceSize = 1

ReDim WaveLoopSequence(0) As WParameters
WaveLoopSequence(0) = WaveDefault
RecalculateTempPeriod
ShowParametersW WaveLoopSequence(0)

FramesBeforeStim = 10
FramesOfStim = 10

BoardN = 0
CounterN = 1
CounterNTrigger = 1
CounterNPol = 1
CounterNTriggerPol = 1


MarkNumber = 1
ReDim MarkArray(MarkNumber - 1) As TTLMark

MarkArray(0).MarkTiming = "Start"
MarkArray(0).MarkCh = 0
MarkArray(0).MarkType = 0
MarkArray(0).MarkPeriod = 0
MarkArray(0).InitialState = 0
MarkArray(0).Position = 0


StimMonitor = 1

loopsizeToAVG = -1

End Sub

Public Sub ReadLUT(LUTFilename As String)
'sFile - name of file contains LUT - three ASCII columns
Dim lFileL As Long
Dim i As Long
Dim OneLine As String
LUTLen = 0
    Open LUTFilename For Input As #1
    'this loop is to find the size of lut file
    Do While Not EOF(1) ' Loop until end of file.
        Line Input #1, OneLine ' Read line into variable.
        LUTLen = LUTLen + 1
    Loop
    Close #1
    ReDim iColorLUT(2, LUTLen - 1) As Integer
    
    Open LUTFilename For Input As #1
    For i = 0 To LUTLen - 1
        Line Input #1, OneLine ' Read line into variable.
    iColorLUT(0, i) = CInt(FirstTocken(OneLine))
    iColorLUT(1, i) = CInt(FirstTocken(OneLine))
    iColorLUT(2, i) = CInt(FirstTocken(OneLine))
    Next i
    Close #1
End Sub

Public Function FirstTocken(ByRef InOutString As String) As String
' this function helps to read space delimited files
Dim localString As String
Dim DelimPos As Long

localString = LTrim(Replace(InOutString, Chr(9), " "))
DelimPos = InStr(localString, " ")
If Not IsNull(DelimPos) Then
    If DelimPos > 1 Then
        FirstTocken = Left(localString, DelimPos - 1)
        InOutString = Right(localString, Len(localString) - DelimPos + 1)
    Else
    FirstTocken = Trim(localString)
    End If
    Else
    FirstTocken = ""
End If
'MsgBox FirstTocken
End Function


Public Sub ShowParametersW(WaveParStruct As WParameters)

With WaveParStruct
frmMain.txtSpPeriod.Text = Format(.WaveSpPeriod, "###0.00")
frmMain.txtOrientation.Text = Format(.Orientation * 180 / cPi, "###0")
frmMain.txtMean.Text = Format(.Mean, "###0.00")
frmMain.txtPhase.Text = Format(.Phase * 180 / cPi, "###0")
frmMain.txtTempPeriodInHz.Text = Format(.WaveTempPeriodInHz, "###0.00")
frmMain.txtTempPeriodInFrames.Text = Format(.WaveTempPeriodInFrames, "###0.00")
frmMain.txtContrast.Text = Format(.Contrast, "###0.00")
frmMain.txtBKGColor.Text = Format(.BKGColor, "###0.00")

frmMain.txtBarW.Text = Format(.BarW, "###0.00")
frmMain.txtRotationPeriodInHz.Text = Format(.RotationPeriodInHz, "#0.0000")
frmMain.txtRotationPeriodInFrames.Text = Format(.RotationPeriodInFrames, "###0.00")
frmMain.txtFlashPeriodInHz.Text = Format(.FlashPeriodInHz, "#0.0000")
frmMain.txtFlashPeriodInFrames.Text = Format(.FlashPeriodInFrames, "###0.00")

frmMain.lblStimPosX.Caption = Format(.StimXposition, "###0.00")
frmMain.lblStimPosY.Caption = Format(.StimYPosition, "###0.00")

End With

End Sub

Public Sub RecalculateTempPeriod()
Dim i As Long
For i = 0 To LoopSequenceSize - 1
    With WaveLoopSequence(i)
        If .WaveTempPeriodInHz <> 0 Then
            .WaveTempPeriodInFrames = CLng(FrameRate / .WaveTempPeriodInHz)
            .WaveTempPeriodInHz = FrameRate / .WaveTempPeriodInFrames
        Else
            .WaveTempPeriodInFrames = 0
            .WaveTempPeriodInHz = 0
        End If
        
        If .RotationPeriodInHz <> 0 Then
            .RotationPeriodInFrames = CLng(FrameRate / .RotationPeriodInHz)
            .RotationPeriodInHz = FrameRate / .RotationPeriodInFrames
        Else
            .RotationPeriodInFrames = 0
            .RotationPeriodInHz = 0
        End If
        
        If .FlashPeriodInHz <> 0 Then
            .FlashPeriodInFrames = CLng(FrameRate / .FlashPeriodInHz)
            .FlashPeriodInHz = FrameRate / .FlashPeriodInFrames
        Else
            .FlashPeriodInFrames = 0
            .FlashPeriodInHz = 0
        End If
    End With
Next i
End Sub


Public Sub WriteLog(XFileName As String)
Dim outFname As String
Dim TextLine As String
'stimtype 0 - M-seq, 1 - wave
On Error GoTo ErrClose:

outFname = App.Path & "\ExpLog.txt"
Open outFname For Append As #1
Print #1, Chr(13) & Chr(10)
Print #1, "Start - " & Now
Print #1, "ini File: " & IniFileName
Print #1, "LUT File: " & LUTFname
Open App.Path & "\" & IniFileName For Input As #3 'copy ini file into log file
Do While Not EOF(3) ' Loop until end of file.
    Line Input #3, TextLine ' Read line into string
    Print #1, TextLine    ' Write to log file
Loop
Close #3

Print #1, "Square wave?: "; frmMain.optWaveType(0).Value
Print #1, "Sin wave?: "; frmMain.optWaveType(1).Value
Print #1, "Flash min-max?: "; frmMain.optFlashType(0).Value
Print #1, "Flash max-mean?: "; frmMain.optFlashType(1).Value
Print #1, "Flash inverse?: "; frmMain.optFlashType(2).Value
Print #1, "Double Stim?: "; frmMain.chkDoubleStim.Value
Print #1, "Use color LUT?: "; frmMain.chkUseLUT.Value
If frmMain.chkUseLUT.Value = 1 Then
  Print #1, "LUT file: " & LUTFname
End If



If Len(XFileName) > 0 Then
    Print #1, "X File: " & XFileName
    Open XFileName For Input As #2 'copy x file into log file
    'Print #1, "Orientation : "; WaveLoopSequence(0).Orientation
    'Print #1, "WaveSpPeriod : "; WaveLoopSequence(0).WaveSpPeriod
    'Print #1, "WaveTempPeriodInHz : "; WaveLoopSequence(0).WaveTempPeriodInHz
    'Print #1, "WaveTempPeriodInFrames : "; WaveLoopSequence(0).WaveTempPeriodInFrames
    'Print #1, "Contrast : "; WaveLoopSequence(0).Contrast
    'Print #1, "Mean : "; WaveLoopSequence(0).Mean
    'Print #1, "Phase : "; WaveLoopSequence(0).Phase
    
    Do While Not EOF(2) ' Loop until end of file.
        Line Input #2, TextLine ' Read line into string
        Print #1, TextLine    ' Write to log file
    Loop
End If
ErrClose:
Close #2
Close #1
End Sub



Public Sub WaitMs(waitTime As Single)
Dim startTime As Currency
Dim stopTime As Currency
'Dim tikPerSec As Currency
'    QueryPerformanceFrequency tikPerSec
    QueryPerformanceCounter startTime

    Do While Not bStop
        DoEvents
        QueryPerformanceCounter stopTime
        If (stopTime - startTime) * 1000# / tikPerSec > waitTime Then
            Exit Do
        End If
    Loop
End Sub


Public Sub WaitMsAfter(waitTime As Double, startTime As Currency)
Dim curTime As Currency
'Dim tikPerSec As Currency
'    QueryPerformanceFrequency tikPerSec
    Do While Not bStop
        DoEvents
        QueryPerformanceCounter curTime
        If (curTime - startTime) * 1000# / tikPerSec > waitTime Then
            Exit Do
        End If
    Loop
End Sub



Public Function ReadCounter(counter As Long) As Long

Call cbCIn32(BoardN, counter, ReadCounter)

'frmMain.Label31.Caption = CStr(ReadCounter)
End Function

' After ResetCounter, counter is not actually reset immediately.
' Reset occurs when the first input arrives.

'may be not for 9513 counter?

Public Sub ResetCounter(counter As Long)
 '   Call cbCLoad32(BoardN, LOADREG1, 0)
    Call cbCLoad32(BoardN, counter, 0)
End Sub

Public Sub SetCounter(counter As Long, co As Integer)
 '   Call cbCLoad32(BoardN, LOADREG1, co)
    Call cbCLoad32(BoardN, counter, co)
End Sub


Public Function CalibratedRed(ByVal co As Double) As Integer
Dim Lutind As Integer
If co > 1 Then
    co = 1
End If

If co < 0 Then
    co = 0
End If

If frmMain.chkUseLUT.Value = 1 Then
    If CurrentLUTStart > LUTLen Then
        CurrentLUTStart = LUTLen
    End If
    If CurrentLUTEnd > LUTLen Then
        CurrentLUTEnd = LUTLen
    End If
    Lutind = CurrentLUTStart - 1 + CInt(co * (CurrentLUTEnd - CurrentLUTStart))
    
    CalibratedRed = CInt(iColorLUT(0, Lutind))
Else
    CalibratedRed = CInt(co * 255)
End If
End Function
Public Function CalibratedGreen(ByVal co As Double) As Integer
Dim Lutind As Integer
If co > 1 Then
    co = 1
End If

If co < 0 Then
    co = 0
End If

If frmMain.chkUseLUT.Value = 1 Then
    If CurrentLUTStart > LUTLen Then
        CurrentLUTStart = LUTLen
    End If
    If CurrentLUTEnd > LUTLen Then
        CurrentLUTEnd = LUTLen
    End If
    Lutind = CurrentLUTStart - 1 + CInt(co * (CurrentLUTEnd - CurrentLUTStart))
    
    CalibratedGreen = CInt(iColorLUT(1, Lutind))
Else
    CalibratedGreen = CInt(co * 255)
End If
End Function
Public Function CalibratedBlue(ByVal co As Double) As Integer
Dim Lutind As Integer
If co > 1 Then
    co = 1
End If

If co < 0 Then
    co = 0
End If

If frmMain.chkUseLUT.Value = 1 Then
    If CurrentLUTStart > LUTLen Then
        CurrentLUTStart = LUTLen
    End If
    If CurrentLUTEnd > LUTLen Then
        CurrentLUTEnd = LUTLen
    End If
    Lutind = CurrentLUTStart - 1 + CInt(co * (CurrentLUTEnd - CurrentLUTStart))
    
    CalibratedBlue = CInt(iColorLUT(2, Lutind))
Else
    CalibratedBlue = CInt(co * 255)
End If
End Function



Public Sub CreatePictureWaveDirect(TrialIndex As Long, dL As Double, Angle As Double, InvInd As Double, maskOnly As Boolean)
Dim i As Long
Dim Trash As Long
Dim HalfCicleNum As Integer
Dim WaveSize As Double 'size in pixels
Dim CenterposX As Double
Dim CenterposY As Double
Dim stimW As Double
Dim stimL As Double
Dim drWidth As Double
Dim zeroN As Long

Dim MaskCenterposX As Double
Dim MaskCenterposY As Double

Dim MaskMaxSize As Double



With WaveLoopSequence(TrialIndex)
    ddsBack.SetForeColor MakeColor(.Mean * (1 + .Contrast * InvInd))
    drWidth = AnglToPixelCoeff * .BarW
    stimW = .StimDiameter * AnglToPixelCoeff
    stimL = .StimLength * AnglToPixelCoeff
    ddsBack.setDrawWidth drWidth
    If Angle < 0 Then
        Angle = Angle + cPi * 2
    End If
    WaveSize = .WaveSpPeriod * AnglToPixelCoeff

    CenterposX = ResolutionX / 2 + dL * Cos(Angle)
    CenterposY = ResolutionY / 2 - dL * Sin(Angle)
    
    MaskCenterposX = ResolutionX / 2 + .StimXposition * AnglToPixelCoeff
    MaskCenterposY = ResolutionY / 2 + .StimYPosition * AnglToPixelCoeff

    Dim dcenX As Double
    Dim dcenY As Double
    
    dcenX = WaveSize * Cos(Angle)
    dcenY = -1 * WaveSize * Sin(Angle)


    
    If frmMain.chkLockPhaseToMask Then
        CenterposX = CenterposX + .StimXposition * AnglToPixelCoeff
        CenterposY = CenterposY + .StimYPosition * AnglToPixelCoeff
    End If
    
    MaskMaxSize = (stimW ^ 2 + stimL ^ 2) ^ 0.5
    
    
    Dim N As Long
    Dim N1 As Long
    Dim projPos As Double
    
    
N = CLng((ResolutionX ^ 2 + ResolutionY ^ 2) ^ 0.5 / (WaveSize * 2)) + 1
N1 = CLng(MaskMaxSize / WaveSize) / 2 + 1
If N1 < N Then
    N = N1
End If
 '       N = 0
If frmMain.chkLockPhaseToMask Then
    zeroN = 0
    projPos = 0
Else
    zeroN = Round((.StimXposition * Cos(Angle) - .StimYPosition * Sin(Angle)) / .WaveSpPeriod)
    projPos = (.StimXposition * Cos(Angle - cPi / 2) - .StimYPosition * Sin(Angle - cPi / 2)) * AnglToPixelCoeff
End If
If maskOnly = False Then
For i = -N + zeroN To N + zeroN
If N1 > N Then
 '       ddsBack.DrawLine CenterposX - ResolutionY * Sin(Angle) + i * dcenX, CenterposY - ResolutionY * Cos(Angle) + i * dcenY, CenterposX + ResolutionY * Sin(Angle) + i * dcenX, CenterposY + ResolutionY * Cos(Angle) + i * dcenY
    ddsBack.DrawLine CenterposX - (-projPos + ResolutionY) * Sin(Angle) + i * dcenX, CenterposY - (-projPos + ResolutionY) * Cos(Angle) + i * dcenY, CenterposX + (projPos + ResolutionY) * Sin(Angle) + i * dcenX, CenterposY + (projPos + ResolutionY) * Cos(Angle) + i * dcenY
Else
    ddsBack.DrawLine CenterposX - (-projPos + MaskMaxSize / 2) * Sin(Angle) + i * dcenX, CenterposY - (-projPos + MaskMaxSize / 2) * Cos(Angle) + i * dcenY, CenterposX + (projPos + MaskMaxSize / 2) * Sin(Angle) + i * dcenX, CenterposY + (projPos + MaskMaxSize / 2) * Cos(Angle) + i * dcenY
End If
 '        ddsBack.DrawLine CenterposX - ((stimW / 2 - .StimXposition * AnglToPixelCoeff)) * Sin(Angle) + i * dcenX, CenterposY - ((stimW / 2 - .StimYPosition * AnglToPixelCoeff)) * Cos(Angle) + i * dcenY, CenterposX + ((stimW / 2 + .StimXposition * AnglToPixelCoeff)) * Sin(Angle) + i * dcenX, CenterposY + ((stimW / 2 + .StimYPosition * AnglToPixelCoeff)) * Cos(Angle) + i * dcenY
 '       ddsBack.DrawLine CenterposX - ((.StimXposition ^ 2 + .StimYPosition ^ 2) ^ 0.5 * AnglToPixelCoeff) * Sin(Angle) + i * dcenX, CenterposY - ((.StimXposition ^ 2 + .StimYPosition ^ 2) ^ 0.5 * AnglToPixelCoeff) * Cos(Angle) + i * dcenY, CenterposX + ResolutionY * Sin(Angle) + i * dcenX, CenterposY + ResolutionY * Cos(Angle) + i * dcenY
 '       ddsBack.DrawLine CenterposX - 10 * Sin(Angle) + i * dcenX, CenterposY - 10 * Cos(Angle) + i * dcenY, CenterposX + 10 * Sin(Angle) + i * dcenX, CenterposY + 10 * Cos(Angle) + i * dcenY
Next i
End If
If N1 <= N Then
ddsBack.SetForeColor MakeColor(.BKGColor)
drWidth = (ResolutionX ^ 2 + ResolutionY ^ 2) ^ 0.5
'drWidth = 5
ddsBack.setDrawWidth drWidth

If frmMain.chkSquareMask.Value = 0 Then
    ddsBack.DrawCircle MaskCenterposX, MaskCenterposY, .StimDiameter * AnglToPixelCoeff / 2 + drWidth / 2
Else
    If frmMain.chkUseStimOri.Value = 0 Then
         ddsBack.DrawLine MaskCenterposX - (stimL + drWidth) / 2, MaskCenterposY - (stimW + drWidth) / 2, MaskCenterposX + (stimL + drWidth) / 2, MaskCenterposY - (stimW + drWidth) / 2
         ddsBack.DrawLine MaskCenterposX - (stimL + drWidth) / 2, MaskCenterposY - (stimW + drWidth) / 2, MaskCenterposX - (stimL + drWidth) / 2, MaskCenterposY + (stimW + drWidth) / 2
         ddsBack.DrawLine MaskCenterposX - (stimL + drWidth) / 2, MaskCenterposY + (stimW + drWidth) / 2, MaskCenterposX + (stimL + drWidth) / 2, MaskCenterposY + (stimW + drWidth) / 2
         ddsBack.DrawLine MaskCenterposX + (stimL + drWidth) / 2, MaskCenterposY - (stimW + drWidth) / 2, MaskCenterposX + (stimL + drWidth) / 2, MaskCenterposY + (stimW + drWidth) / 2
     Else
        ddsBack.DrawLine MaskCenterposX - (stimW + drWidth) / 2 * Sin(Angle) - (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY - (stimW + drWidth) / 2 * Cos(Angle) + (stimL + drWidth) / 2 * Sin(Angle), MaskCenterposX - (stimW + drWidth) / 2 * Sin(Angle) + (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY - (stimW + drWidth) / 2 * Cos(Angle) - (stimL + drWidth) / 2 * Sin(Angle)
        ddsBack.DrawLine MaskCenterposX + (stimW + drWidth) / 2 * Sin(Angle) - (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY + (stimW + drWidth) / 2 * Cos(Angle) + (stimL + drWidth) / 2 * Sin(Angle), MaskCenterposX + (stimW + drWidth) / 2 * Sin(Angle) + (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY + (stimW + drWidth) / 2 * Cos(Angle) - (stimL + drWidth) / 2 * Sin(Angle)
        ddsBack.DrawLine MaskCenterposX - (stimW + drWidth) / 2 * Sin(Angle) - (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY - (stimW + drWidth) / 2 * Cos(Angle) + (stimL + drWidth) / 2 * Sin(Angle), MaskCenterposX + (stimW + drWidth) / 2 * Sin(Angle) - (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY + (stimW + drWidth) / 2 * Cos(Angle) + (stimL + drWidth) / 2 * Sin(Angle)
        ddsBack.DrawLine MaskCenterposX - (stimW + drWidth) / 2 * Sin(Angle) + (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY - (stimW + drWidth) / 2 * Cos(Angle) - (stimL + drWidth) / 2 * Sin(Angle), MaskCenterposX + (stimW + drWidth) / 2 * Sin(Angle) + (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY + (stimW + drWidth) / 2 * Cos(Angle) - (stimL + drWidth) / 2 * Sin(Angle)
    End If
End If
End If
End With
End Sub

Public Sub CreatePictureWaveByBlt(TrialIndex As Long, dL As Double, Angle As Double, InvInd As Double)
Dim Trash As Long
Dim shiftX As Double
Dim shiftY As Double
Dim modshiftX As Double
Dim modshiftY As Double
Dim RectTo As RECT
Dim RectFrom As RECT

Dim stimW As Double
Dim stimL As Double
Dim drWidth As Double
Dim MaskCenterposX As Double
Dim MaskCenterposY As Double


    If Angle < 0 Then
        Angle = Angle + cPi * 2
    End If

    shiftX = -1 * dL * Cos(Angle)
    shiftY = dL * Sin(Angle)
    
'07/21/05 - to improve motion smoothness
If Abs(Cos(Angle)) > Abs(Cos(cPi / 4)) Then
modshiftX = CDbl(CInt(shiftX))
modshiftY = shiftY + (shiftX - modshiftX) * Sin(Angle) / Cos(Angle)
Else
modshiftY = CDbl(CInt(shiftY))
modshiftX = shiftX + (shiftY - modshiftY) * Cos(Angle) / Sin(Angle)
End If

RectTo.Left = 0
RectTo.Top = 0
RectTo.Right = ResolutionX
RectTo.Bottom = ResolutionY

RectFrom.Left = RectTo.Left + XzeroPosition - dL * Cos(Angle)
RectFrom.Top = RectTo.Top + YzeroPosition + dL * Sin(Angle)
RectFrom.Right = RectTo.Right + XzeroPosition - dL * Cos(Angle)
RectFrom.Bottom = RectTo.Bottom + YzeroPosition + dL * Sin(Angle)

'07/21/05
RectFrom.Left = RectTo.Left + XzeroPosition + modshiftX
RectFrom.Top = RectTo.Top + YzeroPosition + modshiftY
RectFrom.Right = RectTo.Right + XzeroPosition + modshiftX
RectFrom.Bottom = RectTo.Bottom + YzeroPosition + modshiftY

'CalculateCompPar

ddsBack.Blt RectTo, ddsBKGImage, RectFrom, DDBLT_WAIT 'SY 04/18/2005

With WaveLoopSequence(TrialIndex)


'    ddsBack.SetForeColor MakeColor(.BKGColor)
'    ddsBack.setDrawWidth 1600
'    If frmMain.chkSquareMask.Value = 0 Then
'        ddsBack.DrawCircle ResolutionX / 2 + .StimXposition * AnglToPixelCoeff, ResolutionY / 2 + .StimYPosition * AnglToPixelCoeff, .StimDiameter * AnglToPixelCoeff / 2 + 800
'    Else
'        ddsBack.DrawLine 0, ResolutionY / 2 + .StimYPosition * AnglToPixelCoeff - (.StimDiameter * AnglToPixelCoeff / 2 + 800), ResolutionX, ResolutionY / 2 + .StimYPosition * AnglToPixelCoeff - (.StimDiameter * AnglToPixelCoeff / 2 + 800)
'        ddsBack.DrawLine ResolutionX / 2 + .StimXposition * AnglToPixelCoeff - (.StimDiameter * AnglToPixelCoeff / 2 + 800), 0, ResolutionX / 2 + .StimXposition * AnglToPixelCoeff - (.StimDiameter * AnglToPixelCoeff / 2 + 800), ResolutionY
'        ddsBack.DrawLine 0, ResolutionY / 2 + .StimYPosition * AnglToPixelCoeff + (.StimDiameter * AnglToPixelCoeff / 2 + 800), ResolutionX, ResolutionY / 2 + .StimYPosition * AnglToPixelCoeff + (.StimDiameter * AnglToPixelCoeff / 2 + 800)
'        ddsBack.DrawLine ResolutionX / 2 + .StimXposition * AnglToPixelCoeff + (.StimDiameter * AnglToPixelCoeff / 2 + 800), 0, ResolutionX / 2 + .StimXposition * AnglToPixelCoeff + (.StimDiameter * AnglToPixelCoeff / 2 + 800), ResolutionY
'    End If
MaskCenterposX = ResolutionX / 2 + .StimXposition * AnglToPixelCoeff
MaskCenterposY = ResolutionY / 2 + .StimYPosition * AnglToPixelCoeff
stimW = .StimDiameter * AnglToPixelCoeff
stimL = .StimLength * AnglToPixelCoeff

    
    
    
ddsBack.SetForeColor MakeColor(.BKGColor)
drWidth = (ResolutionX ^ 2 + ResolutionY ^ 2) ^ 0.5
'drWidth = 5
ddsBack.setDrawWidth drWidth

If frmMain.chkSquareMask.Value = 0 Then
    ddsBack.DrawCircle MaskCenterposX, MaskCenterposY, .StimDiameter * AnglToPixelCoeff / 2 + drWidth / 2
Else
    If frmMain.chkUseStimOri.Value = 0 Then
         ddsBack.DrawLine MaskCenterposX - (stimL + drWidth) / 2, MaskCenterposY - (stimW + drWidth) / 2, MaskCenterposX + (stimL + drWidth) / 2, MaskCenterposY - (stimW + drWidth) / 2
         ddsBack.DrawLine MaskCenterposX - (stimL + drWidth) / 2, MaskCenterposY - (stimW + drWidth) / 2, MaskCenterposX - (stimL + drWidth) / 2, MaskCenterposY + (stimW + drWidth) / 2
         ddsBack.DrawLine MaskCenterposX - (stimL + drWidth) / 2, MaskCenterposY + (stimW + drWidth) / 2, MaskCenterposX + (stimL + drWidth) / 2, MaskCenterposY + (stimW + drWidth) / 2
         ddsBack.DrawLine MaskCenterposX + (stimL + drWidth) / 2, MaskCenterposY - (stimW + drWidth) / 2, MaskCenterposX + (stimL + drWidth) / 2, MaskCenterposY + (stimW + drWidth) / 2
     Else
        ddsBack.DrawLine MaskCenterposX - (stimW + drWidth) / 2 * Sin(Angle) - (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY - (stimW + drWidth) / 2 * Cos(Angle) + (stimL + drWidth) / 2 * Sin(Angle), MaskCenterposX - (stimW + drWidth) / 2 * Sin(Angle) + (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY - (stimW + drWidth) / 2 * Cos(Angle) - (stimL + drWidth) / 2 * Sin(Angle)
        ddsBack.DrawLine MaskCenterposX + (stimW + drWidth) / 2 * Sin(Angle) - (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY + (stimW + drWidth) / 2 * Cos(Angle) + (stimL + drWidth) / 2 * Sin(Angle), MaskCenterposX + (stimW + drWidth) / 2 * Sin(Angle) + (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY + (stimW + drWidth) / 2 * Cos(Angle) - (stimL + drWidth) / 2 * Sin(Angle)
        ddsBack.DrawLine MaskCenterposX - (stimW + drWidth) / 2 * Sin(Angle) - (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY - (stimW + drWidth) / 2 * Cos(Angle) + (stimL + drWidth) / 2 * Sin(Angle), MaskCenterposX + (stimW + drWidth) / 2 * Sin(Angle) - (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY + (stimW + drWidth) / 2 * Cos(Angle) + (stimL + drWidth) / 2 * Sin(Angle)
        ddsBack.DrawLine MaskCenterposX - (stimW + drWidth) / 2 * Sin(Angle) + (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY - (stimW + drWidth) / 2 * Cos(Angle) - (stimL + drWidth) / 2 * Sin(Angle), MaskCenterposX + (stimW + drWidth) / 2 * Sin(Angle) + (stimL + drWidth) / 2 * Cos(Angle), MaskCenterposY + (stimW + drWidth) / 2 * Cos(Angle) - (stimL + drWidth) / 2 * Sin(Angle)
    End If
End If
    
    
    
    
    
End With


End Sub

Public Sub PreparePictureWaveForBlt(TrialIndex As Long, dL As Double, Angle As Double, InvInd As Double)
Dim i As Long
Dim Trash As Long
Dim HalfCicleNum As Integer
Dim WaveSize As Double 'size in pixels
Dim CenterposX As Double
Dim CenterposY As Double

Dim DrawW As Long
Dim DrawStep As Double

'CalculateCompPar

Call FillBKGImageSurf(RGB(255, 0, 0))
With WaveLoopSequence(TrialIndex)
    ddsBKGImage.SetForeColor MakeColor(.Mean * (1 + .Contrast * InvInd))
    DrawW = 3
    ddsBKGImage.setDrawWidth DrawW
    WaveSize = .WaveSpPeriod * AnglToPixelCoeff

    CenterposX = ResolutionX / 2 + XzeroPosition
    CenterposY = ResolutionY / 2 + YzeroPosition
    
    If frmMain.chkLockPhaseToMask Then
        CenterposX = CenterposX + .StimXposition * AnglToPixelCoeff
        CenterposY = CenterposY + .StimYPosition * AnglToPixelCoeff
    End If
    
    Dim dcenX As Double
    Dim dcenY As Double
    DrawStep = 1
    
    dcenX = WaveSize * Cos(.Orientation)
    dcenY = -1 * WaveSize * Sin(.Orientation)
    dcenX = DrawStep * Cos(.Orientation)
    dcenY = -1 * DrawStep * Sin(.Orientation)
    
    Dim N As Long
    N = CLng(ResolutionX * 2 / WaveSize) + 1
    N = CLng(ResolutionX / DrawStep) + 1
     N = CLng((ResolutionX + ResolutionY) / DrawStep) + 1
   ' N = CLng(ResolutionX * 0.7 / WaveSize) + 1
    For i = -N To N
 '       ddsBKGImage.SetForeColor MakeColor(.Mean * (1 + Sin(CDbl(i) * 2 * cPi / 50) * 0.5))
        ddsBKGImage.SetForeColor MakeColor(.Mean * (1 + Sin(i * 2 * cPi * DrawStep / WaveSize) * (.Contrast * InvInd)))
        ddsBKGImage.DrawLine CenterposX - ResolutionY * Sin(.Orientation) + i * dcenX, CenterposY - ResolutionY * Cos(.Orientation) + i * dcenY, CenterposX + ResolutionY * Sin(.Orientation) + i * dcenX, CenterposY + ResolutionY * Cos(.Orientation) + i * dcenY
    Next i
'    ddsBKGImage.SetForeColor MakeColor(.BKGColor)

'    ddsBKGImage.DrawCircle ResolutionX / 2 + XzeroPosition + .StimXposition * AnglToPixelCoeff, ResolutionY / 2 + YzeroPosition + .StimYPosition * AnglToPixelCoeff, 5 * AnglToPixelCoeff / 2 + 400
End With
End Sub


Public Function MakeColor(co As Double)
Dim Lutind As Integer
If co > 1 Then
    co = 1
End If

If co < 0 Then
    co = 0
End If

If frmMain.chkUseLUT.Value = 1 Then
    If CurrentLUTStart > LUTLen Then
        CurrentLUTStart = LUTLen
    End If
    If CurrentLUTEnd > LUTLen Then
        CurrentLUTEnd = LUTLen
    End If
    Lutind = CurrentLUTStart - 1 + CInt(co * (CurrentLUTEnd - CurrentLUTStart))
    
    MakeColor = RGB(CInt(iColorLUT(0, Lutind)), CInt(iColorLUT(1, Lutind)), CInt(iColorLUT(2, Lutind)))
Else
    MakeColor = RGB(CInt(co * 255), CInt(co * 255), CInt(co * 255))
End If

End Function


Public Sub CalculateParamByFrame(TrInd As Long, Frame As Long, ByRef dL As Double, ByRef Angle As Double, ByRef BarIsOn As Boolean, ByRef FlashPhase As Double)
Dim ddL As Double
Dim dAngle As Double

With WaveLoopSequence(TrInd)
    ddL = .WaveSpPeriod * AnglToPixelCoeff * .WaveTempPeriodInHz / FrameRate
    dAngle = 2 * cPi * .RotationPeriodInHz / FrameRate
    If frmMain.optWaveType(1).Value = True Then
        dAngle = 0 'no rotation for sin wave
    End If


    dL = .WaveSpPeriod * AnglToPixelCoeff * .Phase / (2 * cPi) + ddL * Frame
    Angle = .Orientation + dAngle * Frame
      
    'Do While Abs(dL) > .WaveSpPeriod * AnglToPixelCoeff
    Do While Abs(dL) > .WaveSpPeriod * AnglToPixelCoeff / 2
        dL = dL - Sgn(dL) * .WaveSpPeriod * AnglToPixelCoeff
    Loop

    If .FlashPeriodInFrames > 0 Then
        If Frame Mod .FlashPeriodInFrames < 0.5 * .FlashPeriodInFrames Then
            BarIsOn = True
        Else
            BarIsOn = False
        End If
        FlashPhase = (Frame * 360 / .FlashPeriodInFrames) Mod 360
    Else
        BarIsOn = True
        FlashPhase = 0
    End If
End With

End Sub

Public Function UintToLng(inUint As Integer) As Long
    If inUint < 0 Then
        UintToLng = inUint + 65536
    Else
        UintToLng = inUint
    End If
End Function

Public Function LngToUint(inLng As Long) As Integer
    If inLng > 32767 Then
        LngToUint = inLng - 65536
    Else
        LngToUint = inLng
    End If
End Function

Public Sub ChangeRamp(center As Double, chCoeff As Double)

Dim NewR As DDGAMMARAMP
Dim i As Integer
Dim j As Integer
Dim Lutind As Integer
Dim Lutind2 As Integer

If frmMain.chkUseLUT.Value = 1 Then
    If CurrentLUTStart > LUTLen Then
        CurrentLUTStart = LUTLen
    End If
    If CurrentLUTEnd > LUTLen Then
        CurrentLUTEnd = LUTLen
   End If
    For i = CurrentLUTStart To CurrentLUTEnd

        Lutind = i - 1
        Lutind2 = (i - 1 - (CurrentLUTStart - 1 + (CurrentLUTEnd - CurrentLUTStart) * center)) * chCoeff + (CurrentLUTStart - 1 + (CurrentLUTEnd - CurrentLUTStart) * center)
        NewR.red(CInt(iColorLUT(0, Lutind))) = InitialRamp.red(CInt(iColorLUT(0, Lutind2)))
        NewR.green(CInt(iColorLUT(1, Lutind))) = InitialRamp.green(CInt(iColorLUT(1, Lutind2)))
        NewR.blue(CInt(iColorLUT(2, Lutind))) = InitialRamp.blue(CInt(iColorLUT(2, Lutind2)))
    Next i
Else
    For i = 0 To 255
        Lutind = i
        Lutind2 = (i - 255 * center) * chCoeff + 255 * center
        
        NewR.red(Lutind) = InitialRamp.red(Lutind2)
        NewR.green(Lutind) = InitialRamp.green(Lutind2)
        NewR.blue(Lutind) = InitialRamp.blue(Lutind2)
    Next i
End If


GCon.SetGammaRamp DDSGR_DEFAULT, NewR

End Sub
Public Sub ShowAvg()
Dim maxCount As Double
Dim i As Integer

maxCount = 4
For i = 0 To loopsizeToAVG * 2 - 1
    If AvgCounts(i) > maxCount Then
        maxCount = AvgCounts(i)
    End If
Next i
frmMain.pctCounts.Cls

frmMain.pctCounts.ScaleTop = -1
frmMain.pctCounts.ScaleHeight = CInt(maxCount) + 2
frmMain.pctCounts.ScaleLeft = 0
frmMain.pctCounts.ScaleWidth = loopsizeToAVG * 2 + 1
frmMain.pctCounts.Line (0, CInt(maxCount))-(loopsizeToAVG * 2 + 2, CInt(maxCount)), RGB(0, 255, 0)

For i = 1 To CInt(maxCount)
    frmMain.pctCounts.Line (0, CInt(maxCount) - i)-(loopsizeToAVG * 2 + 2, CInt(maxCount) - i), RGB(0, 100, 0)
Next i


For i = 0 To loopsizeToAVG * 2 - 2
    frmMain.pctCounts.Line (i + 1, -1)-(i + 1, CInt(maxCount) + 2), RGB(0, 100, 0)
    frmMain.pctCounts.Line (i + 1, CInt(maxCount) - AvgCounts(i))-(i + 2, CInt(maxCount) - AvgCounts(i + 1)), RGB(255, 0, 0)
Next i
    frmMain.pctCounts.Line (loopsizeToAVG * 2, -1)-(loopsizeToAVG * 2, CInt(maxCount) + 2), RGB(0, 100, 0)


End Sub

Public Sub SaveAvg(FileName)
Dim i As Integer
Open FileName For Output As #1
For i = 0 To loopsizeToAVG * 2 - 1
Print #1, AvgCounts(i)
Next i
Close #1
End Sub
