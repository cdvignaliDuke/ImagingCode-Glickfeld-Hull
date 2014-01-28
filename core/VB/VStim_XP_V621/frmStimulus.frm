VERSION 5.00
Begin VB.Form frmStimulus 
   BackColor       =   &H80000008&
   BorderStyle     =   0  'None
   ClientHeight    =   5700
   ClientLeft      =   7215
   ClientTop       =   5295
   ClientWidth     =   7860
   KeyPreview      =   -1  'True
   LinkTopic       =   "Form1"
   PaletteMode     =   1  'UseZOrder
   ScaleHeight     =   380
   ScaleMode       =   3  'Pixel
   ScaleWidth      =   524
   ShowInTaskbar   =   0   'False
End
Attribute VB_Name = "frmStimulus"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Private Sub Form_KeyDown(KeyCode As Integer, Shift As Integer)
Select Case KeyCode
Case vbKeyEscape
    bStop = True
Case vbKeyUp
    StimYPosition = StimYPosition + 0.1
    WaveDefault.StimYPosition = StimYPosition
Case vbKeyDown
    StimYPosition = StimYPosition - 0.1
    WaveDefault.StimYPosition = StimYPosition
Case vbKeyLeft
    StimXposition = StimXposition - 0.1
    WaveDefault.StimXposition = StimXposition


Case vbKeyRight
    StimXposition = StimXposition + 0.1
    WaveDefault.StimXposition = StimXposition

Case vbKeyN

End Select

End Sub

Private Sub Form_Resize()
If BackSurfIsReady Then
    CheckSurfAndRestore
End If
End Sub
