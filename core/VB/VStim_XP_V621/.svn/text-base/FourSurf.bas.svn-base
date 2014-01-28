Attribute VB_Name = "FourSurf"
Option Explicit
'module provides basic functionality of
'two full-sreen flipable surfaces and
'two additional surfaces for image manipulation

Dim dx As New DirectX7
Public dd As DirectDraw7

Public ddsFront As DirectDrawSurface7
Public ddsBack As DirectDrawSurface7
Public ddsBKGImage As DirectDrawSurface7
Public ddsFRGImage As DirectDrawSurface7

Public XzeroPosition As Long 'coordinate of pixel of ddsBKGImage and ddsFRGImage
Public YzeroPosition As Long ' ... which correspondes to initial 0 coordinate of ddsBack
Public ImgSurfRect As RECT

Public BackSurfIsReady As Boolean
Dim HaveSecondarySurf As Boolean

Dim ddsdFront As DDSURFACEDESC2
Dim ddsdImage As DDSURFACEDESC2

Dim ddClipper As DirectDrawClipper

Dim DDCAPS As DDSCAPS2
Dim winDC As Long
Dim winhDC As Long

Dim guiMon As String
Dim stimMon As String

Public Type MONITORINFO
    cbSize As Long
    rcMonitor As RECT
    rcWork As RECT
    dwFlags As Long
End Type


Public needsRestore As Boolean 'flag shows if we are in the loop waiting to restore surfaces
                                '... see CheckSurfAndRestore

Declare Function RestoreDC Lib "gdi32" (ByVal hdc As Long, ByVal nSavedDC As Long) As Long
Declare Function SaveDC Lib "gdi32" (ByVal hdc As Long) As Long
Declare Function GetDC Lib "user32" (ByVal hwnd As Long) As Long
Declare Function GetMonitorInfo Lib "user32" Alias "GetMonitorInfoA" (ByVal hMonitor As Long, lpmi As MONITORINFO) As Long

Public GCon As DirectDrawGammaControl
Public InitialRamp As DDGAMMARAMP


'call PositionForms before SetDispMode !!!
'function is for version compatability - use PositionForm
Public Sub PositionForms(guiForm As Form, guiMonIndex As Long, stimForm As Form, stimMonIndex As Long)
Dim sDX As DirectDrawEnum

Dim hmGui As Long
Dim hmStim As Long

Set sDX = dx.GetDDEnum
guiMon = sDX.GetGuid(guiMonIndex)
stimMon = sDX.GetGuid(stimMonIndex)

hmGui = sDX.GetMonitorHandle(guiMonIndex)
hmStim = sDX.GetMonitorHandle(stimMonIndex)

Dim mi As MONITORINFO
mi.cbSize = Len(mi)

If hmGui <> 0 Then
    GetMonitorInfo hmGui, mi
    guiForm.Left = (mi.rcMonitor.Left + 10) * Screen.TwipsPerPixelX
    guiForm.Top = (mi.rcMonitor.Top + 10) * Screen.TwipsPerPixelY
End If

If hmStim <> 0 Then
    GetMonitorInfo hmStim, mi
    stimForm.Left = (mi.rcMonitor.Left + 10) * Screen.TwipsPerPixelX
    stimForm.Top = (mi.rcMonitor.Top + 10) * Screen.TwipsPerPixelY
End If


End Sub

'call PositionForm before SetDispMode !!!
'function positions form on specified coordinate of specified monitor
' !!! - important - index of first monitor =2 and so on
Public Sub PositionForm(frm As Form, monIndex As Long, x As Long, y As Long)

Dim sDX As DirectDrawEnum
Dim monH As Long

Set sDX = dx.GetDDEnum
monH = sDX.GetMonitorHandle(monIndex)

Dim mi As MONITORINFO
mi.cbSize = Len(mi)

If monH <> 0 Then
    GetMonitorInfo monH, mi
    frm.Left = (mi.rcMonitor.Left + x) * Screen.TwipsPerPixelX
    frm.Top = (mi.rcMonitor.Top + y) * Screen.TwipsPerPixelY
End If

End Sub

Public Sub SetDispMode()
'old function - use with PositionForms
'new version - SetFullScreenDispMode
winhDC = GetDC(frmStimulus.hwnd)
winDC = SaveDC(winhDC)
Set dd = dx.DirectDrawCreate(stimMon)
Call dd.SetCooperativeLevel(frmStimulus.hwnd, DDSCL_EXCLUSIVE Or DDSCL_FULLSCREEN)
    
    
'(SY)dd.SetDisplayMode ResolutionX, ResolutionY, ResolutionC, 0, DDSDM_DEFAULT

End Sub

Public Sub SetFullScreenDispMode(frm As Form, monInd As Long)

Dim sDX As DirectDrawEnum
Dim stimMon As String

winhDC = GetDC(frm.hwnd) 'will be used in RestoreDispMode
winDC = SaveDC(winhDC)
Set sDX = dx.GetDDEnum
stimMon = sDX.GetGuid(monInd)

Set dd = dx.DirectDrawCreate(stimMon)
Call dd.SetCooperativeLevel(frm.hwnd, DDSCL_EXCLUSIVE Or DDSCL_FULLSCREEN)


Dim hwCaps As DDCAPS
Dim helCaps As DDCAPS
dd.GetCaps hwCaps, helCaps
If (hwCaps.lCaps2 And DDCAPS2_PRIMARYGAMMA) = 0 Then
    MsgBox "Gamma not supported"
End If


End Sub

Public Sub RestoreDispMode()

Dim Trash As Long

dd.RestoreDisplayMode
dd.SetCooperativeLevel 0, DDSCL_NORMAL
Set dd = Nothing
Trash = RestoreDC(winhDC, winDC)

End Sub


Public Sub CreateMainSurfaces()
'old version - now use CreateMainSurfacesWithClipper
'creates two flipable surfaces - ddsFront and ddsBack
    
'ddsdFront.dwSize = Len(ddsdFront)

ddsdFront.lFlags = DDSD_CAPS Or DDSD_BACKBUFFERCOUNT
ddsdFront.ddsCaps.lCaps = DDSCAPS_PRIMARYSURFACE Or DDSCAPS_FLIP Or DDSCAPS_COMPLEX Or DDSCAPS_VIDEOMEMORY
ddsdFront.lBackBufferCount = 1
Set ddsFront = dd.CreateSurface(ddsdFront)
DDCAPS.lCaps = DDSCAPS_BACKBUFFER
Set ddsBack = ddsFront.GetAttachedSurface(DDCAPS)
   
'create clipper for ddsBack for easy blitting to that surface
Set ddClipper = dd.CreateClipper(0)
ddClipper.SetHWnd frmStimulus.hwnd
ddsBack.SetClipper ddClipper

BackSurfIsReady = True



End Sub

Public Sub CreateMainSurfacesWithClipper(ClipperForm As Form)
'old version - noe use CreateMainSurfacesWithClipper
'creates two flipable surfaces - ddsFront and ddsBack
    
'ddsdFront.dwSize = Len(ddsdFront)

ddsdFront.lFlags = DDSD_CAPS Or DDSD_BACKBUFFERCOUNT
ddsdFront.ddsCaps.lCaps = DDSCAPS_PRIMARYSURFACE Or DDSCAPS_FLIP Or DDSCAPS_COMPLEX Or DDSCAPS_VIDEOMEMORY
ddsdFront.lBackBufferCount = 1
Set ddsFront = dd.CreateSurface(ddsdFront)
DDCAPS.lCaps = DDSCAPS_BACKBUFFER
Set ddsBack = ddsFront.GetAttachedSurface(DDCAPS)

Call ddsFront.GetSurfaceDesc(ddsdFront)
   
'create clipper for ddsBack for easy blitting to that surface
Set ddClipper = dd.CreateClipper(0)
ddClipper.SetHWnd ClipperForm.hwnd
ddsBack.SetClipper ddClipper

BackSurfIsReady = True

Set GCon = ddsFront.GetDirectDrawGammaControl

GCon.GetGammaRamp DDSGR_DEFAULT, InitialRamp

End Sub

Public Sub CreateSecondarySurfaces(xScale As Single, yScale As Single)

' create two surfaces for drawing of BKG and FRG images
    Set ddsBKGImage = Nothing
    Set ddsFRGImage = Nothing
   
    ddsdImage.lFlags = DDSD_CAPS Or DDSD_HEIGHT Or DDSD_WIDTH
    ddsdImage.ddsCaps.lCaps = DDSCAPS_OFFSCREENPLAIN Or DDSCAPS_VIDEOMEMORY
    
    ImgSurfRect.Top = 0
    ImgSurfRect.Left = 0
    ImgSurfRect.Right = ddsdFront.lWidth * xScale
    ImgSurfRect.Bottom = ddsdFront.lHeight * yScale

    ddsdImage.lHeight = ImgSurfRect.Bottom - ImgSurfRect.Top
    ddsdImage.lWidth = ImgSurfRect.Right - ImgSurfRect.Left
    
    XzeroPosition = (ddsdImage.lWidth - ddsdFront.lWidth) / 2
    YzeroPosition = (ddsdImage.lHeight - ddsdFront.lHeight) / 2

    Set ddsBKGImage = dd.CreateSurface(ddsdImage)
    Set ddsFRGImage = dd.CreateSurface(ddsdImage)
    
    'set transparant color for ddsFRGImage to be black
    'anything in black on ddsFRGImage will be NOT copied during blit FROM ddsFRGImage
    Dim key As DDCOLORKEY
    key.high = 0
    key.low = 0
    ddsFRGImage.SetColorKey DDCKEY_SRCBLT, key

    HaveSecondarySurf = True
End Sub


Public Sub DestroySurface()
    Set ddsBKGImage = Nothing
    Set ddsFRGImage = Nothing
    Set ddsBack = Nothing
    Set ddsFront = Nothing
    BackSurfIsReady = False
End Sub


Public Sub FillTwoSurf(fillColor As Long)
    FillBkgSurf (fillColor)
    Flip
    WaitFrameStart
    DoEvents
    FillBkgSurf (fillColor)
End Sub

Public Sub FillBkgSurf(fillColor As Long)
'fillColor can be created by RGB(r,g,b)- actually use of color masks may be needed
Dim ree As RECT
Dim oldColor As Long
    CheckSurfAndRestore

    oldColor = ddsBack.GetForeColor
    ddsBack.SetFillColor (fillColor)
    ddsBack.SetFillStyle 0
    ddsBack.SetForeColor (fillColor)
    ddsBack.DrawBox 0, 0, ResolutionX, ResolutionY
    ddsBack.SetForeColor oldColor
    ddsBack.SetFillStyle 1
'    ddsBack.BltColorFill ree, fillColor
End Sub


Public Sub FillBKGImageSurf(fillColor As Long)
'fillColor can be created by RGB(r,g,b)- actually use of color masks may be needed
Dim oldColor As Long
    CheckSurfAndRestore

    oldColor = ddsBKGImage.GetForeColor
    ddsBKGImage.SetFillColor (fillColor)
    ddsBKGImage.SetFillStyle 0
    ddsBKGImage.SetForeColor (fillColor)
    ddsBKGImage.DrawBox ImgSurfRect.Left, ImgSurfRect.Top, ImgSurfRect.Right, ImgSurfRect.Bottom
    ddsBKGImage.SetForeColor oldColor
    ddsBKGImage.SetFillStyle 1
        
End Sub

Public Sub FillFRGImageSurf(fillColor As Long)
'fillColor can be created by RGB(r,g,b)- actually use of color masks may be needed
Dim oldColor As Long
    CheckSurfAndRestore
    
    oldColor = ddsFRGImage.GetForeColor
    ddsFRGImage.SetFillColor (fillColor)
    ddsFRGImage.SetFillStyle 0
    ddsFRGImage.SetForeColor (fillColor)
    ddsFRGImage.DrawBox ImgSurfRect.Left, ImgSurfRect.Top, ImgSurfRect.Right, ImgSurfRect.Bottom
    ddsFRGImage.SetForeColor oldColor
    ddsFRGImage.SetFillStyle 1
    
End Sub

Public Sub Flip()
Err.Clear
    Do
    CheckSurfAndRestore
    
    ddsFront.Flip Nothing, DDFLIP_WAIT
    If Err = DDERR_SURFACELOST Then
        ddsFront.restore
    End If
    Loop Until Err = 0

End Sub
Public Sub FlipNoVSync()
Err.Clear
    Do
    CheckSurfAndRestore
    
    ddsFront.Flip Nothing, DDFLIP_NOVSYNC
    If Err = DDERR_SURFACELOST Then
        ddsFront.restore
    End If
    Loop Until Err = 0

End Sub
Public Sub WaitFrameStart()
    dd.WaitForVerticalBlank DDWAITVB_BLOCKEND, 0
    
End Sub


Public Sub WaitFrameEnd()
    dd.WaitForVerticalBlank DDWAITVB_BLOCKBEGIN, 0
End Sub


Public Sub CheckSurfAndRestore()
'Dim needsRestore As Boolean

needsRestore = False

Do Until dd.TestCooperativeLevel = DD_OK
    DoEvents    'we are sitting here if we are waiting for stim form to be maximized
    needsRestore = True
Loop

If needsRestore Then
    dd.RestoreAllSurfaces
    If HaveSecondarySurf Then
        CreateSecondarySurfaces CSng(ddsdImage.lWidth) / CSng(ddsdFront.lWidth), CSng(ddsdImage.lHeight) / CSng(ddsdFront.lHeight) ' if have secondary surf
    End If
    needsRestore = False
End If

End Sub
