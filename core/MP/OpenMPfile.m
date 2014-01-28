function mpfile=OpenMPfile (fname)
mpfile=actxcontrol('MPFile.Data');
openResult = invoke(mpfile, 'OpenMPFile', fname);

