%emulation of Winodows API function GetProfileString
% retString=GetProfileString(filename,sectionname,paramname)


function retString=GetProfileString(filename,sectionname,paramname)


fid=fopen(filename);
retString=[];
find_flag=0;

sec_name_string=['[' sectionname ']'];
par_name_string=[paramname '='];

while feof(fid)==0
   line=fgetl(fid);
   pos=findstr(line,sec_name_string);
   if ~isempty(pos)
      if pos(1)==1
         break;
      end
   end
end

while feof(fid)==0
   line=fgetl(fid);
   pos=findstr(line,par_name_string);
   if ~isempty(pos)
      if pos(1)==1
         find_flag=1;
         break;
      end
   end
   pos=findstr(line,'[');
   if ~isempty(pos)
      if pos(1)==1
         find_flag=0;
         break;
      end
   end 
end

if find_flag==0
   string=[];
else
   [token,rem] = strtok(line,'=');
   retString=deblank(rem(2:end));
end

fclose(fid);