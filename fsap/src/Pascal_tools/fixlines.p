(* PIR tapes come in fixed records of 80 characters, with NO RETURN CHARACTER
   at the end of the record.  Many tape utilties (eg. UORCOPY on IBM mainframes)
   will correctly reappend such characters, and maybe even strip off trailing
   blanks.  That is what this program does. *)
program FIXLINES(input,output);

  var LASTNONBLANK,
      COLUMN,
      START,LENGTH:integer;
      LINE:array[1..80] of char;

  begin
    while not eof do begin
       (* Read in an 80 character line. *)
       COLUMN:=0;
       while (COLUMN<80) and (not eof) do begin
             COLUMN:=COLUMN+1;
             read(LINE[COLUMN])
             end;
       LENGTH:=COLUMN;

       (* Write the line to output, omitting leading or trailing ## markers 
	  if it is a continuation line.  If it is not a continuation, write
          a newline character.*)
       if (LINE[1]='#') and (LINE[2]='#') then START:=3 else START:=1;

       if (LINE[79]='#') and (LINE[80]='#') and (LENGTH=80) then 
          for COLUMN:= START to 78 do write(LINE[COLUMN])
       else begin
          COLUMN:=LENGTH;
          while (LINE[COLUMN]=' ') and (COLUMN > START) do COLUMN:=COLUMN-1;
          LASTNONBLANK:=COLUMN;
          for COLUMN:= START to LASTNONBLANK do write(LINE[COLUMN]);
          writeln
          end
       end (* not eof *)
  end. (* FIXLINES *)
