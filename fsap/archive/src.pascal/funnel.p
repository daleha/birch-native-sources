(*$DEBUG- *)
(* ********************************************************* *)
(*                                                           *)
(* FUNNEL  Version  5/13/91   Standard Pascal                *)
(*         Brian Fristensky                                  *)
(*         Dept. of Plant Science                            *)
(*         University of Manitoba                            *)
(*         Winnipeg, MB R3T 2N2  CANADA                      *)
(*                                                           *)
(* Copyright (c) 1987, 1990 by Brian Fristensky              *)
(* !!! in comment indicates feature which may need change    *)
(* ********************************************************* *)
program FUNNEL(input,output (*,INFILE,OUTFILE*));
(*!!! Some Pascals require file parameters in program heading *)
 
const VERSION = 'FUNNEL            Version  5/13/91';
      MAXWORD = 25;
      MAXLINE = 10000;
(* BEGIN MODULE INTLIMITS *)
(*!!!  MAXINT =  2147483647; *)
(*!!!  MININT = -2147483647; *)
(* END MODULE INTLIMITS         VERSION= 'SUNMODS     Version  9/12/91'; *)
 
type
(* BEGIN MODULE TYPE.WORD *)
       (*   <word>::= <non-blank char>[<non-blank char>] *)
       WORD    = record
                 LEN:integer;
                 STR:array[1..MAXWORD] of char
                 end;
(* END MODULE TYPE.WORD         VERSION= 'SUNMODS     Version  9/12/91'; *)

(* BEGIN MODULE TYPE.LINE *)
       CHARARRAY = packed array[1..MAXLINE] of char;
       LINE = record
                STR:CHARARRAY;
                LEN:integer
                end;
(* END MODULE TYPE.LINE         VERSION= 'SUNMODS     Version  9/12/91'; *)
 
var INFILE,OUTFILE:text;
    IFN,OFN:LINE;
    CH: char;
    I,LINESIZE: integer;
    COMMENT     : boolean;

(* BEGIN MODULE INPLINE *)
    (*  Read a line from the console.    *)
    procedure INPLINE(var W:LINE);
      var I : integer;
         CH,BACKSPACE : char;
      begin
        BACKSPACE:=  chr(8);
(*!!!   if eoln(input) then readln; *)
(*!!!   get(input); *)
        with W do begin
          LEN:= 0;
          while not eoln(input) do begin
            read(CH);
(*!!!       if not eoln then*)
              if (CH = BACKSPACE) and (LEN > 0) then begin
                write(BACKSPACE); LEN:= LEN-1 end
              else if LEN < MAXLINE then begin
                LEN := LEN + 1;
                STR[LEN]:= CH
                end
            end;
          readln(input);
          for I := LEN+1 to MAXLINE do STR[I] := ' '
          end
      end; (* INPLINE *)
(* END MODULE INPLINE         VERSION= 'SUNMODS     Version  9/12/91'; *)
 
(* BEGIN MODULE GETFILE *)
(* The c function access is found in sys/file.h. access must be declared in
   the outermost scope of the Pascal program. See access (2) manual pages *)
(*!!!*) function access(PATH:CHARARRAY; MODE:integer):integer; external c;

(************************************************************)
(* Open  files, checking for I/O errors.                    *)
(* Assumes user has already been prompted for filename.     *)
(* Syntax: GETFILE(F,'I',<filename>);                       *)
(*              opens file F for input                      *)
(*         GETFILE(F,'O',<filename>);                       *)
(*              opens file F for output                     *)
(* !!! NOTE: this non-Standard procedure must be rewritten  *)
(* for each different version of Pascal.                    *)
(************************************************************)
procedure GETFILE(var F:text; FTYPE:char; var FILENAME:LINE);

  var CH:char;
      OKAY:boolean;

  (* - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
  (* =true if file exists, =false if it doesn't *)
  (* A comparable procedure is provided with ATT System V Pascal *)
  (* A special thank you to Mark Foster, CIS Rearch Computing, Univ. of Penn.
     and Glenn Gribble of Synaptics Inc. for showing me how to use access()
     within Pascal programs. B.F. *)
  function EXISTS(FILENAME:LINE):boolean;
    const F_OK=0; (* checks if file exists *)
    begin (* EXISTS *)
      with FILENAME do begin
        STR[LEN+1]:=chr(0); (* c strings terminate with null *)
        if access(STR,F_OK)=-1 then EXISTS:=false
        else  EXISTS:=true
        end (* with FILENAME *)
    end; (* EXISTS *)

  begin (* GETFILE - - - - - - - - - - - - - - - - - - - - - *)
      (* Read in filename.*)
      INPLINE(FILENAME);
      case FTYPE of
        (* Input file must exist before being opened.*)
        'I':begin (* input file *)
              while not(EXISTS(FILENAME)) do begin
                writeln('***FILE NOT FOUND');
                writeln('Enter another filename:');
                INPLINE(FILENAME);
                end; (* while *) 
(*!!!*)       reset(F,FILENAME.STR)
            end; (* I *)

        (* If output file already exists, ask user if he really wants
           to over write it. *)
        'O':begin
              repeat
                OKAY:=true;
                if EXISTS(FILENAME) then begin
                  repeat  
                    writeln('*** WARNING! File already exists. Overwrite?[Y|N]');
                    readln(CH)
                  until CH in ['Y','y','N','n'];
                  if not (CH in ['Y','y']) then begin
                    writeln('Enter another filename:');
                    INPLINE(FILENAME);
                    OKAY:=false
                    end
                  end (* EXISTS(FILENAME) *)
              until OKAY;
            rewrite(F,FILENAME.STR)  (* output file *)
            end (* O *)
        end (*case*)
  end;  (* GETFILE *)
(* END MODULE GETFILE         VERSION= 'SUNMODS     Version  9/12/91'; *)
 
(* BEGIN MODULE INPWORD *)
  (* Read a WORD from the terminal. *)
   procedure INPWORD(var W:WORD);
      var I : integer;
         CH : char;
         BACKSPACE: char;
      begin
        with W do begin
          BACKSPACE:= chr(8);
          LEN:= 0;
          while input^ = ' ' do read(input,CH);
          while (input^ <> ' ') and (not eoln(input)) do begin
            read(CH);
(*!!!       if not eoln then*)
              if (CH= BACKSPACE) and (LEN > 0) then begin
                write(BACKSPACE);LEN:= LEN-1 end
              else if LEN < MAXWORD then begin
                LEN := LEN + 1;
                STR[LEN]:= CH
                end
            end;
        (*readln;*)
          for I := LEN+1 to MAXWORD do STR[I] := ' '
          end
      end; (* INPWORD *)
(* END MODULE INPWORD         VERSION= 'SUNMODS     Version  9/12/91'; *)
 
(* BEGIN MODULE GETINTEGER *)
  (* Prompts user for an integer, checks whether all characters are digits,*)
  (* and whether number is within desired range; harasses user until valid *)
  (* integer is received. *)
  procedure GETINTEGER(var NUM:integer;LBOUND,HBOUND:integer);
 
   var I,VAL,ORDZERO  : integer;
     LEGAL,INRANGE,NEGATIVE      : boolean;
     CH : char;
     NUMWORD : WORD;
 
  begin
    INRANGE:= false;
    LEGAL := false;
    NEGATIVE := false;
    ORDZERO:= ord('0');
    repeat
      repeat
        INPWORD(NUMWORD);
        with NUMWORD do begin
 
          (* Evaluate sign, if any *)
          I:=1;CH:=STR[I];
          if CH = '-' then begin
             NEGATIVE := true;
             I:=I+1; CH:= STR[I]
             end
          else begin NEGATIVE:= false; if CH='+' then begin
                         I:=I+1; CH:= STR[I] end end;
 
          (* Evaluate unsigned integer *)
          NUM:= 0;
          while CH in ['0'..'9'] do
           begin
             VAL:= ord(CH)-ORDZERO;
             NUM:= (NUM * 10) + VAL;
             I:= I+1; CH:= STR[I]
             end;
           if I > LEN then LEGAL:= true
           else  begin
              LEGAL:= false;
              writeln;
              writeln('Illegal character encountered.');
              writeln('Enter new value:  ');
             end
         end (* with *)
       until LEGAL;
 
       (* If the number entered was negative, multiply *)
       (* NUM by -1. Check range of number.            *)
       if NEGATIVE then NUM:= - NUM;
       if (NUM >= LBOUND) and (NUM <= HBOUND)
         then INRANGE := true
       else begin
         INRANGE := false;
         writeln;
         writeln('Number is out of range.');
         writeln('Please enter new value:')
         end;
     until INRANGE;
   end;   (* GETINTEGER *)
(* END MODULE GETINTEGER         VERSION= 'SUNMODS     Version  9/12/91'; *)
 
begin
(* BEGIN MODULE STARTUP *)
(* Peform operations which must be done at beginning of main
   procedure. *)
(*!!!   TERMIN(input);    Open input for interactive use *)
(*!!!   TERMOUT(output);   "   output "      "        "  *)
     writeln(VERSION:50);
     writeln;
(* END MODULE STARTUP         VERSION= 'SUNMODS     Version  9/12/91'; *)
(*Ask user for INFILE and OUTFILE filenames*)
  writeln('INPUT FILENAME?');
  GETFILE(INFILE,'I',IFN);
  writeln('OUTPUT FILENAME?');
  GETFILE(OUTFILE,'O',OFN);
  writeln('Enter LINESIZE:');
  GETINTEGER(LINESIZE,1,MAXLINE);readln;
  (*READ INPUT FILE CHARACTER BY CHARACTER*)
  (*AND WRITE EACH CHARACTER INTO OUTPUT FILE*)
  I:= 0;COMMENT:=false;
  while not eof(INFILE) do begin
      if COMMENT then I:= 0;
      COMMENT:= false;
      while not eoln(INFILE) do begin
          read(INFILE,CH);
          if CH in [';','>'] then begin
                           COMMENT := true;
                           if I>0 then writeln(OUTFILE) 
                           end;
          case COMMENT of
            true: write(OUTFILE,CH);
            false:if CH in ['A'..'Z','a'..'z','*']
             then begin
                    write(OUTFILE,CH);
                    I:= I + 1;
                    if I = LINESIZE then  begin
                      I:=0; writeln(OUTFILE) end
                    end
              end
          end;
      readln(INFILE);
      if COMMENT then writeln(OUTFILE) 
      end;
 (*!!!*)CLOSE(OUTFILE) 
  end. (* FUNNEL *)
 
