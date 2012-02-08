(* Revised 10/5/90 by Brian Fristensky*)
(* SUN       Pascal *)
(* !!! in comment indicates feature which may need change *)
program MODULES(RAW,MODLIB,SOUT,MODCAT,LIST,output);
(*!!! Some Pascals may require filenames in program heading *)
(* MODULE REPLACEMENT PROGRAM
   THOMAS SCHNEIDER COPYRIGHT (C) 1982
 
DOCUMENTATION AND OPERATION OF THIS PROGRAM ARE DEFINED IN MODDEF 2.09.
 
MODULE LIBRARIES REQUIRED: DELMAN, DELMODS
*)
 
label 1; (* THE END OF THE PROGRAM *)
 
const
(* BEGIN MODULE VERSION *)
VERSION = 'MODULE 5.34  82 JULY 28 TDS  (revised 10/5/90 by BF )';
(* END MODULE VERSION *)
(* BEGIN MODULE DELMAN.DESCRIBE.MODULE *)
(*
NAME
      MODULE: MODULE REPLACEMENT PROGRAM
 
SYNOPSIS
      MODULE(RAW: IN, MODLIB: IN, SOUT: OUT,
             MODCAT: INOUT, LIST: OUT, output: OUT)
 
FILES
      RAW: THE SOURCE PROGRAM OR FILE
      MODLIB: A LIBRARY OF MODULES (IF EMPTY, MODULES OF RAW ARE STRIPPED)
      SOUT: THE SOURCE PROGRAM WITH MODULES REPLACED FROM MODLIB
      MODCAT: AN ALPHABETIC INDEX TO MODLIB THAT IS RECREATED
         IF IT DOES NOT MATCH MODLIB
      LIST: PROGRESS OF THE TRANSFER.  MEANING OF THE LIST COLUMNS:
         NESTING DEPTH:  HOW DEEPLY THE MODULE WAS NESTED INSIDE OTHER MODULES
         ACTION:  WHAT WAS DONE WITH THE MODULE.  IF A MODULE WAS NOT
            TRANSFERRED, A SYMBOL ON THE LEFT FLAGS THE SITUATION:
              (BLANK) SUCCESSFUL TRANSFER
            * MODULE NOT FOUND IN THE SOURCE
            V NO TRANSFER BECAUSE VERSION MODULES CAN'T BE TRANSFERRED
            ? RECURSIVE TRANSFERS WERE ABORTED BECAUSE THE MODULES MAY BE
              INFINITELY NESTED (THE DEPTH AT WHICH THIS HAPPENS CAN BE
              INCREASED BY CHANGING THE PROGRAM - ASK YOUR PROGRAMMER).
              (PROBLEM: CAN YOU CONSTRUCT THIS BIZARRE INFINITE SITUATION?)
         MODULE NAME: THE NAME OF THE MODULE IN THE SOURCE.  IN RECURSIVE
            CASES, THESE ARE FROM THE MODLIB.
      output: MESSAGES TO THE USER
 
DESCRIPTION
      THE MODULE PROGRAM ALLOWS ONE TO CONSTRUCT LIBRARIES OF SPECIAL
      PURPOSE PROGRAM MODULES, WHICH ONE SIMPLY 'PLUGS' INTO THE
      APPROPRIATE PLACE IN A PROGRAM.  THIS SPEEDS UP BOTH PROGRAM DESIGN
      AND ERROR CORRECTION.  MODULE IS MORE GENERAL-PURPOSE THAN THE STANDARD
      'INCLUDE' type PROCESSES BECAUSE IT PERFORMS A REPLACEMENT RATHER THAN
      A SIMPLE INSERTION.  THE OPERATION IS RECURSIVE, SO A MODULE MAY BE
      COMPOSED OF OTHER MODULES.  THE REPLACEMENT MECHANISM ALSO ALLOWS ONE TO
      RUN THE PROGRAM IN 'REVERSE' SO THAT MODULE-LIBRARIES ARE CREATED BY
      EXTRACTING MODULES FROM EXISTING PROGRAMS.  THIS MAKES THE BUILDING OF
      MODULE LIBRARIES EASY, AND HELPS KEEP THEM UPDATED WITH NEW MODULES AND
      IMPROVEMENTS TO OLD ONES.
         FOR A FULL DESCRIPTION, SEE THE DOCUMENTATION.
 
DOCUMENTATION
      MODDEF, DELMAN.ASSEMBLY.MODULES,
      DELMAN.INTRO.ORGANIZATION 'TECHNICAL NOTES'
 
SEE ALSO
      DELMOD, PRGMOD, MATMOD, BREAK, SHOW (ESPECIALLY...)
 
AUTHOR
      THOMAS D. SCHNEIDER
 
BUGS
      NONE KNOWN
 
*)
(* END MODULE DELMAN.DESCRIBE.MODULE DELMAN 3.40 1982 SEP 1 SCHNEIDER-STORMO *)
 
(* MORE CONSTANTS *)
      LASTCHARACTER = ' '; (* THE LAST CHARACTER AFTER A MODULE NAME *)
      MAXNAME = 50; (* ONE PLUS THE LARGEST NAME ALLOWED *)
 
      (* MAXDEPTH IS THE LARGEST NUMBER OF RECURSIVE TRANSFERS ALLOWED BEFORE
THE PROGRAM ASSUMES THAT THERE MUST BE AN INFINITE NUMBER.  THE VALUE CAN BE
SET VERY HIGH, AN INFINITE EXAMPLE RUN (DEBUGGING SO THAT SOUT IS NOT
DESTROYED IN PROCEDURE HALT) AND THE NUMBER OF SUCCESSFUL TRANSFERS FOUND
AS THE NUMBER OF TRANSFERS SEEN IN SOUT.  THEN MAXDEPTH CAN BE SET TO A VALUE
SOMEWHAT UNDER THE TRUE MAXIMUM OF THE COMPUTER MEMORY. *)
      MAXDEPTH = 10;
 
      (* THE PROGRAM WILL CHECK THE CORRESPONDENCE BETWEEN MODLIB AND MODCAT.
         CHECKUPTIMES IS THE NUMBER OF MODULES TO CHECK.
         SEE THE CHECKUP PROCEDURE. *)
      CHECKUPTIMES = 2;
 
      (* IF THE CHECKUP FAILS, THE CONSTANT RECREATE DETERMINES
         WHAT WILL BE DONE:  HALT OR RECREATE THE CATALOGUE AND GO ON. *)
      RECREATE = true;
 
      MAXLINE = 70;
      
      DEBUGGING = false; (* FOR DEBUGGING PURPOSES *)
 
type
      NAME = record (* A MODULE NAME *)
         LETTER: packed array[1..MAXNAME] of char;
         LENGTH: integer; (* THE LAST CHARACTER OF THE NAME *)
      end;
 
      TRIGGER = record (* AN OBJECT TO BE SEARCHED FOR *)
         N: NAME; (* THE CHARACTERS LOOKED FOR *)
         STATE: integer; (* HOW CLOSE TO TRIGGERING WE ARE *)
         SKIP: boolean; (* TRIGGER NOT FOUND- SKIP THE LINE *)
         FOUND: boolean (* THE TRIGGER WAS FOUND *)
      end;
 
      MODCATITEM = record (* AN ITEM IN MODCAT *)
         MODNAME: NAME; (* The name of a module *)
         LINE: integer; (* THE LINE MODULE IS ON IN MODLIB *)
      end;
 
      MODCATFILE = FILE OF MODCATITEM;
 

 var  RAW, (* THE SOURCE IN FILE *)
      MODLIB, (* THE MODULE LIBRARY *)
      SOUT, (* THE SOURCE OUT FILE *)
         LIST: (* PROGRESS OF THE TRANSFER *)
         text;
 
      MODCAT: MODCATFILE; (* THE CATALOGUE FOR MODLIB *)
 
      RAWLINE: integer; (* THE CURRENT LINE IN RAW *)
      MODLIBLINE: integer; (* THE CURRENT LINE IN MODLIB *)
      
      (* THE TRIGGERS FOR THE MODULES *)
      BEGINTRIGGER, ENDTRIGGER: TRIGGER;
 
      RAWNAME: NAME; (* THE NAME OF THE TOP LEVEL *)
      VERMOD: NAME; (* A MODULE NAMED VERSION *)
 
      (* VARIABLES FOR KEEPING TRACK OF MODULE LIBRARY VERSION *)
      SHOWVERSION: boolean; (* IS THERE A VERSION TO SHOW? *)
      VERNAME: NAME; (* THE NAME OF THE VERSION TO SHOW *)
 
      (* VARIABLES USED FOR CHECKUP.
         THE MAIN PURPOSE FOR THESE TWO VARIABLES IS TO PREVENT A HALT WHEN
         PROCEDURE CHECKUP CALLS GETLINE AND THE LINE REFERED TO BY THE
         CATALOGUE IS PAST THE END OF THE MODLIB.  THIS ALLOWS RECOVERY FROM
         A SWITCH TO A SHORTER MODLIB WITHOUT CHANGING MODCAT:  A NEW
         MODCAT CAN BE CREATED. *)
      DONTHALT, (* THE HALT PROCEDURE SHOULD BE SILENT *)
      HALTCALLED: (* IF HALT IS CALLED AND DONTHALT IS true
                     THEN HALT WILL SET THIS TO true *)
                  boolean;
 
 
      (* THE NEXT TWO VARIABLES COUNT MODULES DETECTED AND TRANSFERED AT THE
      DEPTH=0.  (NOTE: ONE COULD COUNT INNER MODULES, BUT IN A TRIAL,
      THE OUTPUT BECAME CLUTTERED WITH ALMOST USELESS DATA.)  *)
      DETECTEDMODULES, (* NUMBER OF MODULES DETECTED *)
      TRANSFERREDMODULES: (* NUMBER OF MODULES TRANSFERRED *)
         integer;
 
(* HALT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: *)
 
procedure HALT;
(* STOP THE PROGRAM *)
begin
      if DONTHALT
      then HALTCALLED := true
      else begin
         (* THE FOLLOWING TWO LINES PREVENT PARTIALLY CREATED FILES FROM
            BEING ACCIDENTLY USED: *)
         rewrite(MODCAT);
(*         if not DEBUGGING then*) (* ALLOW ACCESS TO SOUT *)
         rewrite(SOUT);
         writeln(output, ' ERROR IN MODULE TRANSFER.  SEE LIST');
         writeln(output, ' PROGRAM HALT.');
         writeln(LIST, ' PROGRAM HALT.');
         goto 1
      end
end;
 
(* NONSTANDARD PROCEDURE TO ALLOW UNLIMITED OUTPUT :::::::::::::::::::::::: *)
(* BEGIN MODULE UNLIMITLN *)
procedure UNLIMITLN(var AFILE: text);
(* THIS PROCEDURE REMOVES A STUPID SYSTEM DEPENDENT LIMIT ON THE NUMBER OF
LINES THAT ONE CAN WRITE TO A FILE.  YOU MAY REMOVE IT FROM THE CODE IF
YOU SYSTEM DOES NOT WANT OR NEED THIS.  SUGGESTED METHOD: PLACE COMMENTS
AROUND THE CONTENTS OF THE PROCEDURE. *)
  begin
    (*LINELIMIT(AFILE, MAXINT);*)(* SET 'INFINITE' LINES ALLOWED FOR AFILE *)
  end;
(* END MODULE UNLIMITLN VERSION = 'DELMOD 6.12  82 MAR  8 TDS/GDS'; *)
 
(* CHARACTER AND LINE MANIPULATION :::::::::::::::::::::::::::::::::::::::: *)
 
procedure COPY (var FIN, FOUT: text; var CH: char);
(* COPY ONE CHARACTER (CH) FROM FILE FIN TO FILE FOUT *)
begin
      if not eof(FIN) then begin
         read(FIN, CH);
         write(FOUT, CH);
(*         if DEBUGGING then writeln(LIST,'COPY ',CH);*)
      end
end;
 
procedure FINISHLINE(var RAW, SOUT: text; var RAWLINE: integer);
(* FINISH COPY OF A LINE IN RAW TO SOUT, INCREMENT RAWLINE *)
var
      CH: char; (* ONE OF THE CHARACTERS COPIED *)
begin
      while not eoln(RAW) do COPY(RAW, SOUT, CH);
      readln(RAW);
      writeln(SOUT);
      RAWLINE := SUCC(RAWLINE)
end;
 
procedure GETTOLINE(LINE: integer; var F: text; var CURRENT: integer);
(* GET TO A LINE IN F FROM THE CURRENT PLACE.
   CURRENT = LINE AFTER GETTOLINE IS DONE *)
begin
      if CURRENT > LINE then begin (* THE LINE IS ABOVE WHERE WE ARE *)
         reset(F);
         CURRENT := 1
      end;
 
      while (CURRENT < LINE) and (not eof(F)) do begin
         CURRENT := SUCC(CURRENT);
         readln(F)
      end;
 
      if eof(F) then begin
         writeln(LIST, ' MODCAT REFERS TO A LINE (', LINE:1,
                        ') THAT IS PAST THE END OF MODLIB.');
         HALT
      end
end;
 
(* NAME MANIPULATIONS ::::::::::::::::::::::::::::::::::::::::::::::::::::: *)
 
procedure CLEARNAME(var N: NAME);
(* SET THE NAME TO BLANK *)
var   L: integer; (* A POSITION IN THE NAME *)
begin
      with N do begin
         for L := 1 to MAXNAME do LETTER[L] := ' ';
         LENGTH := 0
      end
end;
 
procedure PRINTNAME (var F: text; N: NAME);
(* PRINT THE NAME OF N TO FILE F *)
begin
      with N do write(F,LETTER:LENGTH)
end;
 
procedure UNTRAIL(var N: NAME);
(* REMOVE ALL TRAILING BLANKS FROM NAME N AND
   END N WITH THE LASTCHARACTER *)
begin
      with N do begin
         LENGTH := MAXNAME;
         while (LENGTH > 0) and (LETTER[LENGTH] = ' ')
            do LENGTH := PRED(LENGTH);
 
         (* PUT LASTCHARACTER AT end *)
         LENGTH := SUCC(LENGTH);
 
         if LENGTH > MAXNAME then begin
            LENGTH:=MAXNAME;
            write(LIST,' THIS NAME WAS FOUND: ');
            PRINTNAME(LIST, N);
            writeln(LIST,'.');
            writeln(LIST, ' NAMES MUST BE ONE CHARACTER SHORTER THAN ',
                           MAXNAME:1, ' CHARACTERS.');
            HALT
         end
         else N.LETTER[LENGTH] := LASTCHARACTER;
      end
end;
 
procedure GETNAME (var SOURCE: text; var N: NAME);
(* RETURN THE NAME N IN SOURCE, UPTO AND INCLUDING THE
GLOBAL CONSTANT LASTCHARACTER *)
var
      CH: char; (* ONE OF THE CHARACTERS IN N *)
begin
      CLEARNAME(N);
      CH := '.'; (* START WITH A CHARACTER THAT IS NOT LASTCHARACTER *)
      with N do begin
         while (not eoln(SOURCE)) and (CH <> LASTCHARACTER) and
               (LENGTH <= MAXNAME) do begin
            LENGTH := SUCC(LENGTH);
            read(SOURCE, CH);
            LETTER[LENGTH] := CH
         end;
         if LETTER[LENGTH] <> LASTCHARACTER then begin
            writeln(LIST, ' THIS MODULE NAME: ');
            PRINTNAME(LIST, N);
            if LENGTH = MAXNAME
            then writeln(LIST, ' IS TOO LONG (>', (MAXNAME - 1):1,
                                 ' CHARACTERS)')
                 (* IF NOT THAT, IT MUST BE eoln: *)
            else writeln(LIST, ' DID NOT END WITH A ', LASTCHARACTER, '.');
            HALT
         end
      end
end;
 
function EQUALNAME(A, B: NAME): boolean;
(* ARE THE NAMES A AND B THE SAME? *)
begin
      if A.LENGTH = B.LENGTH
      then EQUALNAME := (A.LETTER = B.LETTER)
      else EQUALNAME := false
end;
 
function GREATERNAME(A, B: NAME): boolean;
(* IS A ALPHABETICALLY AFTER B? *)
begin
      GREATERNAME := (A.LETTER > B.LETTER)
end;
 
(* MODULE MECHANISMS :::::::::::::::::::::::::::::::::::::::::::::::::::::: *)
 
procedure RESETTRIGGER(var T: TRIGGER);
(* RESET THE TRIGGER TO GROUND STATE *)
begin
      with T do begin
         STATE := 0;
         SKIP := false;
         FOUND := false
      end
end;
 
procedure TESTFORTRIGGER(CH: char; var T: TRIGGER);
(* LOOK AT THE CHARACTER CH.
   IF IT IS PART OF THE TRIGGER (AT THE CURRENT TRIGGER STATE),
       THEN THE TRIGGER STATE GOES HIGHER.
   IF IT IS NOT PART OF THE TRIGGER THEN THE TRIGGER STATE IS RESET,
      SKIP IS true AND ONE SHOULD SKIP ONWARD TO FIND THE TRIGGER.
   IF THE TRIGGER IS FOUND, FOUND IS true. *)
begin
      with T do begin
         STATE := SUCC(STATE);
(*         if DEBUGGING then begin
            PRINTNAME(LIST,N);
            writeln(LIST,'TESTFORTRIGGER N.LETTER[',STATE:1,']:',
                           N.LETTER[STATE],' CH:',CH);
         end;*)
         if N.LETTER[STATE] = CH
         then begin
            SKIP := false;
            if STATE = N.LENGTH then FOUND := true
                                else FOUND := false
         end
         else begin (* RESET TRIGGER *)
            STATE := 0;
            SKIP := true;
            FOUND := false
         end
      end
end;
 
procedure FINDMODULEEND(var RAW: text; MODNAME: NAME; var RAWLINE: integer);
(* FIND (BY READS) THE END OF THE MODULE IN RAW.  INCREMENT
RAWLINE, THE LINE IN RAW *)
var
      FOUND: boolean; (* THE MODULE END WAS FOUND *)
      CH: char; (* A CHARACTER IN RAW *)
      ENDNAME: NAME; (* PERHAPS THE END *)
begin
      FOUND := false;
      while (not FOUND) and (not eof(RAW)) do begin
         RESETTRIGGER(ENDTRIGGER);
         while not( eoln(RAW) or ENDTRIGGER.SKIP or ENDTRIGGER.FOUND) do begin
            read(RAW, CH);
            TESTFORTRIGGER(CH, ENDTRIGGER)
         end;
 
         if ENDTRIGGER.FOUND then begin
            GETNAME(RAW, ENDNAME);
            if EQUALNAME(ENDNAME, MODNAME) then FOUND := true;
         end;
 
         (* CLOSE RAW LINE UP *)
         readln(RAW);
         RAWLINE := SUCC(RAWLINE)
      end;
 
      if eof(RAW) and (not FOUND) then begin
         write(LIST, ' NO END TO MODULE ');
         PRINTNAME(LIST, MODNAME);
         writeln(LIST, ' WHOSE CONTENTS WERE BEING SKIPPED.');
         HALT
      end
end;
 
function COPYTOBOUND(var RAW, SOUT: text; var LINE: integer): char;
(* COPY FROM RAW TO SOUT UNTIL A MODULE BOUNDARY IS FOUND.
RETURN A CHARACTER:
      B  BEGIN MODULE FOUND
      E  END MODULE FOUND
      F  FILE END = eof FOUND
IN B OR E CASES, THE NAME IS TO BE PICKED UP NEXT. *)
var
      FOUND: boolean; (* A BOUNDARY WAS FOUND *)
      CH: char; (* ONE OF THE CHARACTERS IN RAW *)
begin
(*      if DEBUGGING then writeln(LIST,'COPYTOBOUND');*)
      FOUND := false;
 
      while (not FOUND) and (not eof(RAW)) do begin
         RESETTRIGGER(BEGINTRIGGER);
         RESETTRIGGER(ENDTRIGGER);
         while not (eoln(RAW) or ((BEGINTRIGGER.SKIP or BEGINTRIGGER.FOUND)
                   and (ENDTRIGGER.SKIP or ENDTRIGGER.FOUND) ) ) do begin
            COPY(RAW, SOUT, CH);
            TESTFORTRIGGER(CH, BEGINTRIGGER);
            TESTFORTRIGGER(CH, ENDTRIGGER)
         end;
 
         FOUND := BEGINTRIGGER.FOUND or ENDTRIGGER.FOUND;
 
         if not FOUND
            then if BEGINTRIGGER.SKIP or ENDTRIGGER.SKIP
               then while not eoln(RAW)
                  do COPY(RAW, SOUT, CH);  (* COPY REST OF LINE OUT *)
 
         if eoln(RAW) then begin
            readln(RAW);
            writeln(SOUT);
            LINE := SUCC(LINE)
         end
      end;
 
      if FOUND
      then begin
(*         if DEBUGGING then writeln(LIST,'COPYTOBOUND:FOUND');*)
         if BEGINTRIGGER.FOUND then COPYTOBOUND := 'B';
         if ENDTRIGGER.FOUND    then COPYTOBOUND := 'E'
      end
      else COPYTOBOUND := 'F' (* TERMINATION AT FILE END *)
(*      ;if DEBUGGING then writeln(LIST,'COPYTOBOUND') *)
end;
 
procedure COPYTOEND(var RAW, SOUT: text; MODNAME: NAME;
                    var RAWLINE: integer);
(* COPY TO THE END OF THE MODULE FROM RAW TO SOUT WITHOUT TRANSFERING
INNER MODULES, AND OBJECTING TO eof IN RAW.  INCREMENT RAWLINE *)
var
      DONE: boolean; (* DONE COPYING *)
      ENDNAME: NAME; (* A NAME OF A MODULE END, PERHAPS THAT OF MODULE *)
begin
(*      if DEBUGGING then writeln(LIST,'COPYTOEND'); *)
      DONE := false;
      while not DONE do begin
         case COPYTOBOUND(RAW, SOUT, RAWLINE) OF
            'B': ; (* IGNORE BEGINS *)
            'E': begin (* MAYBE THIS IS IT *)
               GETNAME(RAW, ENDNAME);
               if EQUALNAME(ENDNAME, MODNAME) then DONE := true;
               PRINTNAME(SOUT, ENDNAME);
               FINISHLINE(RAW, SOUT, RAWLINE)
            end;
            'F': begin
               write(LIST,' THE END OF MODULE ');
               PRINTNAME(LIST, MODNAME);
               writeln(LIST, ' WAS NOT FOUND DURING COPYING');
               HALT
            end
         end
      end
end;
 
procedure SKIPTOEND(var RAW, SOUT: text; var RAWLINE: integer);
(* SKIP TO THE END OF THE MODULE FOUND IN RAW.
HOWEVER, WE MUST FINISH THE LINE TO SOUT WHILE PICKING UP THE
MODULE NAME.  ALSO, THE LAST LINE OF THE MODULE MUST BE MADE. *)
var
      MODNAME: NAME; (* THE MODULE BEING SKIPPED *)
begin
(*      if DEBUGGING then writeln(LIST,'SKIPTOEND');*)
      (* OBTAIN THE MODULE NAME AND COPY THE LINE TO SOUT *)
      GETNAME(RAW, MODNAME);
      PRINTNAME(SOUT, MODNAME);
      FINISHLINE(RAW, SOUT, RAWLINE);
 
      if EQUALNAME(MODNAME, VERMOD) then
         (* WOAH THERE... WE CAN'T STRIP A VERSION MODULE... *)
         COPYTOEND(RAW, SOUT, MODNAME, RAWLINE)
      else begin
         (* SKIP OVER THE MODULE *)
         FINDMODULEEND(RAW, MODNAME, RAWLINE);
 
         (* AT THIS POINT ENDTRIGGER FOR MODULE MUST HAVE BEEN FOUND,
            BUT THE END OF THE MODULE WAS NOT WRITEN TO SOUT. *)
         PRINTNAME(SOUT, ENDTRIGGER.N); (* THE TRIGGER *)
         PRINTNAME(SOUT, MODNAME);       (* ITS NAME *)
(*         if DEBUGGING then write(SOUT,'(STRIPPED)'); *)
 
         (* PUT A BLANK AT THE END OF THE COMMENT: *)
         if LASTCHARACTER <> ' ' then write(SOUT, ' ');
         writeln(SOUT, '*)') (* END OF COMMENT *)
(*         ;if DEBUGGING then writeln(LIST, 'SKIPTOEND') *)
      end
end;
 
(* CATALOGUE MANIPULATIONS :::::::::::::::::::::::::::::::::::::::::::::::: *)
 
procedure GRAB(var F: MODCATFILE; var ITEM: MODCATITEM);
(* OBTAIN AN ITEM FROM FILE F *)
begin
      ITEM := F^;
      get(F)
end;
 
procedure DROP(var T: MODCATFILE; var ITEM: MODCATITEM);
(* PLACE AN ITEM INTO FILE T *)
begin
      T^ := ITEM;
      put(T)
end;
 
procedure SHOW(var O: text; var C: MODCATFILE);
(* SHOW THE MODCAT C ON FILE O *)
var
      ITEM: MODCATITEM; (* AN ITEM IN C *)
begin
      reset(C);
      writeln(O);
      writeln(O,'   LINE MODULE NAME');
      while not eof(C) do begin
         GRAB(C,ITEM);
         write(O,' ',ITEM.LINE:6,' ');
         PRINTNAME(O,ITEM.MODNAME);
         writeln(O)
      end;
      writeln(O)
end;
 
procedure BUILD(var MODLIB: text;
                var MODCAT: MODCATFILE);
(* BUILD THE MODCAT FROM THE MODLIB *)
var
      LI: integer; (* CURRENT LINE IN MODLIB *)
      CH: char; (* A CHARACTER IN MODLIB *)
      NA: NAME; (* A MODULE NAME *)
      ITEM: MODCATITEM; (* ONE OF THE RECORDS IN MODCAT *)
      NUMBER: integer; (* HOW MANY MODULES THERE ARE IN MODLIB *)
begin
      reset(MODLIB);
      rewrite(MODCAT);
      LI := 1;
      NUMBER := 0;
 
      while not eof(MODLIB) do begin
         RESETTRIGGER(BEGINTRIGGER);
         RESETTRIGGER(ENDTRIGGER);
         while not ( eoln(MODLIB) or
                     ((BEGINTRIGGER.SKIP or BEGINTRIGGER.FOUND) and
                      (  ENDTRIGGER.SKIP or   ENDTRIGGER.FOUND) ) ) do begin
            read(MODLIB, CH);
            TESTFORTRIGGER(CH, BEGINTRIGGER);
            TESTFORTRIGGER(CH, ENDTRIGGER)
         end;
 
         if BEGINTRIGGER.FOUND then with ITEM do begin
            GETNAME(MODLIB, MODNAME);
            LINE := LI;
            DROP(MODCAT, ITEM);
            NUMBER := SUCC(NUMBER); (* COUNT THE MODULES *)
            FINDMODULEEND(MODLIB, MODNAME, LI)
         end
         else if ENDTRIGGER.FOUND then begin
            write(LIST, ' UNEXPECTED MODULE END: ');
            GETNAME(MODLIB, NA);
            PRINTNAME(LIST, NA);
            writeln(LIST, ' AT LINE ', LI:1, ' IN MODLIB.');
            HALT
         end
         else begin
            readln(MODLIB);
            LI := SUCC(LI)
         end
      end;
 
      if NUMBER = 0
      then begin
         writeln(LIST, ' NO MODULES IN MODLIB.');
         HALT
      end
      else begin
         write(LIST, ' ', NUMBER:1, ' MODULE');
         if NUMBER <> 1 then write(LIST, 'S');
         writeln(LIST, ' IN MODLIB.');
         reset(MODLIB) (* CLEAN UP *)
      end
end; (* BUILD *)
 
procedure SORT(var F: MODCATFILE);
(* SORT THE FILE F.  A SIMPLE MULTIPLE PASS BUBBLE SORT IS USED
SINCE THE NUMBER OF ITEMS IN MODCAT IS OFTEN SMALL.  TWO FILES ARE
USED: F AND AN INTERNAL FILE (I) TO AVOID CONSTRAINTS OF AN ARRAY. *)
var
      I: MODCATFILE; (* AN INTERNAL FILE *)
      CHANGES: boolean; (* WHETHER CHANGES WERE MADE IN A PASS *)
procedure BUBBLEPASS(var F, T: MODCATFILE; var CHANGES: boolean);
(* PASS ONCE ACROSS FILE F COPYING TO FILE T.  INDICATE WHETHER
ANY SORTING HAPPENED USING CHANGES.  THE ALGORITHM IS SIMPLE:
PICKUP TWO ITEMS AND ALWAYS DROP THE SMALLER ONE. *)
var
      A, B: MODCATITEM; (* TWO OF THE ITEMS IN F *)
begin (* BUBBLEPASS *)
(*      if DEBUGGING then SHOW(LIST,F); *)
      CHANGES := false;
      reset(F);
      rewrite(T);
      GRAB(F, A);
(*      if DEBUGGING then write(LIST,'GRAB A,'); *)
      while not eof(F) do begin
(*         if DEBUGGING then write(LIST,'GRAB B,');*)
         GRAB(F, B);
 
         (* ALWAYS DROP THE SMALLER ITEM *)
         if GREATERNAME(B.MODNAME, A.MODNAME) or
              EQUALNAME(B.MODNAME, A.MODNAME)
         then begin
(*            if DEBUGGING then write(LIST,'DROP A, A:=B,');*)
            DROP(T, A);
            A := B (* REPLENISH A *)
         end
         else begin
            CHANGES := true;
(*            if DEBUGGING then write(LIST,'DROP B,');*)
            DROP(T, B) (* RETAIN A *)
         end
      end;
(*      if DEBUGGING then writeln(LIST,'DROP A.'); *)
      DROP(T, A) (* THE LAST ONE *)
end; (* BUBBLE PASS *)
begin (* SORT *)
      CHANGES := true;
(*!!!  ASSIGN(I,chr(0));*)(*req'd by IBM-PC to rewrite existing file*)
(*      if DEBUGGING then writeln(LIST,'SORT'); *)
      while CHANGES do begin
(*         if DEBUGGING then writeln(LIST,'PASS 1'); *)
         BUBBLEPASS(F, I, CHANGES);
(*         if DEBUGGING then writeln(LIST,'PASS 2(?)'); *)
         if CHANGES then BUBBLEPASS(I, F, CHANGES)
      end
(*; if DEBUGGING then writeln(LIST,'END SORT') *)
end; (* SORT *)
 
procedure CHECKDUPLICATE(var F: MODCATFILE);
(* CHECK FILE F FOR DUPLICATE NAMES, TAKING ADVANTAGE OF THE FACT
THAT IT IS SORTED *)
var
      A, B: MODCATITEM; (* TWO ITEMS IN F *)
      OK: boolean; (* NO DUPLICATES *)
begin
      reset(F);
      OK := true;
      GRAB(F, A);
      while not eof(F) do begin
         GRAB(F, B);
         if EQUALNAME(A.MODNAME, B.MODNAME) then begin
            OK := false;
            write(LIST, ' DUPLICATE MODULE NAME: ');
            PRINTNAME(LIST, A.MODNAME);
            writeln(LIST);
            writeln(LIST, '           FOUND AT LINES ',
               A.LINE:1, ' AND ',
               B.LINE:1, ' OF MODLIB.');
         end;
         A := B
      end;
      if not OK then begin
         rewrite(MODCAT); (* DESTROY THE BAD COPY *)
         HALT
      end
end;
 
procedure CREATEMODCAT(var MODLIB: text;
                       var MODCAT: MODCATFILE);
(* BUILD SORT AND CHECK THE MODULE CATALOGUE *)
begin
      writeln(LIST,' CREATING MODULE CATALOGUE (MODCAT)');
      BUILD(MODLIB, MODCAT);
      SORT(MODCAT);
      CHECKDUPLICATE(MODCAT);
      SHOW(LIST,MODCAT);
      reset(MODCAT)
end;
 
function INMODCAT(    MODNAME: NAME;
                  var LINE: integer): boolean;
(* IS THE MODULE IN THE MODCAT?  (MODCAT IS PASSED AS A GLOBAL FOR SPEED)
RETURN THE LINE NUMBER IN MODLIB (SIDE EFFECT) *)
var
      N: MODCATITEM; (* AN ITEM IN MODCAT *)
      FOUND: boolean; (* true WHEN MODULE IS FOUND *)
begin
      (* QUICK CHECK TO SEE if WE CAN AVOID A RESET *)
      if eof(MODCAT)
      then begin (* OH WELL... *)
         reset(MODCAT);
         FOUND := false
      end
      else begin (* WE STAND A CHANCE *)
         GRAB(MODCAT, N);
         if GREATERNAME(N.MODNAME, MODNAME)
         then begin (* IT IS ABOVE THIS POINT - WE LOSE *)
            reset(MODCAT);
            FOUND := false
         end
            (* IT IS BELOW THIS POINT - OR WE ARE ON IT *)
         else if EQUALNAME(N.MODNAME, MODNAME)
            then FOUND := true (* ZOOKS...GOT IT...ZOOKS...*)
            else FOUND := false (* IT IS BELOW THIS POINT - WE WIN *)
      end;
(*      if FOUND and DEBUGGING then writeln(LIST, 'ZOOKS...(INMODCAT)');*)
 
      while (not FOUND) and (not eof(MODCAT)) do begin
         GRAB(MODCAT, N);
         if EQUALNAME(N.MODNAME, MODNAME) then FOUND := true
      end;
 
      if FOUND then begin
         INMODCAT := true;
         LINE := N.LINE
      end
      else begin
         INMODCAT := false;
         LINE := -1 (* AN IMPOSSIBLE LINE NUMBER *)
      end
end; (* INMODCAT *)
 
procedure CHECKUP(var MODLIB: text;
                  var MODCAT: MODCATFILE);
(* CHECK THAT MODLIB CORRESPONDS TO MODCAT.
   THE NUMBER OF MODULES TO CHECK IS SET BY THE GLOBAL CONSTANT
   CHECKUPTIMES. *)
var
      TIMES: integer; (*  NUMBER OF CHECKS COMPLETED *)
      FAIL: boolean; (* WHAT MAY WELL HAPPEN DURING THIS CHECKUP *)
      CAT: MODCATITEM; (* ONE ITEM IN MODCAT *)
      CH: char; (* A CHARACTER FROM MODLIB *)
      LIBNAME: NAME; (* A NAME FROM MODLIB *)
      MODLIBLINE: integer; (* THE CURRENT LINE IN MODLIB *)
      MODCATLINE: integer; (* A LINE REFERED TO BY MODCAT *)
begin
      writeln(LIST, ' CHECK MODLIB-MODCAT CORRESPONDENCE:');
      reset(MODLIB);
      reset(MODCAT);
      MODLIBLINE := 1;
      TIMES := 0;
      FAIL := false;
 
      (* FIRST CHECK: DO ITEMS IN THE CATALOGUE POINT
         TO MODULES IN MODLIB *)
      DONTHALT := true; (* PREVENT HALTING DURING THIS CHECK *)
      repeat (* THIS FORCES AT LEAST ONE CHECK. *)
         (* GET AN ITEM FROM THE CATALOGUE *)
         GRAB(MODCAT, CAT);
 
         (* USE THE ITEM TO LOCATE A LINE IN MODLIB *)
         GETTOLINE(CAT.LINE, MODLIB, MODLIBLINE);
 
         if HALTCALLED
         then begin
            (* REFERENCE BY MODCAT TO A LINE PAST THE END OF MODLIB *)
            FAIL := true;
            HALTCALLED := false;
         end
         else begin
            (* FIRST WE MUST DETERMINE THAT A MODULE <BEGIN> IS THERE *)
            RESETTRIGGER(BEGINTRIGGER);
 
            while not (eoln(MODLIB) or
                       BEGINTRIGGER.FOUND or BEGINTRIGGER.SKIP) do begin
               read(MODLIB, CH);
               TESTFORTRIGGER(CH, BEGINTRIGGER)
            end;
 
            if BEGINTRIGGER.SKIP or eoln(MODLIB)
            then FAIL := true
            else (* BEGINTRIGGER.FOUND *) begin (* CHECK THE NAME *)
               GETNAME(MODLIB,LIBNAME);
               if not EQUALNAME(LIBNAME, CAT.MODNAME) then FAIL := true
            end;
 
            TIMES := SUCC(TIMES);
         end;
      until (TIMES >= CHECKUPTIMES) or (eof(MODCAT)) or FAIL;
      DONTHALT := false; (* ALLOW HALTING AGAIN *)
 
      (* SECOND CHECK: DO ITEMS IN MODLIB HAVE
         CORRESPONDING ITEMS IN MODCAT? *)
      if not FAIL then begin
(*         if DEBUGGING then writeln(LIST,' SECOND CHECK'); *)
         reset(MODLIB);
         reset(MODCAT);
         MODLIBLINE := 1;
         TIMES := 0;
 
         repeat
            RESETTRIGGER(BEGINTRIGGER);
            while not (eoln(MODLIB) or
                 BEGINTRIGGER.SKIP or BEGINTRIGGER.FOUND) do begin
               read(MODLIB, CH);
               TESTFORTRIGGER(CH, BEGINTRIGGER)
            end;
 
            if BEGINTRIGGER.FOUND then begin
               TIMES := SUCC(TIMES);
               GETNAME(MODLIB, LIBNAME);
               if not(INMODCAT(LIBNAME, MODCATLINE))
               then FAIL := true
                  (* MAYBE THE LINES DON'T MATCH... *)
               else if MODCATLINE <> MODLIBLINE then FAIL := true;
               if not FAIL then FINDMODULEEND(MODLIB, LIBNAME, MODLIBLINE)
            end
            else begin
               readln(MODLIB);
               MODLIBLINE := SUCC(MODLIBLINE)
            end
         until (TIMES > CHECKUPTIMES) or eof(MODLIB) or FAIL
      end;
 
      if FAIL
      then begin
         write(LIST,' FAILED: ');
         if RECREATE then CREATEMODCAT(MODLIB, MODCAT)
                     else HALT
      end
      else begin
         writeln(LIST,' PASSED.');
         MODLIBLINE:=1;
         reset(MODLIB);
         reset(MODCAT)
      end
end;

(* MAIN CALLS ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: *)
 
procedure INITIALIZE;

(* START UP THE PROGRAM *)

begin (* INITIALIZE *)            
      writeln(OUTPUT, ' ',VERSION);
     
      reset(RAW);
      reset(MODLIB);
      rewrite(SOUT);
      rewrite(LIST);
      reset(MODCAT);
      UNLIMITLN(LIST);
      UNLIMITLN(OUTPUT);
      UNLIMITLN(SOUT);
 
      writeln(LIST, ' ',VERSION);
       
      RAWLINE := 1;
      MODLIBLINE := 1;
 
      (* SET UP TRIGGERS.  THEY MUST BE THE SAME SIZE AS MAXNAME. *)
      with BEGINTRIGGER do begin
         (*                    1         2         3         4         5 *)
         (*           12345678901234567890123456789012345678901234567890 *)
         N.LETTER := '(* BEGIN MODULE                                   ';
         UNTRAIL(N)
      end;
 
      with ENDTRIGGER do begin
         (*                    1         2         3         4         5 *)
         (*           12345678901234567890123456789012345678901234567890 *)
         N.LETTER := '(* END MODULE                                     ';
         UNTRAIL(N)
      end;
 
      (* MAKE NAME OF TOP LEVEL.  THIS NAME MUST HAVE BLANKS IN IT
         TO AVOID DETECTION OF A MODULE BY THE SAME NAME *)
      (*                          1         2         3         4         5 *)
      (*                 12345678901234567890123456789012345678901234567890 *)
      RAWNAME.LETTER := '(SOURCE INPUT)                                    ';
      UNTRAIL(RAWNAME);
 
      (* MAKE NAME OF THE VERSION MODULE *)
      (*                          1         2         3         4         5 *)
      (*                 12345678901234567890123456789012345678901234567890 *)
      VERMOD.LETTER := 'VERSION                                           ';
      UNTRAIL(VERMOD);
 
      (* SET UP HALT VARIABLES *)
      DONTHALT := false;
      HALTCALLED := false;
 
      (* SET UP MODULE COUNTING VARIABLES *)
      DETECTEDMODULES := 0;
      TRANSFERREDMODULES := 0;
end; (* INITIALIZE *)
 
procedure GETVERSION;
(* FIND THE MODULE NAMED VERSION IN MODLIB, AND RETURN ITS
FIRST LINE IN VERNAME.  IF THERE IS NO MODULE, SET SHOWVERSION TO false *)
var
      LINE: integer; (* THE LINE IN MODLIB OF THE VERSION *)
      CH: char; (* A CHARACTER TO GO INTO VERNAME *)
begin (* GETVERSION *)
      if INMODCAT(VERMOD, LINE) then begin
         GETTOLINE(LINE, MODLIB, MODLIBLINE);
 
         (* MOVE TO FIRST LINE OF MODULE VERSION *)
         readln(MODLIB);
         MODLIBLINE := SUCC(MODLIBLINE);
 
         (* CAPTURE THE LINE *)
         CLEARNAME(VERNAME);
         with VERNAME do begin
            while (not eoln(MODLIB)) and
                  (LENGTH <= MAXNAME) do begin
               LENGTH := SUCC(LENGTH);
               read(MODLIB, CH);
               LETTER[LENGTH] := CH
            end
         end;
         readln(MODLIB);
         MODLIBLINE := SUCC(MODLIBLINE);
 
         UNTRAIL(VERNAME);
         SHOWVERSION := true
      end
      else SHOWVERSION := false
end; (* GETVERSION *)
 
procedure STRIP(var RAW, SOUT: text);
(* REMOVE MODULES IN RAW DURING COPY TO SOUT *)
var
      DONE: boolean; (* WHERE TO STOP *)
      RAWLINE: integer; (* LINE OF RAW *)
      ERROR: NAME; (* END NAME OF AN EXTRA MODULE END *)
begin
      writeln(LIST, ' NO MODULE LIBRARY (MODLIB): STRIPPING RAW TO SOUT.');
      writeln(OUTPUT, ' NO MODULE LIBRARY (MODLIB): STRIPPING RAW TO SOUT.');
      DONE := false;
      RAWLINE := 1;
 
      while not DONE do begin
         case COPYTOBOUND(RAW, SOUT, RAWLINE) of
            'B': begin
               SKIPTOEND(RAW, SOUT, RAWLINE);
               DETECTEDMODULES := SUCC(DETECTEDMODULES)
            end;
            'E': begin
               (* THERE MUST BE AN ERROR: THE B CASE DIDN'T CLOSE PROPERLY *)
               write(LIST, ' EXTRA MODULE END NAMED ');
               GETNAME(RAW, ERROR);
               PRINTNAME(LIST, ERROR);
               writeln(LIST, ' DETECTED AT LINE ', RAWLINE:1, ' OF RAW.');
               HALT
            end;
            'F': DONE := true
         end
      end
end;
 
function TRANSFER(   TMODULE: NAME;
                  var RAW, SOUT, MODLIB: text;
                  var RAWLINE: integer;
                  var MODLIBLINE: integer;
                      DEPTH: INTEGER): boolean;
(* COPY THE MODULE (NAMED TMODULE) FROM FILE RAW TO SOUT.
IT IS ASSUMED THAT THE FIRST LINE OF RAW (THE MODULE'S CALL LINE)
IS ALREADY COMPLETELY COPIED.
   IF FURTHER MODULE CALLS ARE SEEN, RECURSIVELY TRANSFER FROM MODLIB.
RETURN true IF THE MODULE END WAS FOUND, false IF END OF FILE WAS FOUND.
DEPTH KEEPS TRACK OF HOW DEEPLY WE HAVE RECURSED. *)
var
      DONE: boolean; (* true WHEN DONE *)
      ENDNAME: NAME; (* THE END OF A MODULE *)
procedure REPORT(var F: text; (* WHERE THE REPORT GOES *)
                     DEPTH: integer; (* NESTING DEPTH *)
                     WHAT: char; (* WHAT THE REPORT IS ABOUT *)
                     MODNAME: NAME); (* THE MODULE *)
(* REPORT TO FILE F WHAT HAPPENED TO THE MODULE AT SOME DEPTH
OF NESTING.  VALUES OF WHAT:
      T  TRANSFERRED
      N  NOT FOUND
      I  INFINITE RECURSION
      V  NO TRANSFER (THIS IS THE VERSION MODULE)
*)
begin (* REPORT *)
      write(F, ' ');
      case WHAT of
         'T': write(F,' ');
         'N': write(F,'*'); (* WARNING MARK FOR THE USER *)
         'I': write(F,'?'); (* INFINITE? *)
         'V': write(F,'V'); (* VERSION MODULE *)
      end;
      write(F, ' ', DEPTH:3, '    ');
      case WHAT of
         'T': write(F,'TRANSFERRED ');
         'N': write(F,'NOT FOUND   ');
         'I': write(F,'INFINITE??  ');
         'V': write(F,'NO TRANSFER ');
      end;
      PRINTNAME(F, MODNAME);
      writeln(F)
end; (* REPORT *)
procedure RECURSE;
(* TRANSFER THE INSIDES OF A MODULE *)
var
      INNER: NAME; (* NAME OF INNER MODULE *)
      LINE: integer; (* LINE NUMBER ON WHICH TO FIND INNER *)
      REMEMBER: integer; (* THE RAWLINE THAT WE MUST GET BACK TO
                            AFTER RECURSION *)
begin (* RECURSE *)
      GETNAME(RAW, INNER);
      PRINTNAME(SOUT, INNER);
      FINISHLINE(RAW, SOUT, RAWLINE);
      REMEMBER := RAWLINE;
      if DEPTH = 0 then DETECTEDMODULES := SUCC(DETECTEDMODULES);
 
      if INMODCAT(INNER, LINE) then begin (* IS A RECURSION POSSIBLE? *)
         if DEPTH >= MAXDEPTH then begin
            (* IT LOOKS LIKE THIS IS AN INFINITE RECURSIVE CALL,
               SO LETS KICK OUT. *)
            REPORT(LIST, DEPTH, 'I', INNER);
            writeln(SOUT,'(* THE MODULES ARE NESTED TO A DEPTH OF ',
                             (DEPTH + 1):1,' AT THIS POINT.');
            writeln(SOUT,'   PERHAPS THE MODLIB HAS AN INFINITE MODULE ',
                             'NESTING.');
            writeln(SOUT,'   FURTHER RECURSIVE TRANSFERS ARE ABORTED. *)');
            write  (output,' A POSSIBLE INFINITELY RECURSIVE NESTING OF');
            writeln(output,' MODULES WAS DETECTED.  SEE LIST.');
            COPYTOEND(RAW, SOUT, INNER, RAWLINE)
         end
         else if EQUALNAME(VERMOD, INNER) then begin
            (* IGNORE THE MODULE BECAUSE IT IS A VERSION MODULE *)
            REPORT(LIST, DEPTH, 'V', INNER);
            COPYTOEND(RAW, SOUT, INNER, RAWLINE)
         end
         else begin (* GO FOR A RECURSIVE TRANSFER *)
            (* LINE IS THE BEGINNING OF THE MODULE.
               SKIP THAT BY USING LINE + 1 *)
            GETTOLINE(LINE + 1, MODLIB, MODLIBLINE);
            if not TRANSFER(INNER, MODLIB, SOUT, MODLIB,
                            MODLIBLINE, MODLIBLINE, SUCC(DEPTH))
            then begin
               write(LIST, ' MISSING END OF MODULE ');
               PRINTNAME(LIST, INNER);
               writeln(LIST, ' IN MODLIB.');
               HALT
            end
            else begin
               (* THE INNER MODULE WAS INSERTED.  NOW WE MUST MOVE
                  BACK TO THE LINE FOLLOWING THE CALLING LINE: *)
               GETTOLINE(REMEMBER, RAW, RAWLINE);
 
               (* NOW SKIP THE REST OF THE CALLING MODULE *)
               FINDMODULEEND(RAW, INNER, RAWLINE);
 
               (*  CHALK ONE UP (TOP DEPTH ONLY) *)
               if DEPTH = 0 then TRANSFERREDMODULES := SUCC(TRANSFERREDMODULES)
            end
         end
      end
      else begin (* IGNORE THE MODULE SINCE IT IS NOT IN THE MODLIB *)
         REPORT(LIST, DEPTH, 'N', INNER);
         COPYTOEND(RAW, SOUT, INNER, RAWLINE)
      end
end; (* RECURSE *)
 
begin (* TRANSFER *)
      DONE := false;
      while not DONE do begin
         case COPYTOBOUND(RAW, SOUT, RAWLINE) of
            'B': RECURSE;
            'E': begin (* CHECK IF IT IS THE REAL END *)
               GETNAME(RAW, ENDNAME);
               PRINTNAME(SOUT, ENDNAME);
 
               (* SHOW VERSION OF MODLIB *)
               if SHOWVERSION then begin
                  PRINTNAME(SOUT, VERNAME);
                  writeln(SOUT, '*)');
                  readln(RAW); (* TOSS AWAY THE PREVIOUS STUFF... *)
                  RAWLINE := SUCC(RAWLINE)
               end
               else FINISHLINE(RAW, SOUT, RAWLINE);
 
               if EQUALNAME(ENDNAME, TMODULE) then begin
                  DONE := true;
                  TRANSFER := true
               end
               else begin
                  if DEPTH = 0 then begin
                     write(LIST, ' RAW MODULE ');
                     PRINTNAME(LIST, ENDNAME);
                     writeln(LIST, 'ENDED AT LINE ',(RAWLINE-1):1, '.');
                     writeln(LIST, ' THE BEGIN IS MISSING OR INCORRECT.')
                  end
                  else begin
                     write(LIST, ' MODULE BEGAN WITH THE NAME ');
                     PRINTNAME(LIST, TMODULE);
                     writeln(LIST, ',');
                     write(LIST, ' BUT ENDED WITH ');
                     PRINTNAME(LIST, ENDNAME);
                     write(LIST, ' AT LINE ', (RAWLINE-1):1, ' IN MODLIB.');
                  end;
                  HALT
               end
            end;
            'F': begin
               DONE := true;
               TRANSFER := false
            end
         end
      end;
      REPORT(LIST, DEPTH, 'T', TMODULE)
end;
 
begin (* MODULE *)
      INITIALIZE;
 
      if eof(RAW)
      then begin
         writeln(LIST, ' NO SOURCE (RAW) FILE.');
         HALT
      end
      else if eof(MODLIB) then STRIP(RAW, SOUT)
      else begin
         if eof(MODCAT) then CREATEMODCAT(MODLIB, MODCAT)
                        else CHECKUP(MODLIB, MODCAT);
 
         (* SET UP VERSION MECHANISM *)
         GETVERSION;
         if SHOWVERSION
         then begin
            write(LIST,' MODULE ');
            PRINTNAME(LIST,VERNAME);
            writeln(LIST)
         end
         else writeln(LIST,' NO VERSION FOR MODLIB.');
 
         writeln(LIST);
         writeln(LIST,' NESTING              MODULE');
         writeln(LIST,'  DEPTH   ACTION      NAME');
 
         (* DO THE TRANSFER OF RAWNAME *)
         if TRANSFER(RAWNAME, RAW, SOUT, MODLIB,
                     RAWLINE, MODLIBLINE, 0) then begin
            write(LIST, ' ZERO DEPTH MODULE NAME ');
            PRINTNAME(LIST, RAWNAME);
            writeln(LIST, ' DETECTED AS A MODULE - PROGRAM ERROR');
            HALT
         end
      end;
      writeln(OUTPUT,' ', DETECTEDMODULES:1,' MODULES DETECTED IN RAW, ',
                          TRANSFERREDMODULES:1, ' MODULES TRANSFERRED');
      writeln(LIST);
      writeln(LIST  ,' ', DETECTEDMODULES:1,' MODULES DETECTED IN RAW, ',
                          TRANSFERREDMODULES:1, ' MODULES TRANSFERRED');
      writeln(LIST);
1: end. (* MODULES *)

 
