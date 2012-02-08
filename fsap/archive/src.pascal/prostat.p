(*$DEBUG- *)
  (* ********************************************************  *)
  (*                                                           *)
  (*    PROSTAT  Version  8/30/93  Standard Pascal             *)
  (*             Brian Fristensky                              *)
  (*             Dept. of Plant Science                        *)
  (*             University of Manitoba                        *)
  (*             Winnipeg, MB  Canada  R3T 2N2                 *)
  (*                                                           *)
  (* Copyright (c) 1984 - 1993 by Brian Fristensky.            *)
  (* !!! in comment indicates feature which may need change.   *)
  (*  *******************************************************  *)
 
  program PROSTAT(input, output (*,SFILE,OUTFILE*));
 (*!!! Some Pascals require file parameters in program heading *)
 
  const MAXSEQ  =   10000;
        MAXRANGE=      0;
        MAXLINE =    150;
        MAXWORD =     25;
        VERSION = 'PROSTAT      Version  8/30/93';
 
  type AMINOACID = (GLY,ALA,VALINE,LEU,ILE,MET,PHE,PRO,SER,THR,CYS,ASN,GLN,TYR,
                    TRP,ASP,GLU,HIS,LYS,ARG,ASX,GLX,TERM,UNKX);
       PROTEIN   = array[1..MAXSEQ] of AMINOACID;
       AAI       = array[AMINOACID] of integer;
       AAR       = array[AMINOACID] of real;

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

  var SFILE,OUTFILE: text;
      SFN,OFN:LINE;
 
      (* Variables associated with the test protein *)
      SEQ:PROTEIN;
      SEQLEN:integer;
      NAME:WORD;
      MW:real;
      AACOMP       : AAI;    (* number of each amino acid in the protein *)
      AAWT         : AAR;    (* molecular weight of each amino acid *)

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

(* BEGIN MODULE READWORD *)
  (*  Read a word from a textfile           *)
  procedure READWORD(var F:text;var W:WORD);
    var I : integer;
        CH: char;
    begin
      with W do begin
        LEN:=0;
        while F^ = ' ' do
          if not eoln(F) then read(F,CH)
          else if not eof(F) then readln(F);
        while (F^ <> ' ') do 
          if LEN < MAXWORD then begin
            LEN := LEN + 1;
            read(F,STR[LEN])
            end
          else read(F,CH);
        for I := LEN+1 to MAXWORD do STR[I]:= ' '
        end
    end; (* READWORD *)
(* END MODULE READWORD         VERSION= 'SUNMODS     Version  9/12/91'; *)

(* BEGIN MODULE WRITEWORD *)
   (*  Write a word to a file using L char, left-justified.  *)
   procedure WRITEWORD(var F:text;W:WORD;L:integer);
     var I :integer;
     begin
       with W do
         for I := 1 to L do
           if I <= LEN then write(F,STR[I])
           else write(F,' ')
     end;
(* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  9/12/91'; *)

    (* **************************************************** *)
    (*              INITIALIZATION  PROCEDURES              *)
    (* **************************************************** *)
    procedure INITIALIZE;
      var AA1:AMINOACID;
      begin
        for AA1:= GLY to UNKX do AACOMP[AA1]:=0;

        (* Molecular weights of amino acids, water subtracted out *)
        AAWT[GLY]:= 57; AAWT[ALA]:= 71; AAWT[VALINE]:= 99; AAWT[LEU]:=113;
        AAWT[ILE]:=113; AAWT[MET]:=131; AAWT[PHE]:=147; AAWT[PRO]:= 97;
        AAWT[SER]:= 87; AAWT[THR]:=101; AAWT[CYS]:=103; AAWT[ASN]:=114;
        AAWT[GLN]:=128; AAWT[TYR]:=163; AAWT[TRP]:=186; AAWT[ASP]:=115;
        AAWT[GLU]:=129; AAWT[HIS]:=137; AAWT[LYS]:=128; AAWT[ARG]:=156;
        AAWT[ASX]:=114.5; AAWT[GLX]:=128.5; AAWT[TERM]:=  0; AAWT[UNKX]:=119
      end; (* INITIALIZE *)

(* BEGIN MODULE READPRO *)
    (*  ******************************************* *)
    (*  Read a protein    sequence from SFILE       *)
    (*   and store it in PR.                        *)
    (*  ******************************************* *)
    procedure READPRO(var SFILE:text; var PR:PROTEIN; var LEN:integer;
                      var NAME:WORD);
      label 86;

      var CH,ANSWER,FILETYPE: char;
          J: integer;
          TOOLONG:boolean;
 
      begin
        (* Prompt user for sequence file type *)
        writeln('The following file formats can be read:');
        writeln('  F:free format    N:NBRF');
        repeat
          writeln('Type letter of format (F|N)');
          readln(ANSWER)
        until ANSWER in ['F','f','N','n'];
        case ANSWER of
          'F','f':FILETYPE:='F';
          'N','n':FILETYPE:='N'
          end;
        writeln('Reading input file...');

        (* NBRF: read sequence name and title line *)
        (* 'Official' NBRF files have a name line followed by a title line.
            In the name line, four data characters preceed the name itself.
            These are deleted by the program.  Short NBRF files have a title
            following the name, on the same line. The name is not preceeded 
            by data characters. *)
        NAME.LEN:=0;
        if FILETYPE='N' then begin
          if SFILE^='>' then read(SFILE,CH);
          READWORD(SFILE,NAME);
          if NAME.STR[3]=';' then with NAME do begin
             for J:= 4 to LEN do STR[J-3]:=STR[J];
             LEN:=LEN-3;
             readln(SFILE)
             end; (* with NAME *)
          readln(SFILE)
          end; (* NBRF *)

        (* Read in the sequence *)
        J := 0;TOOLONG:=false;
        while not eof(SFILE) do begin
          while not eoln(SFILE) do begin
            read(SFILE,CH);
            if CH in ['G','A','V','L','I','M','F','P','S','T','C','N','Q',
                      'Y','W','D','E','H','K','R','B','Z','*','X'] then
              if J < MAXSEQ-MAXRANGE then begin
                J := J + 1;
                case CH of
                  'G':PR[J]:= GLY; 'A':PR[J]:= ALA; 'V':PR[J]:= VALINE;
                  'L':PR[J]:= LEU; 'I':PR[J]:= ILE; 'M':PR[J]:= MET;
                  'F':PR[J]:= PHE; 'P':PR[J]:= PRO; 'S':PR[J]:= SER;
                  'T':PR[J]:= THR; 'C':PR[J]:= CYS; 'N':PR[J]:= ASN;
                  'Q':PR[J]:= GLN; 'Y':PR[J]:= TYR; 'W':PR[J]:= TRP;
                  'D':PR[J]:= ASP; 'E':PR[J]:= GLU; 'H':PR[J]:= HIS;
                  'K':PR[J]:= LYS; 'R':PR[J]:= ARG; 'X':PR[J]:= UNKX; 
                  'B':PR[J]:= ASX; 'Z':PR[J]:= GLX; 
                  '*':begin  
                        PR[J]:= TERM;
                        if FILETYPE='N' then goto 86 (*ignore rest of file*)
                      end
                  end
                end
              else TOOLONG:= true
            else if CH=';' then readln(SFILE) (* comment in input file *)
            end;
          readln(SFILE)
          end;

(* !!! *) 86:(*CLOSE(SFILE); *) (* branch destination for end of NBRF seq. *)
        if TOOLONG then writeln(
          '>>> WARNING! Sequence exceeds MAXSEQ-MAXRANGE. Sequence truncated.');
        LEN := J
      end; (* READPRO *)
(* END MODULE READPRO         VERSION= 'SUNMODS     Version  9/12/91'; *)

    (***********************************************************************)
    (* Calculate amino acid composition and molecular weight of a protein. *)
    (***********************************************************************)
    procedure PROPARAM(var P:PROTEIN; SEQLEN:integer; var AACOMP:AAI;
                       var MW:real);
      const WATER = 18;
      var I:integer;
          AA:AMINOACID;
      begin
        (* Tabulate number of each amino acid in the protein *)
        for I:= 1 to SEQLEN do AACOMP[P[I]]:= AACOMP[P[I]]+1;
        
        (* Calculate molecular weight.*)     
        for AA:= GLY to UNKX do MW:= MW + (AACOMP[AA]*AAWT[AA]);
        MW:= MW + WATER (* accounts for hydration of N & C termini *)
      end; (* MW *)

    (**********************************************)
    (* Print a report of the findings.            *)
    (**********************************************)
    procedure REPORT;
      begin
        writeln('Type output filename:');
        GETFILE(OUTFILE,'O',OFN);
        writeln(OUTFILE,VERSION);
        WRITEWORD(OUTFILE,NAME,MAXWORD); writeln(OUTFILE,SEQLEN:10,' aa');
        writeln(OUTFILE,'Molecular weight:',MW:10:1);
        writeln(OUTFILE);
        writeln(OUTFILE,'Nonpolar side chains:');
        writeln(OUTFILE,'Gly(G)   ',AACOMP[GLY],' (',
                AACOMP[GLY]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Ala(A)   ',AACOMP[ALA],' (',
                AACOMP[ALA]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Val(V)   ',AACOMP[VALINE],' (',
                AACOMP[VALINE]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Leu(L)   ',AACOMP[LEU],' (',
                AACOMP[LEU]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Ile(I)   ',AACOMP[ILE],' (',
                AACOMP[ILE]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Met(M)   ',AACOMP[MET],' (',
                AACOMP[MET]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Phe(F)   ',AACOMP[PHE],' (',
                AACOMP[PHE]/SEQLEN:5:3,')');             
        writeln(OUTFILE,'Pro(P)   ',AACOMP[PRO],' (',
                AACOMP[PRO]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Neutral polar side chains:');
        writeln(OUTFILE,'Ser(S)   ',AACOMP[SER],' (',
                AACOMP[SER]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Thr(T)   ',AACOMP[THR],' (',
                AACOMP[THR]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Cys(C)   ',AACOMP[CYS],' (',
                AACOMP[CYS]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Asn(N)   ',AACOMP[ASN],' (',
                AACOMP[ASN]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Gln(Q)   ',AACOMP[GLN],' (',
                AACOMP[GLN]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Tyr(Y)   ',AACOMP[TYR],' (',
                AACOMP[TYR]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Trp(W)   ',AACOMP[TRP],' (',
                AACOMP[TRP]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Charged polar side chains:');
        writeln(OUTFILE,'Asp(D)   ',AACOMP[ASP],' (',
                AACOMP[ASP]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Glu(E)   ',AACOMP[GLU],' (',
                AACOMP[GLU]/SEQLEN:5:3,')');
        writeln(OUTFILE,'His(H)   ',AACOMP[HIS],' (',
                AACOMP[HIS]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Lys(K)   ',AACOMP[LYS],' (',
                AACOMP[LYS]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Arg(R)   ',AACOMP[ARG],' (',
                AACOMP[ARG]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Other:');
        writeln(OUTFILE,'Asx(B)   ',AACOMP[ASX],' (',
                AACOMP[ASX]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Glx(Z)   ',AACOMP[GLX],' (',
                AACOMP[GLX]/SEQLEN:5:3,')');
        writeln(OUTFILE,'Unk(X)   ',AACOMP[UNKX],' (',
                AACOMP[UNKX]/SEQLEN:5:3,')')
        end; (*REPORT*)
  
  (* ----------------------------------------------------------  *)
  (* ----------------- MAIN  PROCEDURE  -----------------------  *)
    begin
(* BEGIN MODULE STARTUP *)
(* Peform operations which must be done at beginning of main
   procedure. *)
(*!!!   TERMIN(input);    Open input for interactive use *)
(*!!!   TERMOUT(output);   "   output "      "        "  *)
     writeln(VERSION:50);
     writeln;
(* END MODULE STARTUP         VERSION= 'SUNMODS     Version  9/12/91'; *)
      writeln(VERSION);
      writeln;
 
      writeln('Enter sequence filename:');
      GETFILE(SFILE,'I',SFN);
      NAME.LEN:=0;
      READPRO(SFILE,SEQ,SEQLEN,NAME);
      if NAME.LEN=0 then begin
        writeln('Type name for protein:');
        INPWORD(NAME);readln
        end;
      INITIALIZE;
      if SEQ[SEQLEN]=TERM then SEQLEN:=SEQLEN-1;
      PROPARAM(SEQ,SEQLEN,AACOMP,MW);
      REPORT;
(*!!!*)CLOSE(OUTFILE)      
    end. (* PROSTAT *)
