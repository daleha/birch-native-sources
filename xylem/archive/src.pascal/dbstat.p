  (* ********************************************************  *)
  (*                                                           *)
  (*    DBSTAT   Version  2/ 1/90, Standard Pascal             *)
  (*             Brian Fristensky                              *)
  (*             Dept. of Plant Science                        *)
  (*             University of Manitoba                        *)
  (*             Winnipeg, MB R3T 2N2 CANADA                   *)
  (*                                                           *)
  (* Copyright (c) 1990 by Brian Fristensky.                   *)
  (* !!! in comment indicates feature which may need change.   *)
  (*  *******************************************************  *)
 
  program  DBSTAT(input,output);
 (*!!! Some Pascals require file parameters in program heading *)
 
  const VERSION = 'DBSTAT       Version  2/ 1/90';
 
  type AMINOACID = (G,A,V,L,I,M,F,P,S,T,C,N,Q,Y,W,D,E,H,K,R,B,Z,X);
       AAR       = array[AMINOACID] of integer;

  var 

      (* Variables associated with the test protein *)
      DBSIZE       : integer;
      AACOMP       : AAR;     (* number of each amino acid in the protein *)
      
    (* **************************************************** *)
    (*              INITIALIZATION  PROCEDURES              *)
    (* **************************************************** *)
    procedure INITIALIZE;
      var AA1:AMINOACID;
      begin
        for AA1:= G to X do AACOMP[AA1]:=0
      end; (* INITIALIZE *)

    (*  ******************************************* *)
    (*  Read a protein library and tabulate number  *)
    (*  of amino acids and composition of library.  *)
    (*  ******************************************* *)
    procedure READLIB(var SFILE:text; var DBSIZE:integer; var AACOMP:AAR);
      var CH: char;
      begin
        (* Read in the sequence *)
        DBSIZE:=0;
        while not eof(SFILE) do begin
          while not eoln(SFILE) do begin
            read(SFILE,CH);
            if CH in ['G','A','V','L','I','M','F','P','S','T','C','N','Q',
                      'Y','W','D','E','H','K','R','B','Z','X'] then begin
              DBSIZE:=DBSIZE+1;
              case CH of
                  'G':AACOMP[G]:=AACOMP[G]+1; 'A':AACOMP[A]:=AACOMP[A]+1; 
                  'V':AACOMP[V]:=AACOMP[V]+1; 'L':AACOMP[L]:=AACOMP[L]+1;
                  'I':AACOMP[I]:=AACOMP[I]+1; 'M':AACOMP[M]:=AACOMP[M]+1;
                  'F':AACOMP[F]:=AACOMP[F]+1; 'P':AACOMP[P]:=AACOMP[P]+1;
                  'S':AACOMP[S]:=AACOMP[S]+1; 'T':AACOMP[T]:=AACOMP[T]+1;
                  'C':AACOMP[C]:=AACOMP[C]+1; 'N':AACOMP[N]:=AACOMP[N]+1;
                  'Q':AACOMP[Q]:=AACOMP[Q]+1; 'Y':AACOMP[Y]:=AACOMP[Y]+1;
                  'W':AACOMP[W]:=AACOMP[W]+1; 'D':AACOMP[D]:=AACOMP[D]+1;
                  'E':AACOMP[E]:=AACOMP[E]+1; 'H':AACOMP[H]:=AACOMP[H]+1;
                  'K':AACOMP[K]:=AACOMP[K]+1; 'R':AACOMP[R]:=AACOMP[R]+1;
                  'B':AACOMP[B]:=AACOMP[B]+1; 'Z':AACOMP[Z]:=AACOMP[Z]+1;
                  'X':AACOMP[X]:=AACOMP[X]+1 
                  end
                end
            else if CH in [';','>'] then readln(SFILE) (* comment line *)
            end;
          if not eof(SFILE) then readln(SFILE)
          end
      end; (* READLIB *)
 
  
    (**********************************************)
    (* Print a report of the findings.            *)
    (**********************************************)
    procedure REPORT(var OUTFILE:text);
      begin
        writeln(OUTFILE,VERSION);
        writeln(OUTFILE);
        writeln(OUTFILE);
        writeln(OUTFILE,'Nonpolar side chains:');
        writeln(OUTFILE,'Gly(G)',AACOMP[G]:12,' (',AACOMP[G]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Ala(A)',AACOMP[A]:12,' (',AACOMP[A]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Val(V)',AACOMP[V]:12,' (',AACOMP[V]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Leu(L)',AACOMP[L]:12,' (',AACOMP[L]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Ile(I)',AACOMP[I]:12,' (',AACOMP[I]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Met(M)',AACOMP[M]:12,' (',AACOMP[M]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Phe(F)',AACOMP[F]:12,' (',AACOMP[F]/DBSIZE:5:3,')');   
        writeln(OUTFILE,'Pro(P)',AACOMP[P]:12,' (',AACOMP[P]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Neutral polar side chains:');
        writeln(OUTFILE,'Ser(S)',AACOMP[S]:12,' (',AACOMP[S]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Thr(T)',AACOMP[T]:12,' (',AACOMP[T]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Cys(C)',AACOMP[C]:12,' (',AACOMP[C]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Asn(N)',AACOMP[N]:12,' (',AACOMP[N]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Glu(Q)',AACOMP[Q]:12,' (',AACOMP[Q]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Tyr(Y)',AACOMP[Y]:12,' (',AACOMP[Y]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Trp(W)',AACOMP[W]:12,' (',AACOMP[W]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Charged polar side chains:');
        writeln(OUTFILE,'Asp(D)',AACOMP[D]:12,' (',AACOMP[D]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Glu(E)',AACOMP[E]:12,' (',AACOMP[E]/DBSIZE:5:3,')');
        writeln(OUTFILE,'His(H)',AACOMP[H]:12,' (',AACOMP[H]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Lys(K)',AACOMP[K]:12,' (',AACOMP[K]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Arg(R)',AACOMP[R]:12,' (',AACOMP[R]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Other:');
        writeln(OUTFILE,'Asx(B)',AACOMP[B]:12,' (',AACOMP[B]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Arg(Z)',AACOMP[Z]:12,' (',AACOMP[Z]/DBSIZE:5:3,')');
        writeln(OUTFILE,'Unk(X)',AACOMP[X]:12,' (',AACOMP[X]/DBSIZE:5:3,')')
        end; (*REPORT*)
  
  (* ----------------------------------------------------------  *)
  (* ----------------- MAIN  PROCEDURE  -----------------------  *)
    begin
      INITIALIZE;
      READLIB(input,DBSIZE,AACOMP);
      REPORT(output)
      
    end. (* DBSTAT *)

