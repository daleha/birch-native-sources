/*  The fix to makeexpr_hat was suggested by Dave Gillespie, as in the
    accompanying message:

From synaptx!thymus!daveg@uunet.UU.NET Thu Sep  9 19:42:05 1993
To: uunet!pssun1.plants.umanitoba.ca!frist@uunet.UU.NET (Brian Fristensky)
Subject: Re: p2c 
Reply-To: synaptx!daveg@uunet.UU.NET
Date: Thu, 09 Sep 93 17:34:16 -0700
From: Dave Gillespie <synaptx!thymus!daveg@uunet.UU.NET>
Content-Length: 1933
Status: RO
X-Lines: 44

> "eoftest.c", line 22: stdin_BFLAGS undefined

Hmm, this is a bit tricky.  P2c normally just uses a plain C
"FILE *" for file variables, but there are some things about
Pascal file semantics that "FILE *"'s can't express.  So it
then falls back on adding some extra variables with "_BFLAGS"
and "_BUFFER" appended to their names.  The problem is that
the "input^" expression accesses the current buffer of the
file "stdin", which unfortunately is a system file which
doesn't have a corresponding "_BFLAGS".  P2c apparently marks
the file variable as buffered anyway in its symbol table, and
then goes on to write a call to the "BUFEOF" macro as a
translation of "eof".

The irony is that, for text files, the "input^" notation is
translated as something that doesn't actually need the buffers.
The "BUFEOF" call is pointless and vestigial.

The function "makeexpr_hat" in expr.c checks for the case of
a file followed by "^" as its very first case.  It generates
calls to "chargetfbufname" (which will be "P_peek"),
"arraygetfbufname", or plain "getfbufname" depending on the
type of file.  It also calls "requirefilebuffer(a)" before
doing this, even if it was going to use "chargetfbufname".
The fix is to rearrange the code so that it only calls
"requirefilebuffer" for the non-"char" variants.

If that's too ambitious for you, you can approximate it by
commenting out the call to "requirefilebuffer" there.

								-- Dave

*/

Expr *makeexpr_hat(a, check)
Expr *a;
int check;
{
    Expr *ex;

    if (debug>2) { fprintf(outf,"makeexpr_hat("); dumpexpr(a); fprintf(outf,")\n"); }
    if (isfiletype(a->val.type, -1)) {
	if (*chargetfbufname &&
	    filebasetype(a->val.type)->kind == TK_CHAR)
	    return makeexpr_bicall_1(chargetfbufname,
				     filebasetype(a->val.type),
				     filebasename(a));
	else 
        {
	   requirefilebuffer(a);
           if (*arraygetfbufname &&
		 filebasetype(a->val.type)->kind == TK_ARRAY)
	    return makeexpr_bicall_2(arraygetfbufname,
				     filebasetype(a->val.type),
				     filebasename(a),
				     makeexpr_type(filebasetype(a->val.type)));
	else
	    return makeexpr_bicall_2(getfbufname,
				     filebasetype(a->val.type),
				     filebasename(a),
				     makeexpr_type(filebasetype(a->val.type)));
        }
    }
    if (a->kind == EK_PLUS &&
               (ex = a->args[0])->val.type->kind == TK_POINTER &&
               (ex->val.type->basetype->kind == TK_ARRAY ||
                ex->val.type->basetype->kind == TK_STRING ||
                ex->val.type->basetype->kind == TK_SET)) {
        ex->val.type = ex->val.type->basetype;   /* convert *(a+n) to a[n] */
        deletearg(&a, 0);
        if (a->nargs == 1)
            a = grabarg(a, 0);
        return makeexpr_bin(EK_INDEX, ex->val.type->basetype, ex, a);
    }
    if (a->val.type->kind == TK_STRING ||
        a->val.type->kind == TK_ARRAY ||
        a->val.type->kind == TK_SET) {
        if (starindex == 0)
            return makeexpr_bin(EK_INDEX, a->val.type->basetype, a, makeexpr_long(0));
        else
            return makeexpr_un(EK_HAT, a->val.type->basetype, a);
    }
    if (a->val.type->kind != TK_POINTER || !a->val.type->basetype) {
        warning("bad pointer dereference [165]");
        return a;
    }
    if (a->kind == EK_CAST &&
	a->val.type->basetype->kind == TK_POINTER &&
	a->args[0]->val.type->kind == TK_POINTER &&
	a->args[0]->val.type->basetype->kind == TK_POINTER) {
	return makeexpr_cast(makeexpr_hat(a->args[0], 0),
			     a->val.type->basetype);
    }
    switch (a->val.type->basetype->kind) {

      case TK_ARRAY:
      case TK_STRING:
      case TK_SET:
	if (a->kind != EK_HAT || 1 ||
	    a->val.type == a->args[0]->val.type->basetype) {
	    a->val.type = a->val.type->basetype;
	    return a;
	}

      default:
	if (a->kind == EK_ADDR) {
	    ex = a->args[0];
	    FREE(a);
	    return ex;
	} else {
	    if (check)
		ex = checknil(a);
	    else
		ex = a;
	    return makeexpr_un(EK_HAT, a->val.type->basetype, ex);
        }
    }
}

