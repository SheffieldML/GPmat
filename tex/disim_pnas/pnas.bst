% pnas.bst

% BibTex format for Proceedings of the National Academy of Science
% this version matches the current PNAS format.

% history:
% version = 3.00 of pnas.bst 2009 Jun 28
% 2009 Jun 28 revised modernize emphasis function
% 2009 Jun 26 renamed pnas.bst by Tom Schneider to replace the 1999 version.
% 2008 Jan 30 received by Tom Schneider
% From: Ben Bolker <bolker@zoo.ufl.edu>
%%% BMB modified version from http://www.pnas.org/misc/iforc.shtml
%%%  28 Jan 2008 

%% References. References must be in PNAS style. Only published or
%% in-press papers and books may be cited in the reference
%% list. Unpublished abstracts of papers presented at meetings are not
%% permitted. References should be cited in numerical order as they
%% appear in text. Because tables and figures will be inserted in the
%% text where first cited, references in these sections should be
%% numbered accordingly. Manuscripts must include the full title for each
%% cited article. Authors must translate foreign language titles into
%% English, with a notation of the original language. All authors (unless
%% there are more than five) should be named in the citation. If there
%% are more than five, list the first author's name followed by et
%% al. Provide inclusive page ranges for journal articles and book
%% chapters. Databases are cited in the text or as footnotes.

%% The corresponding author must be prepared to provide a signed
%% authorization for the citation of unpublished data and personal
%% communications.
 
%% Journal articles are cited as follows:
 
%% 10.  Neuhaus J-M, Sitcher L, Meins F, Jr, Boller T (1991) A short C-terminal
%%        sequence is necessary and sufficient for the targeting of chitinases
%%        to the plant vacuole. Proc Natl Acad Sci USA 88:10362-10366.
%%  
%% For correct abbreviations of journal titles, refer to Chemical Abstracts Service Source Index (CASSI).
%%  
%% Articles or chapters in books are cited as follows:
%%  
%% 14.   Hill AVS (1991) in Molecular Evolution of the Major Histocompatibility
%%        Complex, eds Klein J, Klein D (Springer, Heidelberg), pp 403-420.
%%  

%%% ====================================================================
%%%  @BibTeX-style-file{
%%%     author          = "Alan R. Rogers",
%%%     version         = "1.0",
%%%     date            = "12 November 1992",
%%%     time            = "16:07:20 MST",
%%%     filename        = "pnas.bst",
%%%     address         = "Department of Anthropology
%%%                        University of Utah
%%%                        Salt Lake City, UT 84112
%%%                        USA",
%%%     telephone       = "+1 801 581 5529",
%%%     FAX             = "+1 801 581 6252",
%%%     checksum        = "37034 967 2316 17724",
%%%     email           = "rogers@anthro.utah.edu (Internet)",
%%%     codetable       = "ISO/ASCII",
%%%     keywords        = "",
%%%     supported       = "yes",
%%%     abstract        = "",
%%%     docstring       = "This file produces bibliographies in the
%%%		           format of ``Proceedings of the National
%%%		           Academy of Sciences, USA.''",
%%%  }
%%% ====================================================================
% This is "pnas.bst", modified from "unsrt.bst"
% BibTeX standard bibliography style `unsrt'
	% version 0.99a for BibTeX versions 0.99a or later, LaTeX version 2.09.
	% Copyright (C) 1985, all rights reserved.
	% Copying of this file is authorized only if either
	% (1) you make absolutely no changes to your copy, including name, or
	% (2) if you do make changes, you name it something other than
	% btxbst.doc, plain.bst, unsrt.bst, alpha.bst, and abbrv.bst.
	% This restriction helps ensure that all standard styles are identical.
	% The file btxbst.doc has the documentation for this style.

ENTRY
  { address
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
    month
    note
    number
    organization
    pages
    publisher
    school
    series
    title
    type
    volume
    year
  }
  {}
  { label }

INTEGERS { output.state before.all mid.sentence after.sentence after.block }

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}

STRINGS { s t }

FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
%    { ", " * write$ }
    { " " * write$ }
    { output.state after.block =
% BMB hack
	{ 
          % add.period$ 
          write$
	  newline$
	  "\newblock " write$
	}
	{ output.state before.all =
	    'write$
%	    { add.period$ " " * write$ }
	    { ", " * write$ }
	  if$
	}
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

FUNCTION {output.bibitem}
{ newline$
  "\bibitem{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {fin.entry}
{ add.period$
  write$
  newline$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}

FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
	'skip$
	{ after.sentence 'output.state := }
      if$
    }
  if$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {new.block.checka}
{ empty$
    'skip$
    'new.block
  if$
}

FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}

FUNCTION {new.sentence.checka}
{ empty$
    'skip$
    'new.sentence
  if$
}

FUNCTION {new.sentence.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.sentence
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
%    { "{\em " swap$ * "}" * }
% 2009 Jun 28: modernize the emphasis call
    { "\emph{" swap$ * "}" * }
  if$
}

FUNCTION {parenthesize}
{ duplicate$ empty$
    { pop$ "" }
    { "(" swap$ * ")" * }
  if$
}

FUNCTION {boldface}
{ duplicate$ empty$
    { pop$ "" }
    { "{\bf " swap$ * "}" * }
  if$
}

INTEGERS { nameptr namesleft numnames }

FUNCTION {format.names}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames #5 >   % et al for >5 names (copied from chicago.bst)
    { s #1 "{vv~}{ll}{, jj}{, f{}}" format.name$ " et~al." * }
    {
       numnames 'namesleft :=
         { namesleft #0 > }
            { s nameptr "{vv~}{ll}{, jj}{, f{}}" format.name$ 't :=
              nameptr #1 >
   	        { namesleft #1 >    %% >1 name left?
 	           { ", " * t * }  %% yes: add comma, next name
	           { t "others" =
		     { " et~al." * }
 		     { ", " * t * } %% add comma, last name
	             if$
                   }
                   if$
                }
               't
             if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
    while$  %% while namesleft>0
  }
  if$
}

FUNCTION {format.publisher}
{ publisher empty$
    { address empty$
      { "" }
      { address }
      if$
    }
    { address empty$
      { publisher}
      { publisher ", " * address * }
      if$
    }
  if$
}

FUNCTION {format.editors.as.authors}
{ editor empty$
    { organization }
    { editor format.names
      editor num.names$ #1 >
        {  ", eds" * } %% BMB removed period
        {  ", ed" * }  %% BMB ditto
      if$
    }
  if$
}

FUNCTION {format.authors}
{ author empty$
    { editor empty$
        { "" }
        { editor format.names
          editor num.names$ #1 >
           {  ", eds" * } %% BMB removed period
           {  ", ed" * }  %% BMB ditto
          if$
        }
      if$
    }
    { author format.names }
  if$
}

FUNCTION {format.editors}
{ editor empty$
    { "" }
    { author empty$ %in this case, see format.authors
        { "" }
        {  editor num.names$ #1 >
	     { "eds{} " }
	     { "ed{} " }
           if$
           editor format.names *
        }
      if$
    }
  if$
}

FUNCTION {format.title}
% BMB hack -- add period
{ title "t" change.case$ add.period$
}

FUNCTION {n.dashify}
{ 't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
	{ t #1 #2 substring$ "--" = not
	    { "--" *
	      t #2 global.max$ substring$ 't :=
	    }
	    {   { t #1 #1 substring$ "-" = }
		{ "-" *
		  t #2 global.max$ substring$ 't :=
		}
	      while$
	    }
	  if$
	}
	{ t #1 #1 substring$ *
	  t #2 global.max$ substring$ 't :=
	}
      if$
    }
  while$
}

FUNCTION {format.date}
{ year empty$
    { "year?" parenthesize    }
    { year parenthesize }
  if$
}

FUNCTION {format.btitle.series}
{ title empty$
    { "" }
    { title emphasize }
  if$
  series empty$
    { "" * }
    { ", " * series * }
  if$
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {comma.connect}
{ duplicate$ empty$
    { pop$ }
    { swap$
      duplicate$ empty$
        { pop$ }
        { swap$  ", " swap$ * * }
      if$
    }
  if$
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

FUNCTION {format.bvolume.number}
{ volume empty$
    { number empty$
        { "" }
        { "No.{}" number tie.or.space.connect }
      if$
    }
    { "Vol.{}" volume tie.or.space.connect
      "volume and number" number either.or.check
    }
  if$
}

FUNCTION {format.number}
{ volume empty$
    { number empty$
	{ "" }
	{ "No.{}" number tie.or.space.connect	}
      if$
    }
    { "" }
  if$
}

FUNCTION {format.edition}
{ edition empty$
    { "" }
    { output.state mid.sentence =
	{ edition "l" change.case$ " edition" * }
	{ edition "t" change.case$ " edition" * }
      if$
    }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
	{ #1 'multiresult := }
	{ t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}

FUNCTION {format.pages}
{ pages empty$
    { "" }
    { pages multi.page.check
	{ "pp" pages n.dashify tie.or.space.connect }
	{ "p" pages tie.or.space.connect }
      if$
    }
  if$
}

FUNCTION {format.vol.num.pages}
{ volume % boldface %% BB -- not boldface 
  pages empty$
    'skip$
    { duplicate$ empty$
	{ pop$ format.pages }
	{ ":" * pages n.dashify * }
      if$
    }
  if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages
    { type empty$
	{ "chapter" }
	{ type "l" change.case$ }
      if$
      chapter tie.or.space.connect
      pages empty$
	'skip$
	{ ", " * format.pages * }
      if$
    }
  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { "in " booktitle emphasize * }
  if$
  series empty$
    { "" * }
    { booktitle empty$
        { "XXmissing booktitle \& seriesXX" * }
        { ", " * }
      if$
      series *
    }
  if$
}

FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  month empty$ year empty$ note empty$
  and and and and and
    { "all relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type "t" change.case$
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { "Technical Report" }
    'type
  if$
  number empty$
    { "t" change.case$ }
    { number tie.or.space.connect }
  if$
}

FUNCTION {format.article.crossref}
{ key empty$
    { journal empty$
	{ "need key or journal for " cite$ * " to crossref " * crossref *
	  warning$
	  ""
	}
	{ "In {\em " journal * "\/}" * }
      if$
    }
    { "In " key * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.crossref.editor}
{ editor #1 "{vv~}{ll}" format.name$
  editor num.names$ duplicate$
  #2 >
    { pop$ " et~al." * }
    { #2 <
	'skip$
	{ editor #2 "{vv }{ll}{ jj}{, f}" format.name$ "others" =
	    { " et~al." * }
	    { " and " * editor #2 "{vv~}{ll}" format.name$ * }
	  if$
	}
      if$
    }
  if$
}

FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      "In "
    }
    { "Vol.{}" volume tie.or.space.connect
      " of " *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ series empty$
	    { "need editor, key, or series for " cite$ * " to crossref " *
	      crossref * warning$
	      "" *
	    }
	    { "{\em " * series * "\/}" * }
	  if$
	}
	{ key * }
      if$
    }
    { format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.incoll.inproc.crossref}
{ editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ booktitle empty$
	    { "need editor, key, or booktitle for " cite$ * " to crossref " *
	      crossref * warning$
	      ""
	    }
	    { "in {\em " booktitle * "\/}" * }
	  if$
	}
	{ "in " key * }
      if$
    }
    { "in " format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date "year" output.check
  format.title "title" output.check
  new.block
  crossref missing$
    { journal emphasize "journal" output.check
      format.vol.num.pages output
    }
    { format.article.crossref output.nonnull
      format.pages output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {book}
{ output.bibitem
  format.authors output.nonnull
  crossref missing$
    { "author and editor" editor either.or.check }
    'skip$
  if$
  new.block
  format.date "year" output.check
  format.btitle.series "title" output.check
  crossref missing$
    { format.editors output
      new.block
      format.publisher parenthesize output
      format.bvolume.number output
      new.sentence
    }
    { new.block
      format.book.crossref output.nonnull
    }
  if$
  format.edition output
  new.sentence
  format.pages output
  new.block
  note output
  fin.entry
}

FUNCTION {booklet}
{ output.bibitem
  format.authors output
  new.block
  format.date output
  format.title "title" output.check
  howpublished address comma.connect parenthesize output
  new.block
  note output
  fin.entry
}

FUNCTION {inbook}
{ output.bibitem
  format.authors output.nonnull
  new.block
  format.date "year" output.check
  format.btitle.series "title" output.check
  new.sentence
  crossref missing$
    { format.editors output
      new.block
      format.publisher parenthesize output
      format.bvolume.number output
      new.sentence
    }
    { format.chapter.pages "chapter and pages" output.check
      new.block
      format.book.crossref output.nonnull
    }
  if$
  format.edition output
  new.sentence
  format.pages output
  new.block
  note output
  fin.entry
}

FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date "year" output.check
  format.in.ed.booktitle "title" output.check
  new.sentence
  crossref missing$
    { format.editors output
      new.block
      format.publisher parenthesize output
      format.bvolume.number output
      new.sentence
    }
    { format.incoll.inproc.crossref output.nonnull
      format.chapter.pages output
    }
  if$
  format.edition output
  new.sentence
  format.pages output
  new.block
  note output
  fin.entry
}

FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date "year" output.check
  crossref missing$
    { format.btitle.series "title" output.check
      format.editors output
      new.block
      organization publisher comma.connect
      address comma.connect
      parenthesize output
      new.sentence
      format.bvolume.number output
      new.sentence
    }
    { format.incoll.inproc.crossref output.nonnull
    }
  if$
  format.pages output
  new.block
  note output
  fin.entry
}

FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  author empty$
    { organization empty$
	'skip$
	{ organization output.nonnull }
      if$
    }
    { format.authors output.nonnull }
  if$
  new.block
  format.date "year" output.check
  format.btitle.series "title" output.check
  author empty$
    { ""  }
    { organization }
  if$
  publisher comma.connect address comma.connect parenthesize output
  new.sentence
  format.edition output
  new.block
  note output
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date "year" output.check
  "Master's thesis" format.thesis.type output.nonnull
  school address comma.connect parenthesize output
  new.block
  note output
  fin.entry
}

FUNCTION {misc}
{ output.bibitem
  format.authors output
  new.block
  format.date output
  format.title output
  howpublished parenthesize output
  new.block
  note output
  fin.entry
  empty.misc.check
}

FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date "year" output.check
  "Ph.D. thesis" format.thesis.type output.nonnull
  school address comma.connect parenthesize output
  new.block
  note output
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  format.editors.as.authors output
  new.block
  format.date output
  format.btitle.series "title" output.check
  organization publisher comma.connect
  address comma.connect parenthesize output
  format.bvolume.number output
  new.sentence
  format.pages output
  new.block
  note output
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date output
  format.title "title" output.check
  new.sentence
  institution address comma.connect parenthesize output
  new.sentence
  format.tr.number output.nonnull
  new.block
  note output
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date output
  format.title "title" output.check
  new.block
  note "note" output.check
  fin.entry
}

FUNCTION {default.type} { misc }

MACRO {jan} {"January"}

MACRO {feb} {"February"}

MACRO {mar} {"March"}

MACRO {apr} {"April"}

MACRO {may} {"May"}

MACRO {jun} {"June"}

MACRO {jul} {"July"}

MACRO {aug} {"August"}

MACRO {sep} {"September"}

MACRO {oct} {"October"}

MACRO {nov} {"November"}

MACRO {dec} {"December"}

MACRO {acmcs} {"ACM Computing Surveys"}

MACRO {acta} {"Acta Informatica"}

MACRO {cacm} {"Communications of the ACM"}

MACRO {ibmjrd} {"IBM Journal of Research and Development"}

MACRO {ibmsj} {"IBM Systems Journal"}

MACRO {ieeese} {"IEEE Transactions on Software Engineering"}

MACRO {ieeetc} {"IEEE Transactions on Computers"}

MACRO {ieeetcad}
 {"IEEE Transactions on Computer-Aided Design of Integrated Circuits"}

MACRO {ipl} {"Information Processing Letters"}

MACRO {jacm} {"Journal of the ACM"}

MACRO {jcss} {"Journal of Computer and System Sciences"}

MACRO {scp} {"Science of Computer Programming"}

MACRO {sicomp} {"SIAM Journal on Computing"}

MACRO {tocs} {"ACM Transactions on Computer Systems"}

MACRO {tods} {"ACM Transactions on Database Systems"}

MACRO {tog} {"ACM Transactions on Graphics"}

MACRO {toms} {"ACM Transactions on Mathematical Software"}

MACRO {toois} {"ACM Transactions on Office Information Systems"}

MACRO {toplas} {"ACM Transactions on Programming Languages and Systems"}

MACRO {tcs} {"Theoretical Computer Science"}

READ

STRINGS { longest.label }

INTEGERS { number.label longest.label.width }

FUNCTION {initialize.longest.label}
{ "" 'longest.label :=
  #1 'number.label :=
  #0 'longest.label.width :=
}

FUNCTION {longest.label.pass}
{ number.label int.to.str$ 'label :=
  number.label #1 + 'number.label :=
  label width$ longest.label.width >
    { label 'longest.label :=
      label width$ 'longest.label.width :=
    }
    'skip$
  if$
}

EXECUTE {initialize.longest.label}

ITERATE {longest.label.pass}

FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{"  longest.label  * "}" * write$ newline$
}

EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}
