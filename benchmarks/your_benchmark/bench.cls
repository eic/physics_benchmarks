%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         BENCH.CLS                              % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}
\LoadClass[11pt]{article}
\ProvidesClass{bench}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making the Title, with definitions for date, abstract etc %
% This part  written by David Cassel                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\@aabuffer{}
\def\author #1{\expandafter\def\expandafter\@aabuffer\expandafter
{\@aabuffer \small\rm      #1\relax \par}}
\def\address#1{\expandafter\def\expandafter\@aabuffer\expandafter
{\@aabuffer \small\it #1\relax 
\\
\Photo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\par\vspace{1em}}}

\def\maketitle{
\begin{center}
   {\bf \@title \par}       
   \vskip 2em                      % Vertical space after title.
   \@aabuffer\relax
\end{center} \par
\gdef\@aabuffer{}
}

\def\abstracts#1{
\begin{center}
{\begin{minipage}{5.2truein}
                 \footnotesize
                 \parindent=0pt #1\par
                 \end{minipage}}\end{center}
                 \vskip 2em \par}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Margins, textwidths, indentations etc %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\fussy
\flushbottom
\parindent 0.25in
%\oddsidemargin 4mm         
%\evensidemargin 4mm
%\topmargin=0mm 
%\headheight=0mm
%\headsep=0mm 
%\footskip=5mm           
%\textheight = 240mm     
%\textwidth = 160mm      
%\RequirePackage[a4paper,hmargin=1.5cm,vmargin={1.5cm,1.5cm}]{geometry}
\RequirePackage[a4paper,margin=2.5cm]{geometry}
\RequirePackage[english]{babel}
\RequirePackage{graphicx,url}
\RequirePackage[colorlinks=true,urlcolor=blue,linkcolor=black,citecolor=black]{hyperref}
% for BibTeX - sorted numerical labels by order of first citation.
\bibliographystyle{unsrt}

\def\section{\@startsection {section}{1}{\z@}{-3.5ex plus -1ex minus 
    -.2ex}{2.3ex plus .2ex}{\bf }}
\def\subsection{\@startsection{subsection}{2}{\z@}{-3.25ex plus -1ex minus 
   -.2ex}{1.5ex plus .2ex}{\it }}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   alpha footnotes, no running heads and silly citations.      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\thefootnote{\alph{footnote}}
\def\@makefnmark{{$\!^{\@thefnmark}$}}

\pagestyle{empty}


\renewenvironment{thebibliography}[1]
        {\begin{list}{\arabic{enumi}.}
        {\usecounter{enumi}\setlength{\parsep}{0pt}
         \setlength{\itemsep}{0pt} 
         \settowidth
        {\labelwidth}{#1.}\sloppy}}{\end{list}}

%---------------------------------------------------------------------------
%FOLLOWING THREE COMMANDS ARE FOR `LIST' COMMAND.
\topsep=0in\parsep=0in\itemsep=0in

\newcounter{arabiclistc}
\newenvironment{arabiclist}
        {\setcounter{arabiclistc}{0}
         \begin{list}{\arabic{arabiclistc}}
        {\usecounter{arabiclistc}
         \setlength{\parsep}{0pt}
         \setlength{\itemsep}{0pt}}}{\end{list}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACKNOWLEDGEMENT: this portion is from John Hershberger %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\@citex[#1]#2{\if@filesw\immediate\write\@auxout
        {\string\citation{#2}}\fi
\def\@citea{}\@cite{\@for\@citeb:=#2\do
        {\@citea\def\@citea{,}\@ifundefined
        {b@\@citeb}{{\bf ?}\@warning
        {Citation `\@citeb' on page \thepage \space undefined}}
        {\csname b@\@citeb\endcsname}}}{#1}}

\newif\if@cghi
\def\cite{\@cghitrue\@ifnextchar [{\@tempswatrue
        \@citex}{\@tempswafalse\@citex[]}}
\def\citelow{\@cghifalse\@ifnextchar [{\@tempswatrue
        \@citex}{\@tempswafalse\@citex[]}}
\def\@cite#1#2{{$\!^{#1}$\if@tempswa\typeout
        {IJCGA warning: optional citation argument 
        ignored: `#2'} \fi}}
\newcommand{\citeup}{\cite}

\setcounter{secnumdepth}{2}

\def\baselinestretch{1.0}
\ifx\selectfont\undefined
\@normalsize\else\let\glb@currsize=\relax\selectfont
\fi

\ifx\selectfont\undefined
\def\@singlespacing{%
\def\baselinestretch{1}\ifx\@currsize\normalsize\@normalsize\else\@currsize\fi%
}
\else
\def\@singlespacing{\def\baselinestretch{1}\let\glb@currsize=\relax\selectfont}
\fi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Footnote size table and figure captions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\long\def\@caption#1[#2]#3{\par\addcontentsline{\csname
  ext@#1\endcsname}{#1}{\protect\numberline{\csname
  the#1\endcsname}{\ignorespaces #2}}\begingroup
    \@parboxrestore
    \footnotesize
    \expandafter\let\expandafter\@tempa\csname @make#1caption\endcsname
    \ifx\@tempa\relax\let\@tempa\@makecaption\fi
    \@tempa{\csname fnum@#1\endcsname}{\ignorespaces #3}\par
  \endgroup}
%
% Here is the content of a .sty file containing definitions using the above
% hook.
%
% The following is the same as the \@makecaption in book.sty:
\long\def\@makefigurecaption#1#2{%
 \vskip 10pt 
 \setbox\@tempboxa\hbox{#1 -- {\footnotesize #2}}%
 \ifdim \wd\@tempboxa >\hsize #1 -- {\footnotesize #2}\par \else
  \hbox to\hsize{\hfil\box\@tempboxa\hfil} %
  \fi
 }
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           %  
%  END OF FILE bench.cls  %
%                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
