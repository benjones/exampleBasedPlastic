%
% acmsiggraph.cls
%
% version history
%
% - borrows *heavily* from Daniel Kartch's 'siggraph' class. 
% - also uses pieces of 'apalike.sty' and 'authordate1-4.cls'
% - modifications include adjustment of citation format to more closely
%   align with new (for SIGGRAPH 2002 conference and beyond) reference
%   and citation formatting requirements.
% - version 0.1 2001/06/01 Stephen Spencer (spencer@acm.org)
% - version 0.2 2001/10/15 Stephen Spencer
%   - the "\onlineid" variable works with the "review" mode, placing a
%     banner across the top of each page "Online Submission ID 'onlineid'
% - version 0.3 2002/01/11 Stephen Spencer
%   - disabled the cover page option.
% - version 0.4 2002/01/23 Stephen Spencer
%   - based on suggestions from James O'Brien, the following changes:
%     - correction if '\ifcamera' and '\ifcameraelse' commands
%     - page numbering in both review and preprint modes
%     have been made.
% - version 0.5 2005/01/22 Stephen Spencer
%   - removed the category markers
% - version 0.6 2006/12/15 Stephen Spencer
%   - tightened up space around section headings
% - version 0.7 2008/10/21 Stephen Spencer
%   - added 'annualconference' parameter for copyright space
%   - added line number support for 'review' mode
% - version 0.8 2011/03/02 Stephen Spencer (with Hugues Hoppe)
%   - removal of 'widereview' mode
%   - added hyperref
%   - added URL Links for paper DOI and author-provided content.
%
% ---------------------------------------------------------------------------

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{acmsiggraph-job}[2011/03/02 ACM SIGGRAPH class with JOB Mods]

% lots of options aren't allowed.

\DeclareOption{onecolumn}{\OptionNotUsed}
\DeclareOption{a4paper}{\OptionNotUsed}
\DeclareOption{12pt}{\OptionNotUsed}
\DeclareOption{14pt}{\OptionNotUsed}

% options that we DO want to set.

\newlength{\ACMpointsize}
\setlength{\ACMpointsize}{9pt}

% ``documentclass'' options.

\newif\if@ACMannual
\@ACMannualfalse
\DeclareOption{annual}{\@ACMannualtrue}

\newif\if@ACMsponsored
\@ACMsponsoredfalse
\DeclareOption{sponsored}{\@ACMsponsoredtrue}

\newif\if@ACMabstract
\@ACMabstractfalse
\DeclareOption{abstract}{\@ACMabstracttrue}

\newif\if@ACMreview
\@ACMreviewfalse
\DeclareOption{review}{\@ACMreviewtrue}

\newif\if@ACMpreprint
\@ACMpreprintfalse
\DeclareOption{preprint}{\@ACMpreprinttrue}

% -- Begin change by JOB
\newif\if@ACMreprint
\@ACMreprintfalse
\DeclareOption{reprint}{\@ACMreprinttrue}
% -- End change by JOB

% pass the options, and load the 'article' class, on which 'acmsiggraph' 
% is based.

\ProcessOptions\relax
\LoadClass[twocolumn,letterpaper,twoside]{article}

% required packages.

\RequirePackage{amsmath}
\RequirePackage{lineno}
\RequirePackage{lastpage}
\RequirePackage[scaled=.92]{helvet}
\RequirePackage{times}
\RequirePackage{graphicx}
\RequirePackage[labelfont=bf,textfont=it]{caption}
\RequirePackage{parskip}

% define things.

\setlength{\textheight}{9.250in}
\setlength{\textwidth}{7.000in}
\setlength{\topmargin}{-0.700in}
\setlength{\oddsidemargin}{-0.250in}
\setlength{\evensidemargin}{-0.250in}
\setlength{\columnsep}{2pc}
\setlength{\headheight}{0.200in}
\setlength{\headsep}{0.250in}
\setlength{\footskip}{0.500in}

% create a new length ``\ACMcopyrightspace'' and set its initial value.

\newlength{\ACMcopyrightspace}
\setlength{\ACMcopyrightspace}{3.3in}

% set the value of ``\ACMcopyrightspace'' based on the kind of document.

\if@ACMannual
\setlength{\ACMcopyrightspace}{1.5in}
\fi

\if@ACMsponsored
\setlength{\ACMcopyrightspace}{1.0in}
\fi

\if@ACMabstract
\setlength{\ACMcopyrightspace}{0in}
\fi

\if@ACMreprint
\setlength{\ACMcopyrightspace}{1.5in}
\fi

% -- Begin change by JOB
% Comment out setting lengths to zero to keep formatting consistent.
%\if@ACMreview
%\setlength{\ACMcopyrightspace}{0in}
%\fi
%
%\if@ACMpreprint
%\setlength{\ACMcopyrightspace}{0in}
%\fi
% -- End change by JOB




% define the ``\copyrightspace'' command.
% (todo) don't do this if ``\ACMcopyrightspace'' is zero.


% -- Begin change by JOB
% Many people like to put their email addresses above the copyright
% space, so the \abovecopyrightspacetext command is added.  ACM also
% wants a copyright notice to appear on submitted and preprint versions
% of the text so I added \insidecopyrightspacetext with appropriate default
% values. (Sorry about the long command names.)
\newcommand{\acm@defaultinsidecopyrightspacetext}{%
  \if@ACMreview
    \textsf{%
      This work has been submitted for publication. Copyright may be
      transferred without further notice and the accepted version may
      then be posted by the publisher.  Unauthorized viewing,
      duplication, and/or distribution are prohibited.
    }%
  \fi
  \if@ACMpreprint
    % This is not the right notice... it should be changed to whatever
    % ACM wants.
    \textsf{%
      This work has been submitted or accepted for publication. Copyright may be
      transferred without further notice and the accepted version may
      then be posted by the publisher. 
    }%
  \fi
  \if@ACMreprint
    \textsf{
      \vfill\vfill\textbf{\@ACMreprinttext}\vfill
      \begin{minipage}{\columnwidth}
        \renewcommand{\baselinestretch}{1}\fontsize{7}{7}\selectfont% 
        Permission to make digital or hard copies of all or part of
        this work for personal or classroom use is granted without fee
        provided that copies are not made or distributed for profit or
        commercial advantage and that copies bear this notice and the
        full citation on the first page. To copy otherwise, to
        republish, to post on servers or to redistribute to lists,
        requires prior specific permission and/or a fee.\\
        \textbf{\@ACMvenue, \@ACMyear\@ACMlocation}\\
        \copyright~Copyright ACM \@ACMyear\\
      \end{minipage}      
    }%
  \fi
}     
\newcommand{\acm@abovecopyrightspacetext}{}
\newcommand{\abovecopyrightspacetext}[1]{\renewcommand{\acm@abovecopyrightspacetext}{\par\noindent #1\\}}
\newcommand{\acm@insidecopyrightspacetext}{\acm@defaultinsidecopyrightspacetext}
\newcommand{\insidecopyrightspacetext}[1]{\renewcommand{\acm@insidecopyrightspacetext}{#1}}
\newcommand{\copyrightspace}{%
  \renewcommand{\thefootnote}{}%
  \footnotetext[0]{%
    \hspace*{-.225in}
    \parbox[t]{\columnwidth}{%
    \par\noindent\vspace{-2ex}%
    \acm@abovecopyrightspacetext%
    \parbox[t][\ACMcopyrightspace][b]{\columnwidth}{%
        \acm@insidecopyrightspacetext
 	\vspace{1.25ex}%
    }%
    }%
  }%
  \renewcommand{\thefootnote}{\arabic{footnote}}%
}
% -- Previous code 
%\newcommand{\copyrightspace}{%
%  \renewcommand{\thefootnote}{}%
%  \footnotetext[0]{\rule[\ACMcopyrightspace]{2.71828in}{0in}}%
%  \renewcommand{\thefootnote}{\arabic{footnote}}
%}
% -- End change by JOB



% in ``review'' mode, the ``\onlineid'' command passes the value given
% to the paper by the online submission system. This gets printed out 
% at the top of each page.

\newcommand{\TOGonlineid}[1]{\gdef\@ACMonlineid{#1}}
\newcommand{\@ACMreviewtext}{Online Submission ID: \@ACMonlineid}

% in ``preprint'' mode, the ``\TOGvolume'' and ``\TOGnumber'' commands 
% pass the volume and number of the TOG issue in which the paper will be
% published. This gets printed out at the top of each page.

\newcommand{\TOGvolume}[1]{\gdef\@ACMvolume{#1}}
\newcommand{\TOGnumber}[1]{\gdef\@ACMnumber{#1}}

% -- Begin change by JOB
\newcommand{\TOGyear}[1]{\gdef\@ACMyear{#1}}
\newcommand{\TOGmonth}[1]{\gdef\@ACMmonth{#1}}
\newcommand{\TOGvenue}[1]{\gdef\@ACMvenue{#1}}
\newcommand{\TOGarticlenum}[1]{\gdef\@ACMarticlenum{#1}}
\TOGyear{20XX}
\TOGmonth{}
\TOGvenue{ACM SIGGRAPH}
\gdef\@ACMlocation{}
\gdef\@ACMarticlenum{}
\newcommand{\TOGlocation}[1]{\@ACMconferencetrue\gdef\@ACMlocation{, #1}}
\newif\if@ACMconference
\@ACMconferencefalse
% -- End change by JOB

% -- Begin change by JOB
\newcommand{\@ACMpreprinttext}{To appear in \@ACMvenue\ \@ACMyear --- Draft version, please do not distribute.}
\newcommand{\@ACMreprinttext}{%
  \if@ACMconference
    {From the conference proceedings of \@ACMvenue\ \@ACMyear.\\%
     Appearing in ACM Transaction on Graphics Vol. \@ACMvolume, No. \@ACMnumber.}%\pageref{LastPage}
  \else
    {From ACM Transaction on Graphics \@ACMvolume(\@ACMnumber).}%\pageref{LastPage}
  \fi
}
% -- End change by JOB


%%%

\newcommand{\ifreview}[1]{\if@ACMreview #1 \fi}
\newcommand{\ifreviewelse}[2]{\if@ACMreview #1 \else #2 \fi}

%%% 

\newcommand{\@ACMarticleDOI}{}	% of the form 1073204.1073261
\newcommand{\@ACMarticleDOIsuffix}{} % of the form 1073261

\def\@ACMgetDOIsuffix#1.#2\@nil{#2}

\newcommand{\TOGarticleDOI}[1]{%
  \renewcommand{\@ACMarticleDOI}{#1}
  \in@{.}{#1}% given DOI of the form 1073204.1073261, if there is a period, extract the second part
  \ifin@
    \renewcommand{\@ACMarticleDOIsuffix}{\@ACMgetDOIsuffix#1\@nil}
  \else
    \renewcommand{\@ACMarticleDOIsuffix}{#1}
  \fi
}

%%%

\newcommand{\pdfauthor}[1]{\gdef\@ACMpdfauthor{#1}}
\newcommand{\pdftitle}[1]{\gdef\@ACMpdftitle{#1}}  % Added by JOB

% if we are in ``review'' or ``preprint'' modes, set up special headers 
% and footers for each page of the paper.

\renewcommand{\ps@plain}%
{%
    \renewcommand{\@oddhead}{}%
    \renewcommand{\@oddfoot}{}%

    \if@ACMpreprint%
      \renewcommand{\@oddhead}{\hfil\textit{\large\@ACMpreprinttext}\hfil}%
      \renewcommand{\@oddfoot}{\hfil\textrm{\thepage}\hfil}%
    \fi%

    \if@ACMreview%
      \renewcommand{\@oddhead}{\hfil\textit{\large\@ACMreviewtext}\hfil}%
      \renewcommand{\@oddfoot}{\hfil\textrm{\thepage}\hfil}%
    \fi%

    \renewcommand{\@evenhead}{\@oddhead}%
    \renewcommand{\@evenfoot}{\@oddfoot}%


    % -- Begin change by JOB
    \if@ACMreprint%
      \renewcommand{\@oddfoot}{\small \hfill ACM Transactions on Graphics, Vol. \@ACMvolume, No. \@ACMnumber, Article \@ACMarticlenum, Publication date: \@ACMmonth\ \@ACMyear.} %<<>>
      \renewcommand{\@evenfoot}{\small ACM Transactions on Graphics, Vol. \@ACMvolume, No. \@ACMnumber, Article \@ACMarticlenum, Publication date: \@ACMmonth\ \@ACMyear.\hfill} %<<>>


      \renewcommand{\@evenhead}{\textsf{\@ACMarticlenum:\thepage ~~~ $\cdot$ ~~~ \@ACMpdfauthor \hfill}}
      \renewcommand{\@oddhead}{\textsf{\hfill \@ACMpdftitle ~~~ $\cdot$ ~~~ \@ACMarticlenum:\thepage}}

    \fi%
    % -- End change by JOB


}

\renewcommand{\ps@empty}%
{%
    \renewcommand{\@oddhead}{}%
    \renewcommand{\@oddfoot}{}%

    % -- Begin change by JOB
    \if@ACMreprint%
      \renewcommand{\@oddfoot}{\small \hfill ACM Transactions on Graphics, Vol. \@ACMvolume, No. \@ACMnumber, Article \@ACMarticlenum, Publication date: \@ACMmonth\ \@ACMyear.} 
    \fi%
    % -- End change by JOB

    \if@ACMpreprint%
      \renewcommand{\@oddhead}{\hfil\textit{\large\@ACMpreprinttext}\hfil}%
      \renewcommand{\@oddfoot}{\hfil\textrm{\thepage}\hfil}%
    \fi%

    \if@ACMreview%
      \renewcommand{\@oddhead}{\hfil\textit{\large\@ACMreviewtext}\hfil}%
      \renewcommand{\@oddfoot}{\hfil\textrm{\thepage}\hfil}%
    \fi%

    \renewcommand{\@evenhead}{\@oddhead}%
    \renewcommand{\@evenfoot}{\@oddfoot}%
}

% no page numbers!

\pagestyle{plain}

% if we are in ``review'' mode, turn on line numbers.

\if@ACMreview%
  \linenumbers
\fi

% Nine-point type is the standard.

\newcommand{\@ACMsetninepoint}{
    \renewcommand\normalsize{%
        \@setfontsize\normalsize\@ixpt\@xpt
        \abovedisplayskip 9\p@ \@plus2\p@ \@minus4\p@
        \abovedisplayshortskip \z@ \@plus3\p@
        \belowdisplayshortskip 6\p@ \@plus3\p@ \@minus3\p@
        \belowdisplayskip \abovedisplayskip
        \let\@listi\@listI}
    \renewcommand\small{%
        \@setfontsize\small\@viipt\@ixpt
        \abovedisplayskip 8.5\p@ \@plus3\p@ \@minus4\p@
        \abovedisplayshortskip \z@ \@plus2\p@
        \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
        \def\@listi{\leftmargin\leftmargini
                    \topsep 4\p@ \@plus2\p@ \@minus2\p@
                    \parsep 2\p@ \@plus\p@ \@minus\p@
                    \itemsep \parsep}%
        \belowdisplayskip \abovedisplayskip}
    \renewcommand\footnotesize{%
        \@setfontsize\footnotesize\@viiipt{9.5}%
        \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
        \abovedisplayshortskip \z@ \@plus\p@
        \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
        \def\@listi{\leftmargin\leftmargini
                    \topsep 3\p@ \@plus\p@ \@minus\p@
                    \parsep 2\p@ \@plus\p@ \@minus\p@
                    \itemsep \parsep}%
        \belowdisplayskip \abovedisplayskip}
    \renewcommand\scriptsize{\@setfontsize\scriptsize\@viipt\@viiipt}
    \renewcommand\tiny{\@setfontsize\tiny\@vpt\@vipt}
    \renewcommand\large{\@setfontsize\large\@xpt\@xiipt}
    \renewcommand\Large{\@setfontsize\Large\@xiipt{14}}
    \renewcommand\LARGE{\@setfontsize\LARGE\@xivpt{18}}
    \renewcommand\huge{\@setfontsize\huge\@xviipt{22}}
    \renewcommand\Huge{\@setfontsize\Huge\@xxpt{25}}
    \selectfont
}

\@ACMsetninepoint

% Section and subsection environment tweaks.

\newcommand{\@ACMsectionfont}{}
\newcommand{\sectionfont}[1]{\renewcommand{\@ACMsectionfont}{#1}}
\renewcommand\section{\@startsection {section}{1}{\z@}%
                {-1.0ex}%
                {1.0ex}%
                {\reset@font\Large\sffamily\bfseries\@ACMsectionfont}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                {-1.0ex}%
                {1.0ex}%
                {\reset@font\large\sffamily\bfseries\@ACMsectionfont}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                {-1.0ex}%
                {1.0ex}%
                {\reset@font\normalsize\sffamily\bfseries\@ACMsectionfont}}

% Paragraph and subparagraph environment tweaks.

\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                {2.05ex}%
                {-1em}%
                {\reset@font\normalsize\sffamily\bfseries\@ACMsectionfont}}
\renewcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                {3.25ex \@plus1ex \@minus .2ex}%
                {-1em}%
                {\reset@font\normalsize\sffamily\bfseries\@ACMsectionfont}}

% no labels in bibliography.

\def\@biblabel#1{}

% bibitems are flush left, with successive lines indented by 1em.

\newlength{\bibhang}
\setlength{\bibhang}{1em}

\def\thebibliography#1{%
  \section*{%
    \refname\@mkboth{\sl\uppercase{\refname}}{\sl\uppercase{\refname}}}
  \list{\relax}{\setlength{\labelsep}{0em}
                \setlength{\itemindent}{-\bibhang}
                \setlength{\leftmargin}{\bibhang}}
  \def\newblock{\hskip .11em plus .33em minus .07em}
  \sloppy\clubpenalty4000\widowpenalty4000
  \sfcode`\.=1000\relax}

\def\@citex[#1]#2{\if@filesw\immediate\write\@auxout{\string\citation{#2}}\fi
  \def\@citea{}\@cite{\@for\@citeb:=#2\do
    {\@citea\def\@citea{; }\@ifundefined
      {b@\@citeb}{{\bf ?}\@warning
      {Citation '\@citeb' on page \thepage \space undefined}}%
{\csname b@\@citeb\endcsname}}}{#1}}

\let\@internalcite\cite




\def\cite{\def\citename##1{##1 }\def\citeyear##1{##1}\@internalcite}
\def\shortcite{\def\citename##1{}\def\citeyear##1{##1}\@internalcite}

\def\nouncite#1{% 
  \def\citename##1{##1 }\def\citeyear##1{}\csname b@#1\endcsname%
  \def\citename##1{}\def\citeyear##1{##1}\@internalcite{#1}}



% title stuff...

\newcommand{\@ACMempty}{}

\newcommand{\@ACMaffiliation}{}
\newcommand{\affiliation}[1]{\gdef\@ACMaffiliation{#1}}

\newcommand{\@ACMteaser}{}
\newcommand{\teaser}[1]{\gdef\@ACMteaser{#1}}

\newlength{\titlespace}
\setlength{\titlespace}{0.25in}

\newlength{\teaserspace}
\setlength{\teaserspace}{0.25in}

\if@ACMreview
    \renewcommand{\@maketitle}{%
        \begin{center}%
            \renewcommand{\thanks}[1]{}
            {\LARGE\sffamily\bfseries\@ACMsectionfont \@title \par}%
            \vspace{1\baselineskip}%
            \vspace{\titlespace}%
            \ifx\@ACMteaser\@ACMempty \else%
                \begingroup%
                    \def\@captype{figure}%
                    \@ACMteaser%
                \endgroup\par%
                \vspace{\teaserspace}%
            \fi%
        \end{center} \par%
    }
\else
    \renewcommand{\@maketitle}{%
        \begin{center}%
            {\LARGE\sffamily\bfseries\@ACMsectionfont \@title \par}%
            \vspace{1\baselineskip}%
            \large
            \begin{tabular}[t]{c}%
                \@author
            \end{tabular}\par%
            \ifx\@ACMaffiliation\@ACMempty \else%
                \vspace{0.5\baselineskip}%
                \@ACMaffiliation\par%
            \fi%
            \vspace{\titlespace}%
            \ifx\@ACMteaser\@ACMempty \else%
                \begingroup%
                    \def\@captype{figure}%
                    \@ACMteaser%
                \endgroup\par%
                \vspace{\teaserspace}%
            \fi%
        \end{center} \par%
    }
\fi

\let\@ACMorigmaketitle\maketitle
\let\@ACMorigand\and
\renewcommand{\maketitle}{%
    \let\@ACMtitle\@title%
    \let\@ACMauthor\@author%
    \@ACMorigmaketitle%
    \thispagestyle{empty}%
}

% Abstract stuff...

\newtoks\@ACMabs
\if@ACMreview
    \long\def\@ACMadd#1{\global\@ACMabs\expandafter{\the\@ACMabs#1}}
    \long\def\@ACMcollect{%
        \global\@ACMabs{}%
        \let\abstract\@ACM@collect
        \abstract
    }
    \long\def\@ACM@collect#1\end#2{%
        \def\@tempa{#2}%
        \ifx\@tempa\@currenvir
            \@ACMadd{#1}%
            \edef\abstract{\noexpand\end{\@tempa}}%
        \else
            \@ACMadd{#1\end{#2}}%
        \fi
        \abstract
    }
    \newcommand{\@ACMmodify}{%
        \let\@ACMorigabs\abstract%
        \let\@ACMorigendabs\endabstract%
        \renewenvironment{abstract}%
            {\@ACMcollect}%
            {\begingroup
                \let\abstract\@ACMorigabs
                \let\endabstract\@ACMorigendabs
                \begin{abstract} \the\@ACMabs \end{abstract}
            \endgroup}%
    }
    \AtBeginDocument{\@ACMmodify}
\fi

%%% Keyword list.

\newcommand{\@ACMkeywords}{}
\newcommand{\keywords}[1]{\renewcommand{\@ACMkeywords}{#1}}

\newcommand{\keywordlist}{%
    \vspace{0.5\baselineskip}%
    \par\noindent \textbf{Keywords: } \@ACMkeywords%
}

%%% Computing Review (CR) categories.

\newenvironment{CRcatlist}{%
    \vspace{0.5\baselineskip}%
    \par\noindent \textbf{CR Categories: }%
}{}

\newcommand{\CRcat}[4]{%
    #1 [#2]%
    \ifx#3\@ACMempty \else : #3\fi%
    \ifx#4\@ACMempty \else ---#4\fi%
}

%%% Icon(s) with hyperlinks.



\newcommand{\@ACMurlWeb}{}
\newcommand{\TOGprojectURL}[1]{\renewcommand{\@ACMurlWeb}{#1}}
\newcommand{\@ACMurlVideo}{}
\newcommand{\TOGvideoURL}[1]{\renewcommand{\@ACMurlVideo}{#1}}
\newcommand{\@ACMurlData}{}
\newcommand{\TOGdataURL}[1]{\renewcommand{\@ACMurlData}{#1}}
\newcommand{\@ACMurlCode}{}
\newcommand{\TOGcodeURL}[1]{\renewcommand{\@ACMurlCode}{#1}}

% -- Begin change by JOB
\newcommand{\@ACMurlPDF}{}
\newcommand{\TOGpaperURL}[1]{\renewcommand{\@ACMurlPDF}{#1}}

\newcommand{\@ACMurlWebData}{}
\newcommand{\TOGwebdataURL}[1]{\renewcommand{\@ACMurlWebData}{#1}}

% -- End change by JOB


\newcommand{\icon@DL}{\raisebox{-.02in}{\includegraphics[width=.13in]{Icons/iconacm}}\hspace{.01in}DL}
\newcommand{\icon@PDF}{\raisebox{-.02in}{\includegraphics[width=.11in]{Icons/iconpdf}}\hspace{.015in}PDF}
\newcommand{\icon@Web}{\raisebox{-.025in}{\includegraphics[width=.11in]{Icons/iconweb}}\hspace{.01in}Web}
\newcommand{\icon@WebData}{\raisebox{-.025in}{\includegraphics[width=.11in]{Icons/iconweb}}\hspace{.01in}Website and Data}
\newcommand{\icon@Video}{\raisebox{-.01in}{\includegraphics[width=.11in]{Icons/iconvideo}}\hspace{.007in}Video}
\newcommand{\icon@Data}{\raisebox{-.01in}{\includegraphics[width=.10in]{Icons/icondata}}\hspace{.015in}Data}
\newcommand{\icon@Code}{\raisebox{-.005in}{\includegraphics[width=.10in]{Icons/iconcode}}\hspace{.017in}Code}

\newcommand{\TOGlinkslist}{%
  \if@ACMannual
  {%
    \ifx\@ACMarticleDOI\@ACMempty%
	  \ClassWarningNoLine{acmsiggraph}{Please set \noexpand\TOGarticleDOI.}%
	\fi%
    {\bf Links:}%
    \sc%
    \ifx\@ACMarticleDOI\@ACMempty\else%
      \hspace{0.09in}\href{http://doi.acm.org/10.1145/\@ACMarticleDOI}{\icon@DL}%
    \ifx\@ACMurlPDF\@ACMempty
       \ifx\@ACMarticleDOIsuffix\@ACMempty\else%
       \hspace{0.09in}\href{http://portal.acm.org/ft_gateway.cfm?id=\@ACMarticleDOIsuffix&type=pdf}{\icon@PDF}\fi%
    \else
       \hspace{0.09in}\href{\@ACMurlPDF}{\icon@PDF}\fi%
     \ifx\@ACMurlVideo\@ACMempty\else%
       \hspace{0.09in}\href{\@ACMurlVideo}{\icon@Video}\fi%
     \ifx\@ACMurlWebData\@ACMempty\else%
       \hspace{0.09in}\href{\@ACMurlWebData}{\icon@WebData}\fi%
     \ifx\@ACMurlWeb\@ACMempty\else%
       \hspace{0.09in}\href{\@ACMurlWeb}{\icon@Web}\fi%
     \ifx\@ACMurlData\@ACMempty\else%
       \hspace{0.09in}\href{\@ACMurlData}{\icon@Data}\fi%
     \ifx\@ACMurlCode\@ACMempty\else%
       \hspace{0.09in}\href{\@ACMurlCode}{\icon@Code}\fi%
    \fi%
  }%
  \fi
  \if@ACMreprint
  {%
    \ifx\@ACMarticleDOI\@ACMempty%
	  \ClassWarningNoLine{acmsiggraph}{Please set \noexpand\TOGarticleDOI.}%
	\fi%
    {\bf Links:}%
    \sc%
    \ifx\@ACMarticleDOI\@ACMempty\else%
      \hspace{0.09in}\href{http://doi.acm.org/10.1145/\@ACMarticleDOI}{\icon@DL}%
    \ifx\@ACMurlPDF\@ACMempty
       \ifx\@ACMarticleDOIsuffix\@ACMempty\else%
       \hspace{0.09in}\href{http://portal.acm.org/ft_gateway.cfm?id=\@ACMarticleDOIsuffix&type=pdf}{\icon@PDF}\fi%
    \else
       \hspace{0.09in}\href{\@ACMurlPDF}{\icon@PDF}\fi%
     \ifx\@ACMurlVideo\@ACMempty\else%
       \hspace{0.09in}\href{\@ACMurlVideo}{\icon@Video}\fi%
     \ifx\@ACMurlWebData\@ACMempty\else%
       \hspace{0.09in}\href{\@ACMurlWebData}{\icon@WebData}\fi%
     \ifx\@ACMurlWeb\@ACMempty\else%
       \hspace{0.09in}\href{\@ACMurlWeb}{\icon@Web}\fi%
     \ifx\@ACMurlData\@ACMempty\else%
       \hspace{0.09in}\href{\@ACMurlData}{\icon@Data}\fi%
     \ifx\@ACMurlCode\@ACMempty\else%
       \hspace{0.09in}\href{\@ACMurlCode}{\icon@Code}\fi%
    \fi%
  }%
  \fi
}

% Package hyperref should come late in class file to properly apply to any modified commands.
\RequirePackage{ifpdf}
\ifpdf
  \usepackage[pdftex,breaklinks,hyperfootnotes=false,bookmarks=false]{hyperref}
\else
  \usepackage[hyperfootnotes=false,bookmarks=false]{hyperref}
\fi
\hypersetup{
  colorlinks=true,linkcolor=black,citecolor=black,urlcolor=black,breaklinks=true
}

% remove newlines ('\\') from \title when assigning \pdftitle
\pdfstringdefDisableCommands{%
  % \let~\textasciitilde
  % \def\url{\pdfstringdefWarn\url}%
  % \let\textcolor\@gobble
  \let\\\@gobble
}

\AtBeginDocument{%
  \if@ACMreview%
  \hypersetup{
    pdfcreator={LaTeX acmsiggraph.cls},
    pdfkeywords={\@ACMkeywords},
    pdftitle={\@ACMpdftitle},
  }%
  \else%
  \hypersetup{
    pdfcreator={LaTeX acmsiggraph.cls},
    pdfkeywords={\@ACMkeywords},
    pdftitle={\@ACMpdftitle},
    pdfauthor={\@ACMpdfauthor},
  }%
  \fi
}

%% end of file
\endinput
