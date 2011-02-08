" Vim syntax extension file
" Language:     ini file
" Maintainer:   Malcolm Augat
" Last Change:  September 30, 2010
" Version:      0.1

syn match iniShortComment   +#.*+
syn region iniLongComment   start=+#/+ end=+/#+ keepend

hi link iniLongComment Comment
hi link iniShortComment Comment

syn keyword Note    NOTE contained containedin=iniShortComment,iniLongComment
syn keyword Malcolm MALCOLM contained containedin=iniShortComment,iniLongComment
syn keyword Todo    TODO contained containedin=iniShortComment,iniLongComment
