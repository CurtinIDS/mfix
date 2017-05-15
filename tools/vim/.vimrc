" Syntax highlighting for mfix.dat
au BufRead,BufNewFile mfix.dat set filetype=mfix
au! Syntax mfix source ~/mfix/tools/vim/mfix.vim
syntax on

" Auto-completion for mfix keywords
set complete=k
set dictionary+=~/mfix/tools/vim/mfixdictionary.txt

" Switch from an unsaved buffer without saving it first. Allows to keep an undo
" history for multiple files. Vim will complain if you try to quit without
" saving, and swap files will keep you safe if your computer crashes.
set hidden
" Better command-line completion

set wildmenu

" Show partial commands in the last line of the screen
set showcmd
" Highlight searches (use <C-L> to temporarily turn off highlighting; see the
" mapping of <C-L> below)
set hlsearch

" Use case insensitive search, except when using capital letters
set ignorecase
set smartcase

" Allow backspacing over autoindent, line breaks and start of insert action
set backspace=indent,eol,start

" When opening a new line and no filetype-specific indenting is enabled, keep
" the same indent as the line you're currently on. Useful for READMEs, etc.
set autoindent

" Display the cursor position on the last line of the screen or in the status
" line of a window
set ruler

" Always display the status line, even if only one window is displayed
set laststatus=2

" Instead of failing a command because of unsaved changes, instead raise a
" dialogue asking if you wish to save changed files.
set confirm

" Use visual bell instead of beeping when doing something wrong
set visualbell

" And reset the terminal code for the visual bell.  If visualbell is set, and
" this line is also included, vim will neither flash nor beep.  If visualbell
" is unset, this does nothing.
set t_vb=

" Set the command window height to 2 lines, to avoid many cases of having to
" "press <Enter> to continue"
set cmdheight=2

" Display line numbers on the left
set number
set numberwidth=5

" Cursor options
set cursorcolumn
set cursorline
set showmatch
set expandtab
set shiftwidth=3
set softtabstop=3


