@ECHO OFF
SETLOCAL
PUSHD %~dp0\..

IF "%1"=="" GOTO all
IF /I "%1"=="html" GOTO html
IF /I "%1"=="pdf" GOTO pdf
IF /I "%1"=="all" GOTO all

ECHO Usage: docs\make.bat [html^|pdf^|all]
POPD
EXIT /B 2

:html
hatch run docs:html
GOTO end

:pdf
hatch run docs:pdf
GOTO end

:all
hatch run docs:build

:end
SET RESULT=%ERRORLEVEL%
POPD
EXIT /B %RESULT%
