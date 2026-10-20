@echo off
rem One PFGW instance per work_* directory (each needs its own pfgw.ini / pfgw.log).
rem PRPs found are appended to work_*\pfgw.log.  Edit PFGW if pfgw64.exe is not on PATH.
set PFGW=pfgw64.exe
cd /d "%~dp0"
for /d %%D in (work_*) do (
  start "pfgw %%D" /low /d "%%D" %PFGW% -f0 -l"pfgw_out.txt" candidates.txt
)
echo started; check work_*\pfgw.log for PRP results
