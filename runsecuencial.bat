@echo off

for /L %%j in (1,1,10) do (
    for %%i in (600 800 1000 2000 4000 8000 10000) do (
        matricesecuencial.exe %%i >> times2.doc
    )
)

echo Finalizado.
pause