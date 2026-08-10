@echo off
cd /d "%~dp0"
echo ==================================================
echo   油井硫化氢(H2S)泄漏扩散预测平台
echo   基于美国 SLAB 重气扩散模型
echo --------------------------------------------------
echo   正在启动计算服务，请稍候...
echo   本窗口请勿关闭，关闭本窗口即停止服务
echo ==================================================
if exist "D:\PYTHONJIESHIQI\python.exe" (
    set "PYEXE=D:\PYTHONJIESHIQI\python.exe"
) else (
    set "PYEXE=python"
)
set "LOGFILE=%~dp0web_platform\server.log"
start /b "" "%PYEXE%" -X utf8 "%~dp0web_platform\app.py" >"%LOGFILE%" 2>&1
echo   正在等待服务就绪...
for /l %%i in (1,1,30) do (
    powershell -NoProfile -Command "try{(New-Object Net.Sockets.TcpClient).Connect('127.0.0.1',5000);exit 0}catch{exit 1}" >nul 2>&1
    if not errorlevel 1 goto ready
    timeout /t 1 /nobreak >nul
)
echo   服务启动失败，请查看 web_platform\server.log
pause
exit /b
:ready
echo.
if exist "%~dp0web_platform\访问地址.txt" type "%~dp0web_platform\访问地址.txt"
echo.
start "" http://127.0.0.1:5000
echo   本机浏览器已自动打开
echo   关闭本窗口即停止服务
pause