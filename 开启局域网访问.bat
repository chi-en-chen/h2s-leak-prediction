@echo off
chcp 936 >nul
cd /d "%~dp0"
echo ==================================================
echo   开启局域网访问设置
echo --------------------------------------------------
echo   本脚本将自动放行防火墙 5000 端口
echo   以便局域网内其他电脑访问本平台
echo ==================================================
net session >nul 2>&1
if errorlevel 1 (
    echo.
    echo   未检测到管理员权限!
    echo   请右键本文件，选择"以管理员身份运行"
    echo.
    pause
    exit /b
)
netsh advfirewall firewall delete rule name="H2S平台局域网访问" >nul 2>&1
netsh advfirewall firewall add rule name="H2S平台局域网访问" dir=in action=allow protocol=TCP localport=5000
echo.
echo   防火墙已放行 TCP 5000 端口
echo.
echo   本机局域网地址(把其中一个发给同事访问):
setlocal enabledelayedexpansion
for /f "tokens=2 delims=:" %%a in ('ipconfig ^| findstr /i "IPv4"') do (
    set "ip=%%a"
    echo   http://!ip: =!:5000
)
echo.
echo   之后双击"启动网页平台.bat"启动平台即可
echo   同事需与本机连接同一个局域网(路由器/WiFi)
pause