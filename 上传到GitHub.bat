@echo off
chcp 936 >nul
cd /d "%~dp0"
echo ==================================================
echo   上传到 GitHub 一键脚本
echo --------------------------------------------------
echo   步骤: 初始化仓库 - 提交文件 - 推送到 GitHub
echo.
echo   开始前请先完成:
echo     1. 注册 GitHub 账号: https://github.com
echo     2. 在 GitHub 网页上新建一个空仓库
echo        (不要勾选自动生成 README)
echo     3. 复制新建仓库的地址
echo ==================================================
git --version >nul 2>&1
if errorlevel 1 (
    echo.
    echo   [错误] 未检测到 git，请先安装:
    echo   https://git-scm.com/download/win
    echo   安装后重新双击本脚本
    pause
    exit /b
)
echo.
echo [1/4] 初始化本地 git 仓库...
if not exist .git (
    git init
    if errorlevel 1 (
        echo   初始化失败
        pause
        exit /b
    )
) else (
    echo   已存在 git 仓库，跳过初始化
)
echo.
echo [2/4] 设置提交身份(首次需要)...
git config user.name >nul 2>&1
if errorlevel 1 set /p GHNAME=  请输入你的 GitHub 用户名: 
if defined GHNAME set "GHNAME=%GHNAME: =%"
if defined GHNAME git config user.name "%GHNAME%"
git config user.email >nul 2>&1
if errorlevel 1 set /p GHEMAIL=  请输入你的 GitHub 邮箱: 
if defined GHEMAIL set "GHEMAIL=%GHEMAIL: =%"
if defined GHEMAIL git config user.email "%GHEMAIL%"
echo.
echo [3/4] 提交全部文件...
git add .
git commit -m "init: H2S leak prediction platform based on EPA SLAB model"
if errorlevel 1 (
    git commit --allow-empty -m "init: H2S leak prediction platform based on EPA SLAB model"
)
git branch -M main
echo.
echo [4/4] 关联并推送到 GitHub...
set "REPOURL="
set /p REPOURL=  请粘贴你的 GitHub 仓库地址(形如 https://github.com/用户名/仓库名.git): 
if defined REPOURL set "REPOURL=%REPOURL: =%"
if not defined REPOURL (
    echo   未输入地址，跳过推送。可稍后重新双击本脚本
    pause
    exit /b
)
git remote | findstr /i "^origin$" >nul
if errorlevel 1 (
    git remote add origin "%REPOURL%"
) else (
    git remote set-url origin "%REPOURL%"
)
git push -u origin main
if errorlevel 1 (
    echo.
    echo   推送失败。常见原因:
    echo     - 仓库地址填错(必须以 .git 结尾)
    echo     - 未登录 GitHub(首次推送会弹出登录窗口，按提示登录)
    echo     - 仓库里已有同名文件(新建仓库时不要勾选 README/.gitignore)
) else (
    echo.
    echo   上传成功! 打开 https://github.com 即可查看你的仓库
)
pause