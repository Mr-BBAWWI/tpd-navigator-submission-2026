param([switch]$NoBrowser)
$ErrorActionPreference = 'Stop'
$taskRoot = $PSScriptRoot
$python = Join-Path $taskRoot '.venv\Scripts\python.exe'
$entry = Join-Path $taskRoot 'run_workbench.py'
$taskData = Join-Path $taskRoot '.localdata\workbench-20260928'
$url = 'http://127.0.0.1:8770'
if (-not (Test-Path -LiteralPath $python)) { throw 'Python environment is missing. See active/docs/WORKBENCH_START.md.' }
if (-not (Test-Path -LiteralPath (Join-Path $taskData 'index.sqlite3'))) { throw 'The case store is missing. Import the internal replay package first.' }
$healthy = $false
try { $reply=Invoke-RestMethod "$url/api/health" -TimeoutSec 2; $healthy=$reply.workbench_version -eq '20260930.design.1' } catch {}
if (-not $healthy) {
    $listener=Get-NetTCPConnection -LocalPort 8770 -State Listen -ErrorAction SilentlyContinue
    if ($listener) { throw 'Port 8770 is in use by another app. Do not stop it automatically.' }
    $arguments=@('-X','utf8','-B',('"'+$entry+'"'),'--port','8770')
    $process=Start-Process -FilePath $python -ArgumentList $arguments -WorkingDirectory $taskRoot -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $taskData 'server.stdout.log') -RedirectStandardError (Join-Path $taskData 'server.stderr.log')
    for ($i=0; $i -lt 40; $i++) {
        Start-Sleep -Milliseconds 500
        try {$reply=Invoke-RestMethod "$url/api/health" -TimeoutSec 2; if($reply.workbench_version -eq '20260930.design.1'){$healthy=$true;break}}catch{}
        if($process.HasExited){break}
    }
    if(-not $healthy){throw 'The server did not start. See server.stderr.log in the case store.'}
    $reply.server_pid | Set-Content -LiteralPath (Join-Path $taskData 'server.pid')
}
if (-not $NoBrowser) { Start-Process $url }
Write-Host 'TPD Navigator is ready. Close the browser at any time; saved results remain available.'

