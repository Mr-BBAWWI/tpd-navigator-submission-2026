$ErrorActionPreference = 'Stop'
$taskPidFile = Join-Path $PSScriptRoot '.localdata\workbench-20260928\server.pid'
if (-not (Test-Path -LiteralPath $taskPidFile)) { Write-Host 'No launcher process record found.'; exit 0 }
$taskPid = [int](Get-Content -LiteralPath $taskPidFile -Raw)
$taskProcess = Get-CimInstance Win32_Process -Filter "ProcessId=$taskPid"
if (-not $taskProcess) { Write-Host 'The server is already stopped.'; exit 0 }
$taskEntry = Join-Path $PSScriptRoot 'run_workbench.py'
if (-not $taskProcess.CommandLine.Contains($taskEntry)) { throw 'The recorded process belongs to another program. It was not stopped.' }
Stop-Process -Id $taskPid
Write-Host 'TPD Navigator stopped. Saved results remain. In-flight work will be marked interrupted on the next launch.'
