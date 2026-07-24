' redock_monitor_launch.vbs
' Hidden, fire-and-forget launcher for the redock monitor cron.
' Window style 0 = no console window (nothing for the user to accidentally close,
' which previously killed the cron mid-run and left Task Scheduler zombied at
' 267009 'still running', blocking subsequent checkpoints).
' bWaitOnReturn = False  -> returns immediately, so the Scheduled Task completes
' instantly (exit 0) and can never appear 'still running' / pile up. The monitor
' itself self-limits via its internal 1200s rebuild timeout.
CreateObject("WScript.Shell").Run "wsl -u owner -e bash -c ""cd /mnt/c/Personal/tickdock && python3 scripts/redock_monitor.py >> logs/redock_monitor_cron.log 2>&1""", 0, False
