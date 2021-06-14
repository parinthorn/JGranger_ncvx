@ECHO OFF
for %f in (".\*.eps") do (epstopdf "%f")
PAUSE